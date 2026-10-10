"""Optional global integration and CellTypist benchmark for a merged cell object."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .config import read_json
from .harmonize import (
    _annotate_celltypist,
    _clear_neighbor_graph,
    _compute_neighbors_and_umap,
    _log_step,
    _normalize_nullable_strings_for_h5ad,
    _predict_celltypist,
)
from .logging_utils import configure_logging
from .merge import _match_csr_index_dtypes


def _thin_benchmark_artifact(adata, *, methods: list[dict]) -> object:
    """Keep labels, embeddings, and method graphs without expression data."""
    import anndata as ad

    embedding_keys = [method["embedding_key"] for method in methods]
    representation_keys = [
        method["representation_key"]
        for method in methods
        if method.get("representation_key") is not None
    ]
    graph_keys = [key for method in methods for key in method["graph_keys"].values()]
    missing_embeddings = set(embedding_keys + representation_keys) - set(adata.obsm)
    missing_graphs = set(graph_keys) - set(adata.obsp)
    if missing_embeddings or missing_graphs:
        raise ValueError(
            "Integration benchmark is missing embeddings or neighbor graphs: "
            f"embeddings={sorted(missing_embeddings)}, graphs={sorted(missing_graphs)}"
        )
    result = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32),
        obs=adata.obs.copy(),
    )
    for key in dict.fromkeys(embedding_keys + representation_keys):
        result.obsm[key] = adata.obsm[key].copy()
    for key in graph_keys:
        result.obsp[key] = adata.obsp[key].copy()
    result.uns["integration_benchmark"] = {
        "artifact": "labels_embeddings_and_graphs",
        # The JSON sidecar retains the structured per-method record.
        "method_names": [method["name"] for method in methods],
        "source_matrix_retained": False,
    }
    _normalize_nullable_strings_for_h5ad(result)
    return result


def _compute_scvi_representation(
    adata, method: dict, *, random_seed: int, exclude_vdj_regex: str
) -> tuple[np.ndarray, dict]:
    """Train scVI on a sparse count-based HVG subset and return its latent space."""
    import scanpy as sc
    import scvi

    batch_key = method.get("batch_key")
    hvg_batch_key = method.get("hvg_batch_key", batch_key)
    if batch_key is not None and batch_key not in adata.obs:
        raise ValueError(f"scVI batch key {batch_key!r} is absent")
    if hvg_batch_key is not None and hvg_batch_key not in adata.obs:
        raise ValueError(f"scVI HVG batch key {hvg_batch_key!r} is absent")
    if not sparse.issparse(adata.X):
        raise ValueError("scVI integration expects the merged raw count matrix to be sparse")
    started = time.perf_counter()
    sc.experimental.pp.highly_variable_genes(
        adata,
        n_top_genes=min(method["n_top_genes"], adata.n_vars),
        batch_key=hvg_batch_key,
        chunksize=method.get("hvg_chunksize", 1000),
    )
    selected = adata.var["highly_variable"].to_numpy()
    selected &= ~adata.var_names.to_series().str.contains(
        exclude_vdj_regex, regex=True
    ).to_numpy()
    if not selected.any():
        raise ValueError("scVI integration selected no non-V(D)J highly variable genes")
    scvi_adata = adata[:, selected].copy()
    scvi_adata.X = scvi_adata.X.astype(np.float32, copy=False)
    scvi.settings.seed = random_seed
    if batch_key is None:
        scvi.model.SCVI.setup_anndata(scvi_adata)
    else:
        scvi.model.SCVI.setup_anndata(scvi_adata, batch_key=batch_key)
    model = scvi.model.SCVI(scvi_adata, n_latent=method["n_latent"])
    model.train(
        max_epochs=method["max_epochs"],
        batch_size=method.get("batch_size", 512),
        accelerator="cpu",
        devices=1,
        early_stopping=True,
    )
    latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
    del model, scvi_adata
    elapsed = time.perf_counter() - started
    return latent, {
        "batch_key": batch_key,
        "hvg_batch_key": hvg_batch_key,
        "hvg_genes": int(selected.sum()),
        "n_latent": int(latent.shape[1]),
        "max_epochs": method["max_epochs"],
        "integration": "scVI batch correction" if batch_key is not None else "scVI latent without batch correction",
        "training_seconds": round(elapsed, 2),
    }


def _compute_scanorama_representation(adata, method: dict) -> None:
    """Run Scanorama on PCA scores, passing a lightweight batch-ordered object."""
    import anndata as ad
    import scanpy.external as sce

    batch_key = method["batch_key"]
    values = adata.obs[batch_key].astype(str).to_numpy()
    order = np.argsort(values, kind="stable")
    ordered_obs = adata.obs.iloc[order].copy()
    ordered_obs[batch_key] = pd.Categorical(ordered_obs[batch_key].astype(str))
    work = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32),
        obs=ordered_obs,
    )
    work.obsm[method["basis"]] = adata.obsm[method["basis"]][order].copy()
    sce.pp.scanorama_integrate(
        work,
        key=batch_key,
        basis=method["basis"],
        adjusted_basis=method["adjusted_basis"],
        batch_size=method.get("batch_size", 5000),
    )
    inverse_order = np.argsort(order)
    adata.obsm[method["adjusted_basis"]] = work.obsm[method["adjusted_basis"]][
        inverse_order
    ].astype(np.float32, copy=False)


def _retain_current_method(adata, name: str, label_prefix: str, representation_key: str) -> dict:
    """Save the current UMAP and graph under stable method-specific keys."""
    embedding_key = f"X_umap_{name}"
    graph_keys = {key: f"{name}_{key}" for key in ("distances", "connectivities")}
    adata.obsm[embedding_key] = adata.obsm.pop("X_umap")
    for key, output_key in graph_keys.items():
        adata.obsp[output_key] = adata.obsp[key].copy()
    return {
        "name": name,
        "label_column": f"{label_prefix}aifi_l2_majority",
        "embedding_key": embedding_key,
        "representation_key": representation_key,
        "graph_keys": graph_keys,
    }


def run_integration_benchmark(
    input_path: Path,
    output_path: Path,
    report_path: Path,
    pipeline: dict,
    root: Path,
) -> dict:
    """Run configured diagnostic graph/label variants outside the core merge.

    Each method's graph is shared by its CellTypist majority voting and UMAP,
    then retained with its UMAP and low-dimensional representation for later
    comparisons. scVI is fitted from the loaded raw count matrix before it is
    normalized for the PCA-based methods; only its selected HVG subset is copied.
    """
    import scanpy as sc

    spec = pipeline["integration_benchmark"]
    methods = spec["methods"]
    adata = sc.read_h5ad(input_path)
    adata.X = _match_csr_index_dtypes(adata.X)
    if "aifi_l2_majority" not in adata.obs:
        raise ValueError(f"{input_path}: per-study AIFI L2 labels are absent")
    adata.obs["aifi_l2_study_majority"] = adata.obs["aifi_l2_majority"].astype(str)

    method_reports: list[dict] = []
    annotation = {"method": "celltypist", "levels": [spec["annotation_level"]]}
    scvi_latents = {}
    for scvi_method in (method for method in methods if method["method"] == "scvi"):
        latent, details = _compute_scvi_representation(
            adata,
            scvi_method,
            random_seed=pipeline["processing"]["random_seed"],
            exclude_vdj_regex=pipeline["processing"]["exclude_vdj_regex"],
        )
        scvi_latents[scvi_method["name"]] = (latent, details)

    for method in methods:
        name = method["name"]
        label_prefix = f"benchmark_{name}_"
        if method["method"] == "harmony":
            embedding = _annotate_celltypist(
                adata,
                annotation,
                pipeline,
                root,
                "integration_benchmark",
                embedding_batch_key=method["batch_key"],
                embedding_genes=adata.var_names,
                harmony_basis=method["adjusted_basis"],
                label_prefix=label_prefix,
                restore_counts=False,
                compute_umap=True,
            )
            method_report = _retain_current_method(
                adata, name, label_prefix, method["adjusted_basis"]
            )
            method_reports.append({**method_report, **embedding})
        elif method["method"] == "none":
            _clear_neighbor_graph(adata)
            _compute_neighbors_and_umap(
                adata,
                pipeline,
                "integration_benchmark",
                use_rep=method["basis"],
                compute_umap=True,
            )
            _predict_celltypist(
                adata,
                annotation,
                pipeline,
                root,
                "integration_benchmark",
                label_prefix=label_prefix,
            )
            method_reports.append({
                **_retain_current_method(adata, name, label_prefix, method["basis"]),
                "basis": method["basis"],
                "integration": None,
            })
        elif method["method"] == "scanorama":
            _compute_scanorama_representation(adata, method)
            _clear_neighbor_graph(adata)
            _compute_neighbors_and_umap(
                adata,
                pipeline,
                "integration_benchmark",
                use_rep=method["adjusted_basis"],
                compute_umap=True,
            )
            _predict_celltypist(
                adata, annotation, pipeline, root, "integration_benchmark",
                label_prefix=label_prefix,
            )
            method_reports.append({
                **_retain_current_method(adata, name, label_prefix, method["adjusted_basis"]),
                "neighbors_use_rep": method["adjusted_basis"],
                "integration": {
                    "method": "Scanorama",
                    "batch_key": method["batch_key"],
                    "basis": method["basis"],
                },
            })
        elif method["method"] == "bbknn":
            import bbknn
            import scanpy as sc

            _clear_neighbor_graph(adata)
            n_pcs = min(method.get("n_pcs", pipeline["processing"]["neighbors"]["n_pcs"]),
                        adata.obsm[method["basis"]].shape[1])
            started = time.perf_counter()
            bbknn.bbknn(
                adata,
                batch_key=method["batch_key"],
                use_rep=method["basis"],
                n_pcs=n_pcs,
                neighbors_within_batch=method["neighbors_within_batch"],
            )
            _log_step(
                "integration_benchmark", "bbknn", started,
                batch_key=method["batch_key"], n_pcs=n_pcs,
                neighbors_within_batch=method["neighbors_within_batch"],
            )
            started = time.perf_counter()
            sc.tl.umap(adata, random_state=pipeline["processing"]["random_seed"])
            _log_step("integration_benchmark", "umap", started, use_rep=method["basis"])
            _predict_celltypist(
                adata, annotation, pipeline, root, "integration_benchmark",
                label_prefix=label_prefix,
            )
            method_reports.append({
                **_retain_current_method(adata, name, label_prefix, method["basis"]),
                "neighbors_use_rep": method["basis"],
                "integration": {
                    "method": "BBKNN",
                    "batch_key": method["batch_key"],
                    "neighbors_within_batch": method["neighbors_within_batch"],
                },
            })
        elif method["method"] == "scvi":
            if name not in scvi_latents:
                raise RuntimeError(f"scVI representation for {name!r} was not computed")
            scvi_latent, scvi_details = scvi_latents[name]
            _clear_neighbor_graph(adata)
            adata.obsm[method["adjusted_basis"]] = scvi_latent
            _compute_neighbors_and_umap(
                adata,
                pipeline,
                "integration_benchmark",
                use_rep=method["adjusted_basis"],
                compute_umap=True,
            )
            _predict_celltypist(
                adata, annotation, pipeline, root, "integration_benchmark",
                label_prefix=label_prefix,
            )
            method_reports.append({
                **_retain_current_method(adata, name, label_prefix, method["adjusted_basis"]),
                "neighbors_use_rep": method["adjusted_basis"],
                "integration": scvi_details,
            })
            del scvi_latents[name]
        else:  # Validation rejects this before a workflow is submitted.
            raise ValueError(f"Unsupported integration benchmark method: {method['method']!r}")

    benchmark = _thin_benchmark_artifact(adata, methods=method_reports)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    benchmark.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    report = {
        "kind": "integration_benchmark",
        "status": "ok",
        "input": str(input_path),
        "output": str(output_path),
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "methods": method_reports,
        "artifact": "labels_embeddings_and_graphs",
        "method_graphs_reused_for_celltypist_and_umap": True,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark global PBMC integration variants")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--aifi-l2-model", type=Path, required=True)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()
    configure_logging(args.log_level)
    root = args.project_root.resolve()
    pipeline = read_json(root / args.config)
    pipeline["models"]["aifi_l2"] = str(args.aifi_l2_model.resolve())
    report = run_integration_benchmark(
        args.input, args.output, args.report_output, pipeline, root,
    )
    print(json.dumps({"kind": report["kind"], "status": report["status"], "output": str(args.output)}))

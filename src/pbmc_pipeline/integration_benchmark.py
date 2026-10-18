"""Optional global integration and CellTypist benchmark for a merged cell object."""
from __future__ import annotations

import argparse
import json
import logging
import sys
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

LOGGER = logging.getLogger("pbmc_pipeline.integration_benchmark")


def _thin_benchmark_artifact(
    adata, *, methods: list[dict], sampling: dict | None = None,
) -> object:
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
        "sampling": sampling or {"strategy": "none", "n_cells": int(adata.n_obs)},
    }
    _normalize_nullable_strings_for_h5ad(result)
    return result


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
    comparisons.
    """
    import scanpy as sc

    spec = pipeline["integration_benchmark"]
    methods = spec["methods"]
    LOGGER.info("benchmark=read_input_started path=%s", input_path)
    source = sc.read_h5ad(input_path, backed="r")
    try:
        source_n_cells = int(source.n_obs)
        max_cells = spec.get("max_cells")
        random_seed = int(pipeline["processing"]["random_seed"])
        if max_cells is not None and source_n_cells > max_cells:
            sample_indices = np.sort(
                np.random.default_rng(random_seed).choice(
                    source_n_cells, size=max_cells, replace=False,
                )
            )
            adata = source[sample_indices].to_memory()
            sampling = {
                "strategy": "uniform_random_without_replacement",
                "source_n_cells": source_n_cells,
                "n_cells": int(adata.n_obs),
                "fraction": float(adata.n_obs / source_n_cells),
                "seed": random_seed,
            }
        else:
            adata = source.to_memory()
            sampling = {
                "strategy": "none",
                "source_n_cells": source_n_cells,
                "n_cells": int(adata.n_obs),
                "fraction": 1.0,
                "seed": random_seed,
            }
    finally:
        source.file.close()
    adata.X = _match_csr_index_dtypes(adata.X)
    if "aifi_l2_majority" not in adata.obs:
        raise ValueError(f"{input_path}: per-study AIFI L2 labels are absent")
    adata.obs["aifi_l2_study_majority"] = adata.obs["aifi_l2_majority"].astype(str)
    LOGGER.info(
        "benchmark=read_input_complete source_cells=%d sampled_cells=%d "
        "sampling=%s seed=%d genes=%d studies=%d methods=%d sparse=%s",
        source_n_cells,
        adata.n_obs,
        sampling["strategy"],
        random_seed,
        adata.n_vars,
        adata.obs["study"].nunique() if "study" in adata.obs else 0,
        len(methods),
        sparse.issparse(adata.X),
    )

    method_reports: list[dict] = []
    annotation = {"method": "celltypist", "levels": [spec["annotation_level"]]}
    for method in methods:
        name = method["name"]
        method_started = time.perf_counter()
        LOGGER.info(
            "benchmark=method_started method=%s type=%s cells=%d genes=%d",
            name, method["method"], adata.n_obs, adata.n_vars,
        )
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
        else:  # Validation rejects this before a workflow is submitted.
            raise ValueError(f"Unsupported integration benchmark method: {method['method']!r}")
        LOGGER.info(
            "benchmark=method_complete method=%s duration_seconds=%.2f",
            name, time.perf_counter() - method_started,
        )

    LOGGER.info("benchmark=thin_artifact_started methods=%d", len(method_reports))
    benchmark = _thin_benchmark_artifact(adata, methods=method_reports, sampling=sampling)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("benchmark=write_artifact_started path=%s", output_path)
    benchmark.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    LOGGER.info("benchmark=write_artifact_complete path=%s", output_path)
    report = {
        "kind": "integration_benchmark",
        "status": "ok",
        "input": str(input_path),
        "output": str(output_path),
        "n_cells": int(adata.n_obs),
        "sampling": sampling,
        "n_genes_input": int(adata.n_vars),
        "methods": method_reports,
        "artifact": "labels_embeddings_and_graphs",
        "method_graphs_reused_for_celltypist_and_umap": True,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info("benchmark=complete report=%s artifact=%s", report_path, output_path)
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
    configure_logging(args.log_level, stream=sys.stdout)
    root = args.project_root.resolve()
    pipeline = read_json(root / args.config)
    pipeline["models"]["aifi_l2"] = str(args.aifi_l2_model.resolve())
    report = run_integration_benchmark(
        args.input, args.output, args.report_output, pipeline, root,
    )
    print(json.dumps({"kind": report["kind"], "status": report["status"], "output": str(args.output)}))

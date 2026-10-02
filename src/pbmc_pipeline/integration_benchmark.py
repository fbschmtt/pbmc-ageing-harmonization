"""Optional global integration and CellTypist benchmark for a merged cell object."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from .config import read_json
from .harmonize import (
    _annotate_celltypist,
    _clear_neighbor_graph,
    _compute_neighbors_and_umap,
    _normalize_nullable_strings_for_h5ad,
    _predict_celltypist,
)
from .logging_utils import configure_logging
from .merge import _match_csr_index_dtypes


def _thin_benchmark_artifact(adata, *, methods: list[dict]) -> object:
    """Keep labels and coordinates without duplicating the raw count matrix."""
    import anndata as ad

    if "X_umap" not in adata.obsm:
        raise ValueError("Integration benchmark did not produce a UMAP embedding")
    result = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32),
        obs=adata.obs.copy(),
    )
    result.obsm["X_umap"] = adata.obsm["X_umap"].copy()
    result.uns["integration_benchmark"] = {
        "artifact": "labels_and_coordinates_only",
        # The JSON sidecar retains the structured per-method record.
        "method_names": [method["name"] for method in methods],
        "source_matrix_retained": False,
    }
    _normalize_nullable_strings_for_h5ad(result)
    return result


def run_integration_benchmark(
    input_path: Path,
    output_path: Path,
    report_path: Path,
    pipeline: dict,
    root: Path,
) -> dict:
    """Run configured diagnostic graph/label variants outside the core merge.

    The first Harmony pass constructs the graph used for both its CellTypist
    majority voting and UMAP. The unintegrated comparator needs its own PCA
    graph, but no Harmony graph is rebuilt.
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
            method_reports.append({"name": name, "label_column": f"{label_prefix}aifi_l2_majority", **embedding})
        elif method["method"] == "none":
            _clear_neighbor_graph(adata)
            _compute_neighbors_and_umap(
                adata,
                pipeline,
                "integration_benchmark",
                use_rep=method["basis"],
                compute_umap=False,
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
                "name": name,
                "label_column": f"{label_prefix}aifi_l2_majority",
                "basis": method["basis"],
                "integration": None,
            })
        else:  # Validation rejects this before a workflow is submitted.
            raise ValueError(f"Unsupported integration benchmark method: {method['method']!r}")

    _clear_neighbor_graph(adata)
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
        "artifact": "labels_and_coordinates_only",
        "harmony_graph_reused_for_umap": True,
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

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
    _compute_neighbors_and_umap,
    _harmony_integrate,
    _log_step,
    _normalize_nullable_strings_for_h5ad,
    _predict_celltypist,
)
from .logging_utils import configure_logging
from .merge import _match_csr_index_dtypes

LOGGER = logging.getLogger("pbmc_pipeline.integration_benchmark")


def _thin_benchmark_artifact(
    adata, *, methods: list[dict], sampling: dict | None = None,
    include_graphs: bool = False,
) -> object:
    """Keep labels and UMAPs, optionally retaining graphs in method checkpoints."""
    import anndata as ad

    embedding_keys = [method["embedding_key"] for method in methods]
    graph_keys = [key for method in methods for key in method.get("graph_keys", {}).values()]
    missing_embeddings = set(embedding_keys) - set(adata.obsm)
    missing_graphs = set(graph_keys) - set(adata.obsp)
    if missing_embeddings or (include_graphs and missing_graphs):
        raise ValueError(
            "Integration benchmark is missing method outputs: "
            f"embeddings={sorted(missing_embeddings)}, graphs={sorted(missing_graphs)}"
        )
    result = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 0), dtype=np.float32),
        obs=adata.obs.copy(),
    )
    for key in dict.fromkeys(embedding_keys):
        result.obsm[key] = adata.obsm[key].copy()
    if include_graphs:
        for key in graph_keys:
            result.obsp[key] = adata.obsp[key].copy()
    result.uns["integration_benchmark"] = {
        "artifact": "method_checkpoint" if include_graphs else "labels_and_umaps",
        # The JSON sidecar retains the structured per-method record.
        "method_names": [method["name"] for method in methods],
        "source_matrix_retained": False,
        "neighbor_graphs_retained": include_graphs,
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
    """Keep each method's graph available for checkpoint and aggregation."""
    embedding_key = f"X_umap_{name}"
    adata.obsm[embedding_key] = adata.obsm.pop("X_umap")
    graph_keys = {key: f"{name}_{key}" for key in ("distances", "connectivities")}
    for key, output_key in graph_keys.items():
        adata.obsp[output_key] = adata.obsp.pop(key)
    return {
        "name": name,
        "label_column": f"{label_prefix}aifi_l2_majority",
        "embedding_key": embedding_key,
        "neighbors_use_rep": representation_key,
        "graph_keys": graph_keys,
    }


def _prepare_method_input(input_path: Path, pipeline: dict):
    """Load the same deterministic cell sample independently for each method."""
    import scanpy as sc

    spec = pipeline["integration_benchmark"]
    random_seed = int(pipeline["processing"]["random_seed"])
    LOGGER.info("benchmark=read_input_started path=%s", input_path)
    source = sc.read_h5ad(input_path, backed="r")
    try:
        source_n_cells = int(source.n_obs)
        max_cells = spec.get("max_cells")
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
        "benchmark=read_input_complete source_cells=%d sampled_cells=%d sampling=%s "
        "seed=%d genes=%d studies=%d sparse=%s",
        source_n_cells, adata.n_obs, sampling["strategy"], random_seed,
        adata.n_vars, adata.obs["study"].nunique() if "study" in adata.obs else 0,
        sparse.issparse(adata.X),
    )
    return adata, sampling


def _prepare_normalized_pca(adata, pipeline: dict) -> None:
    """Build the common normalized PCA input used by every adapter."""
    import scanpy as sc

    started = time.perf_counter()
    sc.pp.normalize_total(adata, target_sum=pipeline["processing"]["target_sum"])
    sc.pp.log1p(adata)
    _log_step("integration_benchmark", "normalize_log1p", started)
    # _compute_embedding owns the shared feature selection and PCA path. Use
    # Harmony's helper with no batch key and skip UMAP; adapters add their own
    # correction/graph after this common representation exists.
    from .harmonize import _compute_embedding

    _compute_embedding(
        adata, pipeline, "integration_benchmark", compute_umap=False,
    )


def _run_harmony_adapter(adata, method: dict, pipeline: dict, root: Path) -> dict:
    adjusted_basis = method["adjusted_basis"]
    integration = _harmony_integrate(
        adata, method["batch_key"], basis="X_pca",
        adjusted_basis=adjusted_basis,
        random_state=pipeline["processing"]["random_seed"],
    )
    _compute_neighbors_and_umap(
        adata, pipeline, "integration_benchmark", use_rep=adjusted_basis,
        compute_umap=True,
    )
    label_prefix = f"benchmark_{method['name']}_"
    _predict_celltypist(
        adata, {"method": "celltypist", "levels": [pipeline["integration_benchmark"]["annotation_level"]]},
        pipeline, root, "integration_benchmark", label_prefix=label_prefix,
    )
    return {
        **_retain_current_method(adata, method["name"], label_prefix, adjusted_basis),
        **_embedding_metadata(adata, integration=integration, basis=adjusted_basis),
    }


def _embedding_metadata(adata, *, integration: dict | None, basis: str) -> dict:
    return {
        "basis": basis,
        "neighbors_use_rep": basis,
        "integration": integration,
    }


def _run_unintegrated_adapter(adata, method: dict, pipeline: dict, root: Path) -> dict:
    basis = method["basis"]
    _compute_neighbors_and_umap(
        adata, pipeline, "integration_benchmark", use_rep=basis, compute_umap=True,
    )
    label_prefix = f"benchmark_{method['name']}_"
    _predict_celltypist(
        adata, {"method": "celltypist", "levels": [pipeline["integration_benchmark"]["annotation_level"]]},
        pipeline, root, "integration_benchmark", label_prefix=label_prefix,
    )
    return {
        **_retain_current_method(adata, method["name"], label_prefix, basis),
        **_embedding_metadata(adata, integration=None, basis=basis),
    }


def _run_scanorama_adapter(adata, method: dict, pipeline: dict, root: Path) -> dict:
    _compute_scanorama_representation(adata, method)
    basis = method["adjusted_basis"]
    _compute_neighbors_and_umap(
        adata, pipeline, "integration_benchmark", use_rep=basis, compute_umap=True,
    )
    label_prefix = f"benchmark_{method['name']}_"
    _predict_celltypist(
        adata, {"method": "celltypist", "levels": [pipeline["integration_benchmark"]["annotation_level"]]},
        pipeline, root, "integration_benchmark", label_prefix=label_prefix,
    )
    return {
        **_retain_current_method(adata, method["name"], label_prefix, basis),
        **_embedding_metadata(
            adata,
            integration={"method": "Scanorama", "batch_key": method["batch_key"], "basis": method["basis"]},
            basis=basis,
        ),
    }


def _run_bbknn_adapter(adata, method: dict, pipeline: dict, root: Path) -> dict:
    import bbknn
    import scanpy as sc

    basis = method["basis"]
    n_pcs = min(method.get("n_pcs", pipeline["processing"]["neighbors"]["n_pcs"]),
                adata.obsm[basis].shape[1])
    started = time.perf_counter()
    bbknn.bbknn(
        adata, batch_key=method["batch_key"], use_rep=basis, n_pcs=n_pcs,
        neighbors_within_batch=method["neighbors_within_batch"],
    )
    _log_step(
        "integration_benchmark", "bbknn", started,
        batch_key=method["batch_key"], n_pcs=n_pcs,
        neighbors_within_batch=method["neighbors_within_batch"],
    )
    started = time.perf_counter()
    sc.tl.umap(adata, random_state=pipeline["processing"]["random_seed"])
    _log_step("integration_benchmark", "umap", started, use_rep=basis)
    label_prefix = f"benchmark_{method['name']}_"
    _predict_celltypist(
        adata, {"method": "celltypist", "levels": [pipeline["integration_benchmark"]["annotation_level"]]},
        pipeline, root, "integration_benchmark", label_prefix=label_prefix,
    )
    return {
        **_retain_current_method(adata, method["name"], label_prefix, basis),
        **_embedding_metadata(
            adata,
            integration={
                "method": "BBKNN", "batch_key": method["batch_key"],
                "neighbors_within_batch": method["neighbors_within_batch"],
            },
            basis=basis,
        ),
    }


# Add a new method by implementing an adapter with this signature and
# registering its method type here. Configuration controls its parameters.
METHOD_ADAPTERS = {
    "harmony": _run_harmony_adapter,
    "none": _run_unintegrated_adapter,
    "scanorama": _run_scanorama_adapter,
    "bbknn": _run_bbknn_adapter,
}


def run_integration_method(
    input_path: Path, output_path: Path, status_path: Path,
    method: dict, pipeline: dict, root: Path,
) -> dict:
    """Run one configured adapter and persist its result/status independently."""
    name = method["name"]
    started = time.perf_counter()
    status = {"name": name, "method": method["method"], "status": "failed"}
    try:
        adata, sampling = _prepare_method_input(input_path, pipeline)
        _prepare_normalized_pca(adata, pipeline)
        adapter = METHOD_ADAPTERS[method["method"]]
        method_report = adapter(adata, method, pipeline, root)
        n_genes_input = int(adata.n_vars)
        method_report.update({
            "name": name,
            "method": method["method"],
            "sampling": sampling,
            "n_cells": int(adata.n_obs),
            "n_genes_input": n_genes_input,
        })
        # Drop expression data before constructing/writing the checkpoint. The
        # per-method graph and UMAP are all the aggregator needs from this point.
        adata.X = sparse.csr_matrix(adata.shape, dtype=np.float32)
        adata.layers.clear()
        adata.raw = None
        umap_key = method_report["embedding_key"]
        umap = adata.obsm[umap_key].copy()
        adata.obsm.clear()
        adata.obsm[umap_key] = umap
        adata.varm.clear()
        adata.varp.clear()
        benchmark = _thin_benchmark_artifact(
            adata, methods=[method_report], sampling=sampling, include_graphs=True,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        benchmark.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
        status.update({"status": "ok", "method_report": method_report, "artifact": output_path.name})
        LOGGER.info(
            "benchmark=method_complete method=%s duration_seconds=%.2f",
            name, time.perf_counter() - started,
        )
    except Exception as error:  # Save ordinary method errors and continue the workflow.
        import traceback

        status.update({"error": f"{type(error).__name__}: {error}", "traceback": traceback.format_exc()})
        LOGGER.exception("benchmark=method_failed method=%s", name)
    status["duration_seconds"] = round(time.perf_counter() - started, 2)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    return status


def aggregate_integration_methods(
    input_path: Path, status_paths: list[Path], artifact_paths: list[Path], output_path: Path,
    report_path: Path, pipeline: dict,
) -> dict:
    """Combine successful method artifacts and retain failures in the report."""
    import anndata as ad
    import scanpy as sc

    configured = pipeline["integration_benchmark"]["methods"]
    artifacts_by_name = {
        path.name.removeprefix("method_").removesuffix(".h5ad"): path
        for path in artifact_paths
    }
    results_by_name: dict[str, tuple[dict, Path | None]] = {}
    for status_path in status_paths:
        status = json.loads(status_path.read_text())
        results_by_name[status["name"]] = (status, artifacts_by_name.get(status["name"]))

    methods: list[dict] = []
    failed_methods: list[dict] = []
    successful_artifacts: list[tuple[dict, Path]] = []
    for method in configured:
        name = method["name"]
        status, artifact_path = results_by_name.get(name, ({"name": name, "status": "missing"}, None))
        if status.get("status") == "ok" and artifact_path is not None:
            method_report = status["method_report"]
            methods.append(method_report)
            successful_artifacts.append((method_report, artifact_path))
        else:
            failed_methods.append({
                "name": name,
                "method": method["method"],
                "status": status.get("status", "missing"),
                "error": status.get("error", "Method task did not produce a result status."),
            })

    source = sc.read_h5ad(input_path, backed="r")
    source_n_cells = int(source.n_obs)
    source_n_genes = int(source.n_vars)
    source.file.close()
    n_cells = 0
    sampling = {
        "strategy": "not_available", "source_n_cells": source_n_cells,
        "n_cells": 0, "fraction": 0.0,
        "seed": pipeline["processing"]["random_seed"],
    }
    benchmark = None
    graph_comparisons = []
    baseline_name = next(
        (method["name"] for method in configured if method["name"] == "unintegrated_pca"),
        None,
    )
    baseline_graph = None
    if baseline_name and baseline_name in results_by_name:
        baseline_status, baseline_path = results_by_name[baseline_name]
        if baseline_status.get("status") == "ok" and baseline_path is not None:
            baseline_artifact = ad.read_h5ad(baseline_path, backed="r")
            try:
                baseline_graph = _upper_sparse_graph(
                    baseline_artifact.obsp[f"{baseline_name}_connectivities"]
                )
            finally:
                baseline_artifact.file.close()

    for method_report, artifact_path in successful_artifacts:
        artifact = ad.read_h5ad(artifact_path, backed="r")
        try:
            if benchmark is None:
                n_cells = int(artifact.n_obs)
                sampling = artifact.uns["integration_benchmark"]["sampling"]
                benchmark = ad.AnnData(
                    X=sparse.csr_matrix((artifact.n_obs, 0), dtype=np.float32),
                    obs=artifact.obs.copy(),
                )
            elif not benchmark.obs_names.equals(artifact.obs_names):
                raise ValueError(f"Method {method_report['name']} used a different cell sample")
            label_column = method_report["label_column"]
            if label_column not in benchmark.obs:
                benchmark.obs[label_column] = artifact.obs[label_column]
            embedding_key = method_report["embedding_key"]
            if embedding_key in benchmark.obsm:
                raise ValueError(f"Duplicate method embedding key: {embedding_key}")
            # Read only the configured UMAP. Older cached method checkpoints
            # may also contain PCA arrays under shared keys such as X_pca.
            benchmark.obsm[embedding_key] = artifact.obsm[embedding_key].copy()

            if baseline_graph is not None and method_report["name"] != baseline_name:
                graph = _upper_sparse_graph(
                    artifact.obsp[f"{method_report['name']}_connectivities"]
                )
                shared_edges = baseline_graph.multiply(graph).nnz
                union_edges = baseline_graph.nnz + graph.nnz - shared_edges
                graph_comparisons.append({
                    "graphs": f"{baseline_name} vs {method_report['name']}",
                    "shared_undirected_edges": int(shared_edges),
                    "union_undirected_edges": int(union_edges),
                    "edge_jaccard": shared_edges / union_edges if union_edges else 1.0,
                })
                del graph
        finally:
            artifact.file.close()

    del baseline_graph
    if benchmark is not None:
        benchmark.uns["integration_benchmark"] = {
            "artifact": "labels_and_umaps",
            "method_names": [method["name"] for method in methods],
            "failed_methods": [item["name"] for item in failed_methods],
            "source_matrix_retained": False,
            "neighbor_graphs_retained": False,
            "sampling": sampling,
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        benchmark.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    report = {
        "kind": "integration_benchmark",
        "status": "ok" if not failed_methods else ("partial" if methods else "failed"),
        "input": str(input_path),
        "output": str(output_path) if methods else None,
        "n_cells": n_cells,
        "source_n_cells": source_n_cells,
        "sampling": sampling,
        "n_genes_input": source_n_genes,
        "methods": methods,
        "failed_methods": failed_methods,
        "graph_comparisons": graph_comparisons,
        "artifact": "labels_and_umaps",
        "method_graphs_reused_for_celltypist_and_umap": True,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info(
        "benchmark=aggregate_complete succeeded=%d failed=%d artifact=%s",
        len(methods), len(failed_methods), output_path if methods else "none",
    )
    return report


def _upper_sparse_graph(graph):
    """Return a binary CSR upper triangle for memory-bounded edge comparisons."""
    from scipy import sparse as scipy_sparse

    cleaned = graph.tocsr(copy=True)
    cleaned.setdiag(0)
    cleaned.eliminate_zeros()
    upper = scipy_sparse.triu(cleaned, k=1, format="csr")
    upper.data[:] = 1
    return upper


def main() -> None:
    parser = argparse.ArgumentParser(description="Run or aggregate global PBMC integration benchmarks")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    commands = parser.add_subparsers(dest="command", required=True)
    method_parser = commands.add_parser("method", help="Run one configured integration method")
    method_parser.add_argument("--method-name", required=True)
    method_parser.add_argument("--aifi-l2-model", type=Path, required=True)
    method_parser.add_argument("--output", type=Path, required=True)
    method_parser.add_argument("--status-output", type=Path, required=True)
    aggregate_parser = commands.add_parser("aggregate", help="Combine completed method outputs")
    aggregate_parser.add_argument("--method-statuses", type=Path, nargs="*", default=[])
    aggregate_parser.add_argument("--method-artifacts", type=Path, nargs="*", default=[])
    aggregate_parser.add_argument("--output", type=Path, required=True)
    aggregate_parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    configure_logging(args.log_level, stream=sys.stdout)
    root = args.project_root.resolve()
    pipeline = read_json(root / args.config)
    if args.command == "method":
        try:
            method = next(
                item for item in pipeline["integration_benchmark"]["methods"]
                if item["name"] == args.method_name
            )
        except StopIteration as error:
            raise SystemExit(f"Unknown integration benchmark method: {args.method_name}") from error
        pipeline["models"]["aifi_l2"] = str(args.aifi_l2_model.resolve())
        status = run_integration_method(
            args.input, args.output, args.status_output, method, pipeline, root,
        )
        print(json.dumps(status))
    else:
        report = aggregate_integration_methods(
            args.input, args.method_statuses, args.method_artifacts,
            args.output, args.report_output, pipeline,
        )
        print(json.dumps({"kind": report["kind"], "status": report["status"], "output": report["output"]}))

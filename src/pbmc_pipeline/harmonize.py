from __future__ import annotations

import json
import logging
import platform
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from .config import config_digest
from .metadata import apply_joins, build_homogeneous_obs
from .validation import validate_counts, validate_output

LOGGER = logging.getLogger(__name__)


def _log_step(study_id: str, step: str, started: float, **details: object) -> None:
    """Log a completed processing step with elapsed time and compact context."""
    context = " ".join(f"{key}={value}" for key, value in details.items())
    LOGGER.info(
        "study=%s step=%s duration_seconds=%.2f%s",
        study_id,
        step,
        time.perf_counter() - started,
        f" {context}" if context else "",
    )


def harmonize_study(
    root: Path,
    study_id: str,
    study: dict,
    pipeline: dict,
    schema: dict,
    *,
    test: bool = False,
    validate_only: bool = False,
    input_path: Path | None = None,
    output_path: Path | None = None,
    report_path: Path | None = None,
) -> dict:
    import scanpy as sc

    total_started = time.perf_counter()
    if input_path is None:
        input_path = (root / pipeline["test_input_root"] / study["test_input"] if test
                      else root / study["input"])
    else:
        input_path = input_path.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input for {study_id} does not exist: {input_path}")

    warnings: list[str] = []
    LOGGER.info("study=%s step=start input=%s validate_only=%s", study_id, input_path, validate_only)
    started = time.perf_counter()
    adata = sc.read_h5ad(input_path)
    _log_step(study_id, "read_input", started, cells=adata.n_obs, genes=adata.n_vars)
    if study["counts_source"] == "raw":
        started = time.perf_counter()
        if adata.raw is None:
            raise ValueError(f"{study_id}: configured raw counts are absent")
        original_obs = adata.obs.copy()
        adata = adata.raw.to_adata()
        adata.obs = original_obs.loc[adata.obs_names].copy()
        _log_step(study_id, "select_counts", started, source="raw")
    elif study["counts_source"] != "X":
        raise ValueError(f"Unsupported counts_source: {study['counts_source']}")
    else:
        LOGGER.info("study=%s step=select_counts source=X", study_id)

    started = time.perf_counter()
    if hasattr(adata.X, "sum_duplicates"):
        adata.X.sum_duplicates()
        adata.X.sort_indices()
    validate_counts(adata)
    _log_step(study_id, "validate_input", started)
    input_cells, input_genes = adata.shape
    started = time.perf_counter()
    adata, dropped_features = _repair_features(adata, study.get("features", {}))
    _log_step(
        study_id, "repair_features", started, dropped_duplicate_features=dropped_features,
        cells=adata.n_obs, genes=adata.n_vars,
    )
    started = time.perf_counter()
    joined_obs = apply_joins(adata.obs, study.get("joins", []), root, warnings)
    adata.obs = build_homogeneous_obs(joined_obs, study_id, study["metadata"])
    _log_step(study_id, "build_metadata", started, warnings=len(warnings))

    annotation = study["annotation"]
    if annotation["method"] == "celltypist" and not validate_only:
        _annotate_celltypist(adata, annotation, pipeline, root, study_id)
    elif annotation["method"] == "celltypist":
        LOGGER.info("study=%s step=annotate_celltypist status=skipped_validate_only", study_id)
        for level in annotation["levels"]:
            adata.obs[f"aifi_{level}_majority"] = "not_run"
    elif annotation["method"] != "existing":
        raise ValueError(f"Unknown annotation method: {annotation['method']}")
    else:
        LOGGER.info("study=%s step=annotate status=existing_labels", study_id)

    started = time.perf_counter()
    validate_counts(adata)
    summary = validate_output(adata, schema)
    _log_step(study_id, "validate_output", started, cells=adata.n_obs, genes=adata.n_vars)
    report = {
        "study": study_id,
        "status": "validated" if validate_only else "complete",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "input_cells": int(input_cells), "input_genes": int(input_genes),
        "dropped_duplicate_features": dropped_features,
        **summary,
        "warnings": warnings,
        "provenance": study.get("provenance", {}),
        "configuration_sha256": config_digest(study, pipeline, schema),
        "versions": {
            "python": platform.python_version(),
            "anndata": version("anndata"),
            "scanpy": version("scanpy"),
        },
    }
    suffix = ".test" if test else ""
    if report_path is None:
        report_path = root / pipeline["report_dir"] / f"{study_id}{suffix}.json"
    started = time.perf_counter()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    _log_step(study_id, "write_run_report", started, path=report_path)
    if not validate_only:
        if output_path is None:
            output_path = root / pipeline["output_dir"] / f"{study_id}{suffix}.h5ad"
        started = time.perf_counter()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        _normalize_nullable_strings_for_h5ad(adata)
        adata.write_h5ad(
            output_path,
            compression=pipeline["processing"]["output_compression"],
        )
        _log_step(study_id, "write_h5ad", started, path=output_path)
    _log_step(study_id, "complete", total_started, status=report["status"])
    return report


def _normalize_nullable_strings_for_h5ad(adata) -> None:
    """Use broadly compatible object strings instead of Pandas nullable strings."""
    for frame in (adata.obs, adata.var):
        if isinstance(frame.index.dtype, pd.StringDtype):
            index = frame.index.astype(object).where(~frame.index.isna(), None)
            frame.index = pd.Index(index, name=frame.index.name, dtype=object)
        for column in frame.columns:
            series = frame[column]
            if isinstance(series.dtype, pd.StringDtype):
                frame[column] = series.astype(object).where(series.notna(), None)


def _repair_features(adata, spec: dict):
    import pandas as pd

    names = (adata.var[spec["source_column"]].astype(str)
             if spec.get("source_column") else pd.Series(adata.var_names, index=adata.var_names))
    if "split" in spec:
        rule = spec["split"]
        names = names.str.split(rule["separator"]).str[rule["index"]]
    duplicated = names.duplicated(keep=False).to_numpy()
    duplicate_count = int(duplicated.sum())
    policy = spec.get("duplicate_policy", "error")
    if duplicate_count and policy == "drop_all":
        adata = adata[:, ~duplicated].copy()
        names = names.loc[~duplicated]
    elif duplicate_count:
        raise ValueError(f"Found {duplicate_count} ambiguous feature names")
    adata.var_names = names.to_numpy()
    return adata, duplicate_count


def _annotate_celltypist(
    adata, annotation: dict, pipeline: dict, root: Path, study_id: str
) -> None:
    import celltypist
    import scanpy as sc

    raw_counts = adata.X.copy()
    started = time.perf_counter()
    sc.pp.normalize_total(adata, target_sum=pipeline["processing"]["target_sum"])
    sc.pp.log1p(adata)
    _log_step(study_id, "normalize_log1p", started)
    _compute_embedding(adata, pipeline, study_id)
    for level in annotation["levels"]:
        started = time.perf_counter()
        model_path = root / pipeline["models"][f"aifi_{level}"]
        prediction = celltypist.annotate(
            adata, model=celltypist.models.Model.load(str(model_path)), majority_voting=True
        )
        labels = prediction.predicted_labels["majority_voting"]
        if not labels.index.equals(adata.obs_names):
            raise ValueError(f"CellTypist {level} prediction index differs from input cells")
        adata.obs[f"aifi_{level}_majority"] = labels
        _log_step(study_id, "annotate_celltypist", started, level=level, model=model_path.name)
    adata.X = raw_counts


def _compute_embedding(adata, pipeline: dict, study_id: str) -> None:
    import scanpy as sc

    hvg = pipeline["processing"]["hvg"]
    started = time.perf_counter()
    sc.pp.highly_variable_genes(adata, **hvg)
    vdj = adata.var_names.to_series().str.contains(
        pipeline["processing"]["exclude_vdj_regex"], regex=True
    )
    selected = adata.var["highly_variable"].to_numpy() & ~vdj.to_numpy()
    embedding = adata[:, selected].copy()
    _log_step(study_id, "select_embedding_features", started, selected_features=int(selected.sum()))
    started = time.perf_counter()
    sc.pp.scale(embedding)
    _log_step(study_id, "scale_embedding", started)
    started = time.perf_counter()
    sc.tl.pca(embedding, svd_solver="arpack", random_state=pipeline["processing"]["random_seed"])
    _log_step(study_id, "pca", started, components=embedding.obsm["X_pca"].shape[1])
    neighbors = pipeline["processing"]["neighbors"]
    started = time.perf_counter()
    sc.pp.neighbors(embedding, **neighbors, random_state=pipeline["processing"]["random_seed"])
    _log_step(study_id, "neighbors", started, n_neighbors=neighbors["n_neighbors"])
    started = time.perf_counter()
    sc.tl.umap(embedding, random_state=pipeline["processing"]["random_seed"])
    _log_step(study_id, "umap", started)
    for key in ("uns", "obsm", "obsp"):
        setattr(adata, key, getattr(embedding, key))

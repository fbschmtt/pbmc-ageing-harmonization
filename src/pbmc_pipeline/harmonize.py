from __future__ import annotations

import json
import logging
import platform
import re
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

from .config import config_digest
from .preparation import normalize_metadata_frame, read_prepared_cells
from .source import read_study_input
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
    prepared_obs_path: Path | None = None,
    output_path: Path | None = None,
    report_path: Path | None = None,
) -> dict:
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
    adata = read_study_input(root, input_path, study)
    _log_step(study_id, "read_input", started, cells=adata.n_obs, genes=adata.n_vars)
    if study["counts_source"] == "raw":
        started = time.perf_counter()
        if adata.raw is None:
            raise ValueError(f"{study_id}: configured raw counts are absent")
        original_obs = adata.obs.copy()
        adata = adata.raw.to_adata()
        adata.obs = original_obs.loc[adata.obs_names].copy()
        _log_step(study_id, "select_counts", started, source="raw")
    elif study["counts_source"] == "X":
        LOGGER.info("study=%s step=select_counts source=X", study_id)
    elif study["counts_source"].startswith("layer:"):
        layer = study["counts_source"].removeprefix("layer:")
        if not layer or layer not in adata.layers:
            raise ValueError(f"{study_id}: configured count layer is absent: {layer!r}")
        started = time.perf_counter()
        adata.X = adata.layers[layer].copy()
        _log_step(study_id, "select_counts", started, source=f"layer:{layer}")
    else:
        raise ValueError(f"Unsupported counts_source: {study['counts_source']}")

    # The harmonized artifact is a count matrix, not a source-object archive.
    # A source H5AD can retain normalized X/raw/layer matrices alongside the
    # configured counts; remove every alternate expression representation once
    # the selected counts are in X so none can be mistaken for analysis input.
    discarded_layers = list(adata.layers)
    had_raw = adata.raw is not None
    adata.raw = None
    for layer in discarded_layers:
        del adata.layers[layer]
    LOGGER.info(
        "study=%s step=discard_alternate_expression raw=%s layers=%s",
        study_id,
        had_raw,
        ",".join(discarded_layers) if discarded_layers else "none",
    )

    started = time.perf_counter()
    if hasattr(adata.X, "sum_duplicates"):
        adata.X.sum_duplicates()
        adata.X.sort_indices()
    validate_counts(adata)
    _log_step(study_id, "validate_input", started)
    source_input_cells = adata.n_obs
    if prepared_obs_path is None:
        raise ValueError("A prepared per-cell metadata artifact is required")
    started = time.perf_counter()
    prepared = read_prepared_cells(
        prepared_obs_path,
        adata.obs_names,
        schema,
        allow_cell_subset=study["preparation"].get("allow_cell_subset", False),
    )
    if len(prepared) < source_input_cells:
        adata = adata[prepared.index].copy()
        LOGGER.info(
            "study=%s step=filter_input_cells retained=%s excluded=%s",
            study_id,
            adata.n_obs,
            source_input_cells - adata.n_obs,
        )
    adata.obs = prepared
    _log_step(study_id, "attach_prepared_metadata", started, cells=adata.n_obs, warnings=len(warnings))
    input_cells, input_genes = adata.shape
    started = time.perf_counter()
    feature_spec = study.get("features", {})
    feature_name_metrics = _feature_name_metrics(adata, feature_spec)
    adata, aggregated_features = _repair_features(adata, feature_spec)
    _log_step(
        study_id, "repair_features", started, aggregated_duplicate_features=aggregated_features,
        cells=adata.n_obs, genes=adata.n_vars,
        feature_symbols=feature_name_metrics["n_feature_symbols"],
        feature_ids=feature_name_metrics["n_feature_ids"],
    )
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
        "source_input_cells": int(source_input_cells),
        "input_cells": int(input_cells), "input_genes": int(input_genes),
        "excluded_input_cells": int(source_input_cells - input_cells),
        "aggregated_duplicate_features": aggregated_features,
        "feature_name_metrics": feature_name_metrics,
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
        adata.uns["pipeline_provenance"] = {
            "artifact_kind": "harmonized",
            "study": study_id,
            "configuration_sha256": report["configuration_sha256"],
            "input": str(input_path),
            "pipeline_version": "0.1.0",
        }
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


def _normalize_missing_metadata(adata, schema: dict) -> None:
    """Apply the output contract's missing-value representation before validation."""
    adata.obs = normalize_metadata_frame(adata.obs, schema)


ENSEMBL_FEATURE_RE = re.compile(r"^ENS[A-Z]*\d+(?:\.\d+)?$")


def _feature_names(adata, spec: dict) -> pd.Series:
    """Return the configured feature labels before duplicate aggregation."""
    names = (adata.var[spec["source_column"]].astype("string")
             if spec.get("source_column") else pd.Series(adata.var_names, index=adata.var_names, dtype="string"))
    if "split" in spec:
        rule = spec["split"]
        names = names.str.split(rule["separator"]).str[rule["index"]]
    return names


def _feature_name_metrics(adata, spec: dict) -> dict[str, int | str]:
    """Summarize symbol-like and identifier-like labels supplied by the source."""
    names = _feature_names(adata, spec)
    nonempty = names.notna() & names.str.strip().ne("")
    feature_ids = names.str.match(ENSEMBL_FEATURE_RE, na=False)
    duplicated = names[nonempty].duplicated(keep=False)
    duplicate_names = names[nonempty][duplicated].nunique()
    return {
        "source_column": spec.get("source_column", "var_names"),
        "n_features": len(names),
        "n_feature_symbols": int((nonempty & ~feature_ids).sum()),
        "n_feature_ids": int((nonempty & feature_ids).sum()),
        "n_missing_feature_names": int((~nonempty).sum()),
        "n_duplicate_feature_entries": int(duplicated.sum()),
        "n_duplicate_feature_names": int(duplicate_names),
    }


def _repair_features(adata, spec: dict):
    import pandas as pd
    from scipy import sparse

    names = _feature_names(adata, spec).astype(str)
    duplicated = names.duplicated(keep=False).to_numpy()
    duplicate_count = int(duplicated.sum())
    policy = spec.get("duplicate_policy", "error")
    if duplicate_count and policy == "sum":
        codes, unique_names = pd.factorize(names, sort=False)
        membership = sparse.csr_matrix(
            (np.ones(adata.n_vars, dtype=np.int8), (np.arange(adata.n_vars), codes)),
            shape=(adata.n_vars, len(unique_names)),
        )
        counts = (adata.X @ membership if sparse.issparse(adata.X)
                  else sparse.csr_matrix(adata.X) @ membership)
        representatives = ~names.duplicated(keep="first")
        adata = adata[:, representatives].copy()
        adata.X = counts
        names = pd.Series(unique_names, index=adata.var_names)
    elif duplicate_count:
        raise ValueError(f"Found {duplicate_count} ambiguous feature names")
    adata.var_names = names.to_numpy()
    # Feature names now differ from the source gene-ID index; retaining its
    # ``gene_id`` name conflicts with the preserved ``var["gene_id"]`` column
    # when AnnData serializes the object.
    adata.var_names.name = None
    return adata, duplicate_count


def _annotate_celltypist(
    adata,
    annotation: dict,
    pipeline: dict,
    root: Path,
    study_id: str,
    *,
    embedding_batch_key: str | None = None,
    embedding_genes: pd.Index | None = None,
    harmony_basis: str = "X_pca_harmony",
    label_prefix: str = "",
) -> dict[str, object]:
    import scanpy as sc

    raw_counts = adata.X.copy()
    started = time.perf_counter()
    sc.pp.normalize_total(adata, target_sum=pipeline["processing"]["target_sum"])
    sc.pp.log1p(adata)
    _log_step(study_id, "normalize_log1p", started)
    embedding = _compute_embedding(
        adata,
        pipeline,
        study_id,
        batch_key=embedding_batch_key,
        embedding_genes=embedding_genes,
        harmony_basis=harmony_basis,
    )
    _predict_celltypist(adata, annotation, pipeline, root, study_id, label_prefix=label_prefix)
    adata.X = raw_counts
    return embedding


def _predict_celltypist(
    adata,
    annotation: dict,
    pipeline: dict,
    root: Path,
    study_id: str,
    *,
    label_suffix: str = "",
    label_prefix: str = "",
) -> None:
    """Run CellTypist majority voting using the graph currently on ``adata``."""
    import celltypist

    for level in annotation["levels"]:
        started = time.perf_counter()
        model_path = root / pipeline["models"][f"aifi_{level}"]
        prediction = celltypist.annotate(
            adata, model=celltypist.models.Model.load(str(model_path)), majority_voting=True
        )
        labels = prediction.predicted_labels["majority_voting"]
        if not labels.index.equals(adata.obs_names):
            raise ValueError(f"CellTypist {level} prediction index differs from input cells")
        adata.obs[f"{label_prefix}aifi_{level}{label_suffix}_majority"] = labels
        _log_step(
            study_id, "annotate_celltypist", started, level=level,
            model=model_path.name, graph=adata.uns["neighbors"]["params"].get("use_rep"),
        )


def _compute_embedding(
    adata,
    pipeline: dict,
    study_id: str,
    *,
    batch_key: str | None = None,
    embedding_genes: pd.Index | None = None,
    harmony_basis: str = "X_pca_harmony",
) -> dict[str, object]:
    """Compute an embedding, optionally correcting its PCs with Harmony.

    When ``embedding_genes`` is supplied, feature selection is deliberately
    confined to that intersection before HVG calculation. Only the selected
    HVGs are copied for scaling/PCA, keeping the full raw count matrix out of
    the embedding path.
    """
    import scanpy as sc

    if embedding_genes is None or adata.var_names.equals(embedding_genes):
        embedding_input = adata
    else:
        missing = embedding_genes.difference(adata.var_names)
        if not missing.empty:
            raise ValueError(f"{study_id}: embedding genes absent from expression matrix")
        embedding_input = adata[:, embedding_genes].copy()
    if batch_key is not None and batch_key not in embedding_input.obs:
        raise ValueError(f"{study_id}: Harmony batch key {batch_key!r} is absent")
    embedding_input_genes = int(embedding_input.n_vars)

    hvg = pipeline["processing"]["hvg"]
    started = time.perf_counter()
    sc.pp.highly_variable_genes(embedding_input, **hvg)
    vdj = embedding_input.var_names.to_series().str.contains(
        pipeline["processing"]["exclude_vdj_regex"], regex=True
    )
    selected = embedding_input.var["highly_variable"].to_numpy() & ~vdj.to_numpy()
    embedding = embedding_input[:, selected].copy()
    if embedding.n_vars == 0:
        raise ValueError(f"{study_id}: no non-VDJ highly variable genes were selected")
    # Use float32 and only reach this point after feature selection to bound the
    # embedding working set.
    embedding.X = embedding.X.astype(np.float32, copy=False)
    _log_step(
        study_id,
        "select_embedding_features",
        started,
        input_features=embedding_input.n_vars,
        selected_features=int(selected.sum()),
    )
    if embedding is not adata and embedding_input is not adata:
        del embedding_input
    started = time.perf_counter()
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="zero-centering a sparse array/matrix densifies it.",
            category=UserWarning,
        )
        sc.pp.scale(embedding)
    _log_step(study_id, "scale_embedding", started)
    started = time.perf_counter()
    sc.tl.pca(embedding, svd_solver="arpack", random_state=pipeline["processing"]["random_seed"])
    _log_step(study_id, "pca", started, components=embedding.obsm["X_pca"].shape[1])
    neighbor_rep = "X_pca"
    if batch_key is not None:
        import scanpy.external as sce

        started = time.perf_counter()
        sce.pp.harmony_integrate(
            embedding,
            batch_key,
            basis="X_pca",
            adjusted_basis=harmony_basis,
            random_state=pipeline["processing"]["random_seed"],
            verbose=False,
        )
        neighbor_rep = harmony_basis
        _log_step(study_id, "harmony_integrate", started, batch_key=batch_key, basis=harmony_basis)
    neighbors = pipeline["processing"]["neighbors"]
    started = time.perf_counter()
    neighbor_args = dict(neighbors)
    effective_n_pcs = min(neighbor_args["n_pcs"], embedding.obsm[neighbor_rep].shape[1])
    neighbor_args["n_pcs"] = effective_n_pcs
    sc.pp.neighbors(
        embedding,
        **neighbor_args,
        use_rep=neighbor_rep,
        random_state=pipeline["processing"]["random_seed"],
    )
    _log_step(
        study_id,
        "neighbors",
        started,
        n_neighbors=neighbors["n_neighbors"], n_pcs=effective_n_pcs,
        use_rep=neighbor_rep,
    )
    started = time.perf_counter()
    sc.tl.umap(embedding, random_state=pipeline["processing"]["random_seed"])
    _log_step(study_id, "umap", started)
    for key in ("uns", "obsm", "obsp"):
        setattr(adata, key, getattr(embedding, key))
    return {
        "basis": "X_umap",
        "input_genes": embedding_input_genes,
        "hvg_genes": int(embedding.n_vars),
        "n_pcs": int(embedding.obsm["X_pca"].shape[1]),
        "neighbors_n_pcs": effective_n_pcs,
        "neighbors_use_rep": neighbor_rep,
        "integration": None if batch_key is None else {
            "method": "harmony", "batch_key": batch_key, "adjusted_basis": harmony_basis,
        },
    }

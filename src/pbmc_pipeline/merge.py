"""Cross-study pseudobulk and optional single-cell merge operations."""
from __future__ import annotations

import json
import logging
import platform
from collections.abc import Iterable
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .config import config_digest
from .harmonize import _annotate_celltypist, _normalize_nullable_strings_for_h5ad

LOGGER = logging.getLogger(__name__)


def _write_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")


def _matrix_nnz_by_row(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.diff(matrix.tocsr().indptr)
    return np.count_nonzero(matrix, axis=1)


def pseudobulk_study(
    input_path: Path, output_path: Path, report_path: Path, pipeline: dict
) -> dict:
    """Sum raw counts by sample and AIFI-L2, retaining legacy bulk metadata."""
    import anndata as ad
    import scanpy as sc

    spec = pipeline["merge"]["pseudobulk"]
    groupby = list(spec["groupby"])
    adata = sc.read_h5ad(input_path)
    missing = set(groupby) - set(adata.obs)
    if missing:
        raise ValueError(f"{input_path}: pseudobulk grouping columns absent: {sorted(missing)}")
    # Convert categoricals before filling: a missing label must not require a
    # category mutation in source study metadata.
    groups = adata.obs[groupby].astype("string").fillna("not_provided").astype(str)
    group_index = pd.MultiIndex.from_frame(groups)
    codes, unique = pd.factorize(group_index, sort=True)
    membership = sparse.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.int64), (codes, np.arange(adata.n_obs))),
        shape=(len(unique), adata.n_obs),
    )
    counts = membership @ adata.X
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    obs = _legacy_bulk_metadata(adata.obs, codes, len(unique))
    for position, column in enumerate(groupby):
        obs[column] = unique.get_level_values(position).astype(str)
    study_id = str(adata.obs["study"].iloc[0])
    identifier_values = unique.to_list()
    if "study" not in groupby:
        identifier_values = [(study_id, *values) for values in identifier_values]
    obs.index = pd.Index(
        ["::".join(values) for values in identifier_values], name="pseudobulk_id"
    )
    obs["n_cells"] = np.bincount(codes, minlength=len(unique)).astype(int)
    obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    obs["library_size"] = obs["total_counts"]
    obs["detected_genes"] = _matrix_nnz_by_row(counts).astype(int)
    output = ad.AnnData(X=counts, obs=obs, var=adata.var.copy())
    output.uns["pseudobulk"] = {
        "groupby": groupby,
        "metadata_policy": "legacy_homogeneous_or_heterogeneous_during_bulk",
        "source": str(input_path),
    }
    output.uns["pipeline_provenance"] = {
        "artifact_kind": "pseudobulk",
        "study": study_id,
        "configuration_sha256": config_digest(pipeline),
        "input": str(input_path),
        "pipeline_version": "0.1.0",
    }
    _normalize_nullable_strings_for_h5ad(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    report = {
        "kind": "pseudobulk_study", "status": "complete",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "input": str(input_path),
        "study": str(adata.obs["study"].iloc[0]), "groupby": groupby,
        "input_cells": int(adata.n_obs), "input_genes": int(adata.n_vars),
        "n_pseudobulks": int(output.n_obs), "n_genes": int(output.n_vars),
        "cells_per_pseudobulk": _distribution(obs["n_cells"]),
        "configuration_sha256": config_digest(pipeline),
    }
    _write_report(report_path, report)
    return report


def merge_pseudobulks(
    input_paths: Iterable[Path], output_path: Path, report_path: Path, pipeline: dict
) -> dict:
    """Outer-join per-study pseudobulk matrices, treating missing genes as zero."""
    import anndata as ad
    import scanpy as sc

    paths = list(input_paths)
    parts = [sc.read_h5ad(path) for path in paths]
    if not parts:
        raise ValueError("No pseudobulk inputs supplied")
    merged = ad.concat(parts, join="outer", merge="same", index_unique=None, fill_value=0)
    if not merged.obs_names.is_unique:
        raise ValueError("Pseudobulk identifiers are not unique across studies")
    _annotate_gene_availability(merged, parts)
    merged.uns["pipeline_provenance"] = _merge_provenance(
        "pseudobulk_merge", paths, pipeline, gene_join="outer"
    )
    _normalize_nullable_strings_for_h5ad(merged)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    report = _merge_report("pseudobulk_merge", merged, paths, pipeline)
    _write_report(report_path, report)
    return report


def _shared_gene_names(parts: list) -> pd.Index:
    """Return common genes in the first study's order without materializing matrices."""
    shared = parts[0].var_names
    for part in parts[1:]:
        shared = shared.intersection(part.var_names, sort=False)
    return shared


def merge_single_cells(
    input_paths: Iterable[Path], output_path: Path, report_path: Path, pipeline: dict, root: Path
) -> dict:
    """Inner-join genes, embed all cells jointly, and obtain fresh AIFI L2 labels.

    This intentionally loads every input in memory and is only invoked by the
    explicitly enabled Nextflow pathway.
    """
    import anndata as ad
    import scanpy as sc

    paths = list(input_paths)
    parts = [sc.read_h5ad(path) for path in paths]
    if not parts:
        raise ValueError("No single-cell inputs supplied")
    spec = pipeline["merge"]["single_cell"]
    shared_genes = _shared_gene_names(parts)
    if shared_genes.empty:
        raise ValueError("Studies have no shared genes for single-cell merge")
    merged = ad.concat(
        parts, join=spec["gene_join"], merge="same", index_unique="::", label="source_study"
    )
    if merged.n_vars == 0:
        raise ValueError("Studies have no shared genes for single-cell merge")
    merged.obs["aifi_l2_study_majority"] = merged.obs["aifi_l2_majority"].astype(str)
    # The shipped configuration uses an inner join.  Keep provenance available
    # if a caller explicitly opts into an outer join, where fill_value=0 creates
    # synthetic rather than observational zeros.
    if spec["gene_join"] == "outer":
        _annotate_gene_availability(merged, parts)
    integration = spec["integration"]
    # CellTypist detects ``obsp['connectivities']`` and uses it for its
    # over-clustering/majority-voting stage. The helper therefore installs the
    # Harmony-derived graph before CellTypist is invoked.
    embedding = _annotate_celltypist(
        merged,
        {"method": "celltypist", "levels": [spec["annotation_level"]]},
        pipeline,
        root,
        "merged_single_cell",
        embedding_batch_key=integration["batch_key"],
        embedding_genes=shared_genes,
        harmony_basis=integration["adjusted_basis"],
    )
    merged.uns["pipeline_provenance"] = _merge_provenance(
        "single_cell_merge", paths, pipeline, gene_join=spec["gene_join"]
    )
    # The helper restores raw counts after computing the Harmony graph/UMAP.
    merged.uns["embedding"] = embedding
    _normalize_nullable_strings_for_h5ad(merged)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    report = _merge_report("single_cell_merge", merged, paths, pipeline)
    report.update({"embedding": embedding})
    _write_report(report_path, report)
    return report


def _distribution(values: pd.Series) -> dict[str, float]:
    return {key: float(values.quantile(q)) for key, q in {"min": 0, "p50": .5, "p95": .95, "max": 1}.items()}


def _legacy_bulk_metadata(obs: pd.DataFrame, codes: np.ndarray, n_groups: int) -> pd.DataFrame:
    """Retain homogeneous metadata and mark mixed values as in the legacy notebook."""
    values_by_column: dict[str, list[object]] = {}
    for column in obs.columns:
        aggregated: list[object] = []
        series = obs[column]
        for group_code in range(n_groups):
            values = series.iloc[np.flatnonzero(codes == group_code)]
            present = values[values.notna()]
            if present.empty:
                aggregated.append(np.nan)
            elif len(present) == len(values) and present.nunique(dropna=True) == 1:
                aggregated.append(present.iloc[0])
            else:
                aggregated.append("Heterogeneous During Bulk")
        if "Heterogeneous During Bulk" in aggregated:
            values_by_column[column] = [
                value if isinstance(value, str) else str(value) for value in aggregated
            ]
        else:
            values_by_column[column] = aggregated
    return pd.DataFrame(values_by_column)


def _annotate_gene_availability(merged, parts: list) -> None:
    """Expose which outer-join zeros are synthetic for every source study."""
    study_to_part = {str(part.obs["study"].iloc[0]): part for part in parts}
    if len(study_to_part) != len(parts):
        raise ValueError("Pseudobulk inputs must contain one unique study each")
    availability: dict[str, np.ndarray] = {}
    for study_id, part in sorted(study_to_part.items()):
        available = merged.var_names.isin(part.var_names)
        merged.var[f"available_in_{study_id}"] = available
        availability[study_id] = available
    availability_frame = pd.DataFrame(availability, index=merged.var_names)
    merged.var["n_studies_with_gene"] = availability_frame.sum(axis=1).astype(int)
    merged.var["studies_with_gene"] = availability_frame.apply(
        lambda row: ",".join(row.index[row.to_numpy()]), axis=1
    )
    merged.var["synthetic_zero_filled_studies"] = availability_frame.apply(
        lambda row: ",".join(row.index[~row.to_numpy()]), axis=1
    )
    merged.var["has_synthetic_zeros"] = merged.var["n_studies_with_gene"] < len(parts)


def _merge_provenance(kind: str, paths: list[Path], pipeline: dict, *, gene_join: str) -> dict:
    return {
        "artifact_kind": kind,
        "configuration_sha256": config_digest(pipeline),
        "inputs": [str(path) for path in paths],
        "gene_join": gene_join,
        "pipeline_version": "0.1.0",
    }


def _merge_report(kind: str, adata, paths: list[Path], pipeline: dict) -> dict:
    per_study = adata.obs.groupby("study", observed=True).size().sort_index()
    celltype = (adata.obs.groupby(["study", "aifi_l2_majority"], observed=True).size()
                .rename("n").reset_index())
    return {
        "kind": kind, "status": "complete", "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": [str(path) for path in paths], "n_observations": int(adata.n_obs),
        "n_genes": int(adata.n_vars), "observations_by_study": per_study.astype(int).to_dict(),
        "samples_by_study": adata.obs.groupby("study", observed=True)["sample"].nunique().astype(int).to_dict(),
        "subjects_by_study": adata.obs.groupby("study", observed=True)["subject"].nunique().astype(int).to_dict(),
        "aifi_l2_by_study": celltype.to_dict(orient="records"),
        "genes_with_synthetic_zeros": int(adata.var.get("has_synthetic_zeros", pd.Series(False)).sum()),
        "configuration_sha256": config_digest(pipeline),
        "versions": {"python": platform.python_version(), "anndata": version("anndata"), "scanpy": version("scanpy")},
    }

"""Cross-study pseudobulk and optional single-cell merge operations."""
from __future__ import annotations

import json
import logging
import platform
import tempfile
import time
from collections.abc import Iterable
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .config import config_digest
from .harmonize import (
    _annotate_celltypist,
    _normalize_nullable_strings_for_h5ad,
    _predict_celltypist,
)

LOGGER = logging.getLogger(__name__)


def _start_step(study_id: str, step: str, **details: object) -> float:
    """Log the start of a potentially long stage and return its timer."""
    context = " ".join(f"{key}={value}" for key, value in details.items())
    LOGGER.info("study=%s step=%s start%s", study_id, step, f" {context}" if context else "")
    return time.perf_counter()


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


def _write_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")


def _matrix_nnz_by_row(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.diff(matrix.tocsr().indptr)
    return np.count_nonzero(matrix, axis=1)


def _match_csr_index_dtypes(matrix):
    """Make CSR index arrays compatible with Scanpy's Numba normalizer.

    ``concat_on_disk`` can write an ``int32`` ``indices`` array alongside an
    ``int64`` ``indptr`` array. Scanpy requires both arrays to share a dtype
    when normalizing a CSR matrix. Prefer downcasting the much smaller pointer
    array; retain 64-bit pointers when the matrix has too many nonzeros.
    """
    if not sparse.isspmatrix_csr(matrix) or matrix.indices.dtype == matrix.indptr.dtype:
        return matrix
    if int(matrix.indptr[-1]) <= np.iinfo(matrix.indices.dtype).max:
        matrix.indptr = matrix.indptr.astype(matrix.indices.dtype, copy=False)
    else:
        matrix.indices = matrix.indices.astype(matrix.indptr.dtype, copy=False)
    return matrix


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
    # Pseudobulks are small, dense downstream interchange inputs. Keep them
    # uncompressed; all cell-level H5AD artifacts use the configured compression.
    output.write_h5ad(output_path, compression=None)
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
    # See ``pseudobulk_study``: merged pseudobulk is deliberately uncompressed.
    merged.write_h5ad(output_path, compression=None)
    report = _merge_report("pseudobulk_merge", merged, paths, pipeline)
    report["gene_join_accounting"] = _gene_join_count_accounting(
        parts, merged.var_names, "outer",
    )
    _write_report(report_path, report)
    return report


def _shared_gene_names(parts: list) -> pd.Index:
    """Return common genes in the first study's order without materializing matrices."""
    shared = parts[0].var_names
    for part in parts[1:]:
        shared = shared.intersection(part.var_names, sort=False)
    return shared


def _gene_join_count_accounting(parts: list, joined_genes: pd.Index, gene_join: str) -> dict:
    """Count raw UMIs retained or excluded by the merge's gene-name join."""
    study_rows = []
    total_counts = 0
    retained_counts = 0
    for part in parts:
        gene_names = pd.Index(part.var_names.astype(str))
        retained_mask = gene_names.isin(joined_genes)
        gene_counts = np.zeros(part.n_vars, dtype=np.int64)
        for start in range(0, part.n_obs, 2048):
            block = part.X[start:min(start + 2048, part.n_obs)]
            if hasattr(block, "to_memory"):
                block = block.to_memory()
            block_totals = (
                np.asarray(block.sum(axis=0)).ravel()
                if sparse.issparse(block)
                else np.asarray(block).sum(axis=0, dtype=np.int64)
            )
            gene_counts += block_totals.astype(np.int64, copy=False)

        study_total = int(gene_counts.sum(dtype=np.int64))
        study_retained = int(gene_counts[retained_mask].sum(dtype=np.int64))
        discarded_mask = ~retained_mask
        study_discarded = study_total - study_retained
        top_discarded_gene = None
        top_discarded_gene_counts = 0
        if study_discarded:
            discarded_positions = np.flatnonzero(discarded_mask)
            top_position = discarded_positions[np.argmax(gene_counts[discarded_mask])]
            top_discarded_gene = str(gene_names[top_position])
            top_discarded_gene_counts = int(gene_counts[top_position])
        study_rows.append({
            "study": str(part.obs["study"].iloc[0]),
            "input_gene_symbols": int(part.n_vars),
            "retained_gene_symbols": int(retained_mask.sum()),
            "discarded_gene_symbols": int(discarded_mask.sum()),
            "input_counts": study_total,
            "retained_counts": study_retained,
            "discarded_counts": study_discarded,
            "discarded_fraction": study_discarded / study_total if study_total else 0.0,
            "top_discarded_gene": top_discarded_gene,
            "top_discarded_gene_counts": top_discarded_gene_counts,
            "top_gene_fraction_of_discarded_counts": (
                top_discarded_gene_counts / study_discarded if study_discarded else 0.0
            ),
        })
        total_counts += study_total
        retained_counts += study_retained
    discarded_counts = total_counts - retained_counts
    return {
        "gene_join": gene_join,
        "joined_gene_symbols": len(joined_genes),
        "input_counts": total_counts,
        "retained_counts": retained_counts,
        "discarded_counts": discarded_counts,
        "discarded_fraction": discarded_counts / total_counts if total_counts else 0.0,
        "studies": sorted(study_rows, key=lambda row: row["study"]),
    }


def _link_h5ad_for_merge(source: Path, destination: Path) -> None:
    """Create a lightweight H5AD exposing only slots used by merge.

    Single-cell merging uses X, obs, and var. Source embeddings, graphs, and
    unstructured analysis results are discarded and recomputed on the merged
    object. HDF5 external links expose the required source slots without
    copying the large count matrix or making concat_on_disk process embeddings.
    """
    import h5py

    with h5py.File(source, "r") as input_file, h5py.File(destination, "w") as output_file:
        for key, value in input_file.attrs.items():
            output_file.attrs[key] = value
        for key in ("X", "obs", "var"):
            output_file[key] = h5py.ExternalLink(str(source.resolve()), f"/{key}")
        for key in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
            group = output_file.create_group(key)
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"


def merge_single_cells(
    input_paths: Iterable[Path], output_path: Path, report_path: Path, pipeline: dict, root: Path
) -> dict:
    """Inner-join genes and compute experimental merged AIFI-L2 predictions.

    Expression matrices are concatenated on disk to avoid retaining all source
    matrices alongside the merged matrix in memory.
    """
    import anndata as ad
    import scanpy as sc

    paths = list(input_paths)
    if not paths:
        raise ValueError("No single-cell inputs supplied")
    spec = pipeline["merge"]["single_cell"]
    study_id = "merged_single_cell"
    LOGGER.info(
        "study=%s step=start inputs=%d gene_join=%s",
        study_id, len(paths), spec["gene_join"],
    )
    # Read only the small obs/var metadata in backed mode. concat_on_disk
    # streams the large expression matrices into one temporary H5AD, avoiding
    # the previous peak of all source matrices plus a full in-memory concat.
    started = _start_step(study_id, "read_input_metadata", inputs=len(paths))
    parts = [ad.read_h5ad(path, backed="r") for path in paths]
    try:
        shared_genes = _shared_gene_names(parts)
        if shared_genes.empty:
            raise ValueError("Studies have no shared genes for single-cell merge")
        if spec["gene_join"] == "inner":
            joined_genes = shared_genes
        else:
            joined_genes = parts[0].var_names
            for part in parts[1:]:
                joined_genes = joined_genes.union(part.var_names, sort=False)
        accounting_started = _start_step(
            study_id, "gene_join_count_accounting", join=spec["gene_join"]
        )
        gene_join_accounting = _gene_join_count_accounting(
            parts, joined_genes, spec["gene_join"],
        )
        _log_step(study_id, "gene_join_count_accounting", accounting_started)
        _log_step(
            study_id, "read_input_metadata", started,
            inputs=len(paths), shared_genes=len(shared_genes),
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        # Keep this temporary concat alive as the raw-count checkpoint.
        concat_temp = tempfile.TemporaryDirectory(
            prefix=".single-cell-concat-", dir=output_path.parent
        )
        concat_path = Path(concat_temp.name) / "merged.h5ad"
        concat_inputs = []
        csc_inputs = [
            index for index, part in enumerate(parts)
            if getattr(part.X, "format", None) == "csc"
        ]
        # Close the metadata-only handles before creating lean input copies.
        for part in parts:
            part.file.close()
        for index, source_path in enumerate(paths):
            started = _start_step(study_id, "prepare_concat_input", input=source_path.name)
            lean_path = Path(concat_temp.name) / f"input_{index}.h5ad"
            if index not in csc_inputs:
                _link_h5ad_for_merge(source_path, lean_path)
            else:
                # AnnData's on-disk concat requires a common sparse format.
                # Convert one CSC source at a time and omit unused slots.
                conversion_started = _start_step(
                    study_id, "convert_input_to_csr", input=source_path.name
                )
                part = sc.read_h5ad(source_path)
                part.X = part.X.tocsr()
                lean = ad.AnnData(X=part.X, obs=part.obs.copy(), var=part.var.copy())
                lean.write_h5ad(lean_path)
                del lean, part
                _log_step(
                    study_id, "convert_input_to_csr", conversion_started,
                    input=source_path.name,
                )
            concat_inputs.append(lean_path)
            _log_step(study_id, "prepare_concat_input", started, input=source_path.name)
        started = _start_step(
            study_id, "concat_on_disk", inputs=len(paths), gene_join=spec["gene_join"]
        )
        ad.experimental.concat_on_disk(
            concat_inputs,
            concat_path,
            join=spec["gene_join"],
            merge="same",
            index_unique="::",
            label="source_study",
            # The AnnData default can load roughly 400 MB of sparse data at
            # once. Smaller chunks reduce concat's working set on RAM-limited
            # production machines.
            max_loaded_elems=10_000_000,
        )
        _log_step(study_id, "concat_on_disk", started)
        started = _start_step(study_id, "load_merged_matrix", path=concat_path)
        merged = sc.read_h5ad(concat_path)
        merged.X = _match_csr_index_dtypes(merged.X)
        _log_step(
            study_id, "load_merged_matrix", started,
            cells=merged.n_obs, genes=merged.n_vars,
        )
    finally:
        for part in parts:
            part.file.close()
    if merged.n_vars == 0:
        raise ValueError("Studies have no shared genes for single-cell merge")
    # Per-study embeddings are not used: the merged embedding below is
    # recomputed from counts. Drop them before allocating the global versions.
    merged.obsm.clear()
    # The shipped configuration uses an inner join.  Keep provenance available
    # if a caller explicitly opts into an outer join, where fill_value=0 creates
    # synthetic rather than observational zeros.
    if spec["gene_join"] == "outer":
        started = _start_step(study_id, "gene_availability")
        _annotate_gene_availability(merged, parts)
        _log_step(study_id, "gene_availability", started)
    del parts
    started = _start_step(study_id, "qc_metrics", cells=merged.n_obs, genes=merged.n_vars)
    _add_single_cell_qc_metrics(merged)
    _log_step(study_id, "qc_metrics", started)
    integration = spec["integration"]
    # Keep each study's own L2 annotation as the downstream ground truth. The
    # merged prediction is diagnostic only and must never overwrite it.
    merged.obs["aifi_l2_study_majority"] = merged.obs["aifi_l2_majority"].astype(str)
    embedding_genes = merged.var_names if spec["gene_join"] == "inner" else shared_genes
    # CellTypist detects ``obsp['connectivities']`` and uses it for its
    # over-clustering/majority-voting stage. The helper therefore installs the
    # Harmony-derived graph before CellTypist is invoked.
    started = _start_step(study_id, "harmony_celltypist")
    embedding = _annotate_celltypist(
        merged,
        {"method": "celltypist", "levels": [spec["annotation_level"]]},
        pipeline,
        root,
        "merged_single_cell",
        embedding_batch_key=integration["batch_key"],
        embedding_genes=embedding_genes,
        harmony_basis=integration["adjusted_basis"],
        label_prefix="experimental_",
        restore_counts=False,
        compute_umap=False,
    )
    _log_step(study_id, "harmony_celltypist", started)
    # The Harmony-graph prediction is experimental. Repeat the same diagnostic
    # with an unintegrated-PC graph, then restore Harmony's graph/UMAP as the
    # primary embedding stored in the merged artifact.
    _clear_neighbor_graph(merged)
    # The Harmony/CellTypist pass already left X normalized and log1p
    # transformed, so reuse it for the unintegrated-PC prediction.
    _compute_neighbors_and_umap(
        merged, pipeline, "merged_single_cell", use_rep="X_pca", compute_umap=False
    )
    _predict_celltypist(
        merged,
        {"method": "celltypist", "levels": [spec["annotation_level"]]},
        pipeline,
        root,
        "merged_single_cell",
        label_suffix="_unintegrated",
        label_prefix="experimental_",
    )
    # The unintegrated graph has served its prediction pass. Remove it before
    # rebuilding the final Harmony graph to avoid holding both full graphs.
    _clear_neighbor_graph(merged)
    # Release normalized counts before reading the raw checkpoint back. The
    # backed read materializes only X, not the old per-study embeddings.
    started = _start_step(study_id, "restore_raw_counts")
    merged.X = None
    backed_counts = ad.read_h5ad(concat_path, backed="r")
    try:
        backed_matrix = backed_counts.X
        raw_counts = (
            backed_matrix.to_memory()
            if hasattr(backed_matrix, "to_memory")
            else np.asarray(backed_matrix)
        )
    finally:
        backed_counts.file.close()
    merged.X = raw_counts
    del raw_counts, backed_matrix, backed_counts
    _log_step(study_id, "restore_raw_counts", started)
    started = _start_step(study_id, "final_harmony_graph_umap")
    _compute_neighbors_and_umap(
        merged, pipeline, "merged_single_cell", use_rep=integration["adjusted_basis"]
    )
    _log_step(study_id, "final_harmony_graph_umap", started)
    merged.uns["pipeline_provenance"] = _merge_provenance(
        "single_cell_merge", paths, pipeline, gene_join=spec["gene_join"]
    )
    # The raw-count checkpoint has been restored for the published artifact.
    merged.uns["embedding"] = embedding
    _normalize_nullable_strings_for_h5ad(merged)
    started = _start_step(study_id, "write_merged_h5ad", output=output_path)
    merged.write_h5ad(output_path, compression=pipeline["processing"]["output_compression"])
    _log_step(study_id, "write_merged_h5ad", started)
    started = _start_step(study_id, "write_merge_report", output=report_path)
    report = _merge_report("single_cell_merge", merged, paths, pipeline)
    report["gene_join_accounting"] = gene_join_accounting
    report.update({"embedding": embedding})
    _write_report(report_path, report)
    _log_step(study_id, "write_merge_report", started)
    concat_temp.cleanup()
    return report


def _clear_neighbor_graph(adata) -> None:
    """Release graph matrices once their CellTypist/UMAP consumer is done."""
    for key in ("distances", "connectivities"):
        if key in adata.obsp:
            del adata.obsp[key]
    adata.uns.pop("neighbors", None)


def _add_single_cell_qc_metrics(adata) -> None:
    """Add raw-count depth and mitochondrial fraction for merged-cell QC."""
    counts = adata.X
    total_counts = np.asarray(counts.sum(axis=1)).ravel()
    mitochondrial = adata.var_names.str.upper().str.startswith("MT-")
    mito_counts = np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel()
    adata.obs["UMIs_per_cell"] = total_counts
    adata.obs["percent_mito"] = np.divide(
        mito_counts * 100,
        total_counts,
        out=np.zeros(adata.n_obs, dtype=float),
        where=total_counts != 0,
    )


def _compute_neighbors_and_umap(
    adata, pipeline: dict, study_id: str, *, use_rep: str, compute_umap: bool = True
) -> None:
    """Install a graph and UMAP derived from one existing PCA representation."""
    import scanpy as sc

    neighbors = pipeline["processing"]["neighbors"]
    neighbor_args = dict(neighbors)
    effective_n_pcs = min(neighbor_args["n_pcs"], adata.obsm[use_rep].shape[1])
    neighbor_args["n_pcs"] = effective_n_pcs
    started = _start_step(study_id, "neighbors", use_rep=use_rep)
    sc.pp.neighbors(
        adata,
        **neighbor_args,
        use_rep=use_rep,
        random_state=pipeline["processing"]["random_seed"],
    )
    LOGGER.info(
        "study=%s step=neighbors duration_seconds=%.2f n_neighbors=%s n_pcs=%s use_rep=%s",
        study_id, time.perf_counter() - started, neighbors["n_neighbors"], effective_n_pcs, use_rep,
    )
    if compute_umap:
        started = _start_step(study_id, "umap", use_rep=use_rep)
        sc.tl.umap(adata, random_state=pipeline["processing"]["random_seed"])
        LOGGER.info(
            "study=%s step=umap duration_seconds=%.2f use_rep=%s",
            study_id, time.perf_counter() - started, use_rep,
        )


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
    if "n_cells" in adata.obs:
        celltype = (adata.obs.groupby(["study", "aifi_l2_majority"], observed=True)["n_cells"]
                    .sum().rename("n_cells").reset_index())
    else:
        celltype = (adata.obs.groupby(["study", "aifi_l2_majority"], observed=True).size()
                    .rename("n_cells").reset_index())
    return {
        "kind": kind, "status": "complete", "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": [str(path) for path in paths], "n_observations": int(adata.n_obs),
        "n_genes": int(adata.n_vars), "observations_by_study": per_study.astype(int).to_dict(),
        "samples_by_study": adata.obs.groupby("study", observed=True)["sample"].nunique().astype(int).to_dict(),
        "subjects_by_study": adata.obs.groupby("study", observed=True)["subject"].nunique().astype(int).to_dict(),
        "aifi_l2_cells_by_study": celltype.to_dict(orient="records"),
        "genes_with_synthetic_zeros": int(adata.var.get("has_synthetic_zeros", pd.Series(False)).sum()),
        "configuration_sha256": config_digest(pipeline),
        "versions": {"python": platform.python_version(), "anndata": version("anndata"), "scanpy": version("scanpy")},
    }

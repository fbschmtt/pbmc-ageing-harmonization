"""Split a merged single-cell H5AD and compute per-cell-type embeddings."""
from __future__ import annotations

import argparse
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import read_json
from .harmonize import _harmony_integrate, _normalize_nullable_strings_for_h5ad
from .utils import slugify


def add_unintegrated_umap(adata, pipeline: dict) -> str:
    """Add an unintegrated UMAP from the retained native PCA scores.

    A separate neighbor graph key keeps the analysis artifact's Harmony graph
    and its existing integrated UMAP unchanged. The temporary graph is removed
    after coordinates are computed; only the additional embedding is retained.
    """
    import scanpy as sc

    if "X_pca" not in adata.obsm:
        raise ValueError("An unintegrated UMAP requires native PCA scores in obsm['X_pca']")
    if adata.n_obs < 3:
        raise ValueError("An unintegrated UMAP requires at least three cells")

    neighbors_key = "unintegrated_neighbors"
    embedding_key = "X_umap_unintegrated"
    neighbors = pipeline["processing"]["neighbors"]
    n_pcs = min(int(neighbors["n_pcs"]), adata.obsm["X_pca"].shape[1])
    if n_pcs < 2:
        raise ValueError("An unintegrated UMAP requires at least two native PCA components")

    sc.pp.neighbors(
        adata,
        n_neighbors=min(int(neighbors["n_neighbors"]), adata.n_obs - 1),
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=pipeline["processing"]["random_seed"],
        key_added=neighbors_key,
    )
    sc.tl.umap(
        adata,
        random_state=pipeline["processing"]["random_seed"],
        neighbors_key=neighbors_key,
        key_added=embedding_key,
    )

    graph_metadata = adata.uns.pop(neighbors_key)
    for graph_field in ("distances_key", "connectivities_key"):
        graph_key = graph_metadata.get(graph_field)
        if graph_key in adata.obsp:
            del adata.obsp[graph_key]
    adata.uns.pop(embedding_key, None)
    return embedding_key


def split_cell_types(input_path: Path, output_dir: Path, manifest_path: Path, pipeline: dict) -> dict:
    """Write raw-count H5AD splits while keeping the merged matrix on disk."""
    import anndata as ad

    split_by = pipeline["cell_type_analysis"]["split_by"]
    # Backed mode keeps the large sparse expression matrix in the merged file.
    # Only the current type's rows are materialized before writing its split.
    adata = ad.read_h5ad(input_path, backed="r")
    output_dir.mkdir(parents=True, exist_ok=True)
    records, used_slugs = [], set()
    try:
        if split_by not in adata.obs:
            raise ValueError(f"{input_path}: split column {split_by!r} is absent")
        # Keep the original column rather than making a second full-size string
        # array for large merges. Check empty labels against the small unique set.
        labels = adata.obs[split_by]
        if labels.isna().any():
            raise ValueError(f"{input_path}: split column {split_by!r} contains missing labels")
        unique_labels = labels.unique().tolist()
        if any(not str(label).strip() for label in unique_labels):
            raise ValueError(f"{input_path}: split column {split_by!r} contains missing labels")
        for label in sorted(unique_labels, key=str):
            slug = slugify(str(label))
            if slug in used_slugs:
                raise ValueError(f"Cell-type labels produce the same output name: {label!r} -> {slug!r}")
            used_slugs.add(slug)
            selected = (labels == label).to_numpy(dtype=bool, na_value=False, copy=False)
            subset = adata[selected].to_memory()
            del selected
            subset.obsm.clear()
            subset.obsp.clear()
            subset.uns = {
                "pipeline_provenance": {
                    "artifact_kind": "cell_type_split",
                    "input": str(input_path),
                    "split_by": split_by,
                    "cell_type": str(label),
                }
            }
            _normalize_nullable_strings_for_h5ad(subset)
            path = output_dir / f"{slug}.h5ad"
            subset.write_h5ad(path, compression=pipeline["processing"]["output_compression"])
            records.append({
                "cell_type": str(label), "slug": slug,
                "n_cells": int(subset.n_obs), "path": path.name,
            })
            del subset
    finally:
        adata.file.close()
    document = {"kind": "cell_type_split", "input": str(input_path), "split_by": split_by, "cell_types": records}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(document, indent=2) + "\n")
    return document


def analyse_cell_type(input_path: Path, output_dir: Path, pipeline: dict) -> dict:
    """Create a type-specific embedding, Leiden clusters, and descriptive tables."""
    import scanpy as sc

    adata = sc.read_h5ad(input_path)
    spec = pipeline["cell_type_analysis"]
    output_dir.mkdir(parents=True, exist_ok=True)
    cell_types = adata.obs["aifi_l2_majority"].astype(str).unique()
    if len(cell_types) != 1:
        raise ValueError(f"{input_path}: expected exactly one configured L2 type and parent L1")
    parent = spec.get("l2_parent_l1", {}).get(cell_types[0])
    if parent is None:
        raise ValueError(
            f"{input_path}: no configured AIFI L2 → L1 parent mapping for {cell_types[0]!r}"
        )
    # Parent L1 is report metadata derived from the shipped taxonomy. Keep it
    # out of both the raw-count split and the analysis AnnData observations.
    report = {
        "kind": "cell_type_analysis", "status": "complete", "input": str(input_path),
        "cell_type": cell_types[0], "aifi_l1_parent": parent,
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    if adata.n_obs < spec["min_cells"]:
        report.update({"status": "insufficient_cells", "minimum_cells": spec["min_cells"]})
        _normalize_nullable_strings_for_h5ad(adata)
        adata.write_h5ad(output_dir / "analysis.h5ad", compression=pipeline["processing"]["output_compression"])
        (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        return report

    raw_counts = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=pipeline["processing"]["target_sum"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, **pipeline["processing"]["hvg"])
    excluded = adata.var_names.to_series().str.contains(pipeline["processing"]["exclude_vdj_regex"], case=False, regex=True)
    selected = adata.var["highly_variable"].to_numpy() & ~excluded.to_numpy()
    if selected.sum() < 2:
        selected = ~excluded.to_numpy()
    working = adata[:, selected].copy()
    # Native PCA is retained for the exploratory PC--age correlations below.
    # A future analysis may reasonably prefer the Harmony-adjusted PCs there,
    # but that choice should be made deliberately rather than implicitly.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="zero-centering a sparse array/matrix densifies it.",
            category=UserWarning,
        )
        sc.pp.scale(working, max_value=10)
    n_comps = min(pipeline["processing"]["neighbors"]["n_pcs"], working.n_obs - 1, working.n_vars - 1)
    if n_comps < 2:
        raise ValueError(f"{input_path}: too few cells or genes for PCA")
    sc.tl.pca(working, n_comps=n_comps, random_state=pipeline["processing"]["random_seed"])
    integration = spec["integration"]
    batch_key = integration["batch_key"]
    neighbor_rep = "X_pca"
    integration_report = {"method": "harmony", "batch_key": batch_key, "adjusted_basis": integration["adjusted_basis"]}
    if working.obs[batch_key].nunique() > 1:
        # Integrate studies for the local neighbour graph, UMAP, and Leiden
        # clustering. This should be revisited as the downstream analysis
        # matures: correction can also remove biological age-associated signal.
        integration_report.update(_harmony_integrate(
            working,
            batch_key,
            basis="X_pca",
            adjusted_basis=integration["adjusted_basis"],
            random_state=pipeline["processing"]["random_seed"],
        ))
        neighbor_rep = integration["adjusted_basis"]
    else:
        integration_report["status"] = "skipped_single_study"
    sc.pp.neighbors(
        working,
        n_neighbors=min(pipeline["processing"]["neighbors"]["n_neighbors"], working.n_obs - 1),
        n_pcs=min(n_comps, working.obsm[neighbor_rep].shape[1]),
        use_rep=neighbor_rep,
        random_state=pipeline["processing"]["random_seed"],
    )
    sc.tl.umap(working, random_state=pipeline["processing"]["random_seed"])
    sc.tl.leiden(
        working,
        key_added="cluster",
        resolution=spec["clustering"]["resolution"],
        random_state=pipeline["processing"]["random_seed"],
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    # Wilcoxon log fold-changes require non-negative log-normalized data.
    # ``working`` is scaled for PCA, so rank markers on the matching unscaled
    # normalized values to avoid undefined log-fold changes.
    marker_input = adata[:, selected].copy()
    marker_input.obs["cluster"] = working.obs["cluster"].astype(str).to_numpy()
    sc.tl.rank_genes_groups(marker_input, groupby="cluster", method="wilcoxon", use_raw=False)
    markers = _marker_table(sc.get.rank_genes_groups_df(marker_input, group=None), marker_input.obs["cluster"])
    markers.to_csv(output_dir / "markers.tsv", sep="\t", index=False)
    adata.obsm["X_pca"] = working.obsm["X_pca"]
    adata.uns["pca"] = working.uns["pca"].copy()
    # Retain loadings in the full-gene coordinate system so reports can use
    # Scanpy's pca_loadings plot. Excluded genes have zero loading because they
    # were intentionally outside the local PCA input.
    pca_loadings = np.zeros((adata.n_vars, n_comps), dtype=working.varm["PCs"].dtype)
    pca_loadings[selected] = working.varm["PCs"]
    adata.varm["PCs"] = pca_loadings
    adata.var["used_for_local_pca"] = selected
    if neighbor_rep != "X_pca":
        adata.obsm[integration["adjusted_basis"]] = working.obsm[integration["adjusted_basis"]]
    adata.obsm["X_umap"] = working.obsm["X_umap"]
    adata.obsp["distances"] = working.obsp["distances"]
    adata.obsp["connectivities"] = working.obsp["connectivities"]
    adata.uns["neighbors"] = working.uns["neighbors"]
    adata.obs["cluster"] = working.obs["cluster"].astype(str)
    adata.X = raw_counts
    composition = adata.obs.groupby(["cluster", "study"], observed=True).size().rename("n_cells").reset_index()
    composition.to_csv(output_dir / "cluster_study_composition.tsv", sep="\t", index=False)
    pc_age = _pc_age_correlations(adata)
    pc_age.to_csv(output_dir / "pc_age_correlations.tsv", sep="\t", index=False)
    report.update({"n_clusters": int(adata.obs["cluster"].nunique()), "clustering": spec["clustering"], "embedding": {"n_hvg": int(selected.sum()), "n_vdj_genes_excluded": int(excluded.sum()), "n_pcs": n_comps, "neighbors_use_rep": neighbor_rep, "integration": integration_report}})
    _normalize_nullable_strings_for_h5ad(adata)
    adata.write_h5ad(output_dir / "analysis.h5ad", compression=pipeline["processing"]["output_compression"])
    (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def _marker_table(markers: pd.DataFrame, clusters: pd.Series) -> pd.DataFrame:
    """Make Scanpy's single-cluster marker table conform to the multi-cluster schema."""
    if "group" in markers:
        return markers
    unique_clusters = clusters.astype(str).unique()
    if len(unique_clusters) != 1:
        raise ValueError("Scanpy marker table has no group column for multiple clusters")
    return markers.assign(group=unique_clusters[0])


def _pc_age_correlations(adata) -> pd.DataFrame:
    age = pd.to_numeric(adata.obs.get("age"), errors="coerce")
    valid = age.notna().to_numpy()
    rows = []
    for index in range(min(10, adata.obsm["X_pca"].shape[1])):
        if valid.sum() < 3:
            correlation, pvalue = np.nan, np.nan
        else:
            correlation, pvalue = spearmanr(adata.obsm["X_pca"][valid, index], age[valid])
        rows.append({"pc": f"PC{index + 1}", "spearman_r": correlation, "pvalue": pvalue, "n_cells": int(valid.sum())})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Split or analyse per-cell-type PBMC H5AD artifacts")
    subparsers = parser.add_subparsers(dest="command", required=True)
    split = subparsers.add_parser("split")
    split.add_argument("--input", type=Path, required=True)
    split.add_argument("--output-dir", type=Path, required=True)
    split.add_argument("--manifest", type=Path, required=True)
    split.add_argument("--config", type=Path, required=True)
    analyse = subparsers.add_parser("analyse")
    analyse.add_argument("--input", type=Path, required=True)
    analyse.add_argument("--output-dir", type=Path, required=True)
    analyse.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    pipeline = read_json(args.config)
    if args.command == "split":
        split_cell_types(args.input, args.output_dir, args.manifest, pipeline)
    else:
        analyse_cell_type(args.input, args.output_dir, pipeline)


if __name__ == "__main__":
    main()

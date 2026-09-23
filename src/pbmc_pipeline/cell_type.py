"""Split a merged single-cell H5AD and compute per-cell-type embeddings."""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import read_json
from .harmonize import _normalize_nullable_strings_for_h5ad


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "unnamed"


def split_cell_types(input_path: Path, output_dir: Path, manifest_path: Path, pipeline: dict) -> dict:
    """Load the global merge once and write one raw-count H5AD per AIFI L2 type."""
    import scanpy as sc

    split_by = pipeline["cell_type_analysis"]["split_by"]
    adata = sc.read_h5ad(input_path)
    if split_by not in adata.obs:
        raise ValueError(f"{input_path}: split column {split_by!r} is absent")
    labels = adata.obs[split_by].astype("string")
    if labels.isna().any() or (labels.str.strip() == "").any():
        raise ValueError(f"{input_path}: split column {split_by!r} contains missing labels")
    _add_sample_cell_counts(adata, pipeline["cell_type_analysis"]["l2_parent_l1"], split_by)
    output_dir.mkdir(parents=True, exist_ok=True)
    records, used_slugs = [], set()
    for label in sorted(labels.unique().tolist()):
        slug = _slug(str(label))
        if slug in used_slugs:
            raise ValueError(f"Cell-type labels produce the same output name: {label!r} -> {slug!r}")
        used_slugs.add(slug)
        subset = adata[labels == label].copy()
        # Retain only the global UMAP coordinates: the report uses them before
        # calculating the type-specific embedding.
        global_umap = subset.obsm.get("X_umap")
        subset.obsm.clear()
        if global_umap is not None:
            subset.obsm["X_umap"] = global_umap
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
        records.append({"cell_type": str(label), "slug": slug, "n_cells": int(subset.n_obs), "path": path.name})
    document = {"kind": "cell_type_split", "input": str(input_path), "split_by": split_by, "cell_types": records}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(document, indent=2) + "\n")
    return document


def _add_sample_cell_counts(adata, parent_map: dict[str, str], split_by: str) -> None:
    """Record retained-PBMC and L2-derived parent-compartment denominators before splitting."""
    groupby = ["study", "sample"]
    missing = set(groupby) - set(adata.obs)
    if missing:
        raise ValueError(f"Sample cell counts require columns: {sorted(missing)}")
    adata.obs["n_cells_in_sample"] = (
        adata.obs.groupby(groupby, observed=True)["study"].transform("size").astype(int)
    )
    # The merged L2 annotation is the sole taxonomy source for this analysis.
    # Do not mix in ``aifi_l1_majority``: it comes from a separate classifier
    # pass and can disagree with the refreshed merged L2 label. Map *every*
    # retained cell from L2 to its configured theoretical L1 parent first, then
    # count study × sample × parent on that derived assignment.
    parents = adata.obs[split_by].astype(str).map(parent_map)
    if parents.isna().any():
        unknown = sorted(adata.obs.loc[parents.isna(), split_by].astype(str).unique())
        raise ValueError(f"Missing configured AIFI L2 → L1 parent mapping for: {unknown}")
    adata.obs["aifi_l1_parent_for_l2"] = parents.astype(str)
    parent_groupby = [*groupby, "aifi_l1_parent_for_l2"]
    adata.obs["n_cells_in_sample_l1_parent"] = (
        adata.obs.groupby(parent_groupby, observed=True)["study"].transform("size").astype(int)
    )
    # Future: add other L2-derived n_cells_in_sample_compartment denominators
    # before splitting to support additional hierarchy views without external
    # inputs. Keep all such denominators on this one derived taxonomy.


def analyse_cell_type(input_path: Path, output_dir: Path, pipeline: dict) -> dict:
    """Create a type-specific embedding, Leiden clusters, and descriptive tables."""
    import scanpy as sc

    adata = sc.read_h5ad(input_path)
    spec = pipeline["cell_type_analysis"]
    output_dir.mkdir(parents=True, exist_ok=True)
    cell_types = adata.obs["aifi_l2_majority"].astype(str).unique()
    parents = adata.obs["aifi_l1_parent_for_l2"].astype(str).unique()
    if len(cell_types) != 1 or len(parents) != 1:
        raise ValueError(f"{input_path}: expected exactly one configured L2 type and parent L1")
    report = {
        "kind": "cell_type_analysis", "status": "complete", "input": str(input_path),
        "cell_type": cell_types[0], "aifi_l1_parent": parents[0],
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
        import scanpy.external as sce

        # Integrate studies for the local neighbour graph, UMAP, and Leiden
        # clustering. This should be revisited as the downstream analysis
        # matures: correction can also remove biological age-associated signal.
        sce.pp.harmony_integrate(
            working,
            batch_key,
            basis="X_pca",
            adjusted_basis=integration["adjusted_basis"],
            random_state=pipeline["processing"]["random_seed"],
            verbose=False,
        )
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
    sc.tl.leiden(working, key_added="cluster", resolution=spec["clustering"]["resolution"], random_state=pipeline["processing"]["random_seed"])
    # Wilcoxon log fold-changes require non-negative log-normalized data.
    # ``working`` is scaled for PCA, so rank markers on the matching unscaled
    # normalized values to avoid undefined log-fold changes.
    marker_input = adata[:, selected].copy()
    marker_input.obs["cluster"] = working.obs["cluster"].astype(str).to_numpy()
    sc.tl.rank_genes_groups(marker_input, groupby="cluster", method="wilcoxon", use_raw=False)
    markers = _marker_table(sc.get.rank_genes_groups_df(marker_input, group=None), marker_input.obs["cluster"])
    markers.to_csv(output_dir / "markers.tsv", sep="\t", index=False)
    adata.obsm["X_pca"] = working.obsm["X_pca"]
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
    for index in range(adata.obsm["X_pca"].shape[1]):
        if valid.sum() < 3:
            correlation, pvalue = np.nan, np.nan
        else:
            correlation, pvalue = spearmanr(adata.obsm["X_pca"][valid, index], age[valid])
        rows.append({"pc": f"PC{index + 1}", "spearman_r": correlation, "pvalue": pvalue, "n_cells": int(valid.sum())})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Split a merged PBMC H5AD into primitive cell-type artifacts")
    subparsers = parser.add_subparsers(dest="command", required=True)
    split = subparsers.add_parser("split")
    split.add_argument("--input", type=Path, required=True)
    split.add_argument("--output-dir", type=Path, required=True)
    split.add_argument("--manifest", type=Path, required=True)
    split.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    pipeline = read_json(args.config)
    split_cell_types(args.input, args.output_dir, args.manifest, pipeline)


if __name__ == "__main__":
    main()

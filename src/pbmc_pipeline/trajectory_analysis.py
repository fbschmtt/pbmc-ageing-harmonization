"""Clustering and cross-cell-type summaries for fitted age trajectories."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from .config import AgeTrajectorySettings


def cluster_trajectory_profiles(
    profiles: pd.DataFrame,
    settings: AgeTrajectorySettings,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cluster complete trajectory vectors and return labels, UMAP coordinates, and means."""
    values = profiles.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    values = values.dropna(axis=0, how="any")
    annotations = pd.DataFrame(index=values.index)
    annotations.index.name = profiles.index.name
    if values.empty:
        annotations["trajectory_cluster"] = pd.Series(dtype="int64")
        annotations["umap_1"] = pd.Series(dtype="float64")
        annotations["umap_2"] = pd.Series(dtype="float64")
        return annotations, pd.DataFrame(columns=[*values.columns, "n_trajectories"])

    if len(values) == 1:
        cluster_ids = np.ones(1, dtype=int)
    else:
        tree = linkage(
            values.to_numpy(dtype=float),
            method=settings.linkage_method,
            metric=settings.distance_metric,
        )
        cluster_ids = fcluster(
            tree, t=min(settings.max_clusters, len(values)), criterion="maxclust"
        ).astype(int)
    annotations["trajectory_cluster"] = cluster_ids

    if len(values) >= settings.minimum_umap_trajectories:
        import umap

        embedding = umap.UMAP(
            n_components=2,
            n_neighbors=min(settings.umap_neighbors, len(values) - 1),
            min_dist=settings.umap_min_dist,
            metric=settings.distance_metric,
            random_state=settings.random_state,
            transform_seed=settings.random_state,
            n_jobs=1,
        ).fit_transform(values.to_numpy(dtype=float))
        annotations["umap_1"] = embedding[:, 0]
        annotations["umap_2"] = embedding[:, 1]
    else:
        annotations["umap_1"] = np.nan
        annotations["umap_2"] = np.nan

    clustered = values.assign(trajectory_cluster=cluster_ids)
    means = clustered.groupby("trajectory_cluster", sort=True).agg(
        {**{column: "mean" for column in values.columns},
         "trajectory_cluster": "size"}
    ).rename(columns={"trajectory_cluster": "n_trajectories"})
    return annotations, means


def _read_cell_type_trajectory_results(
    differential_expression_dir: Path,
) -> tuple[dict[str, pd.DataFrame], dict[str, dict]]:
    results: dict[str, pd.DataFrame] = {}
    metadata_by_type: dict[str, dict] = {}
    for metadata_path in sorted(differential_expression_dir.glob("*/run_metadata.json")):
        metadata = json.loads(metadata_path.read_text())
        fit = next(
            (model for model in metadata.get("models", [])
             if model.get("purpose") == "age_trajectory"),
            None,
        )
        if fit is None:
            continue
        metadata_by_type[str(metadata["cell_type"])] = metadata
        result_path = metadata_path.parent / "combined" / "age_bin.csv"
        if fit.get("status") == "complete" and result_path.is_file():
            results[str(metadata["cell_type"])] = pd.read_csv(result_path)
    return results, metadata_by_type


def build_cross_cell_type_trajectories(
    results_by_type: dict[str, pd.DataFrame],
    metadata_by_type: dict[str, dict],
    *,
    settings: AgeTrajectorySettings,
) -> dict[str, pd.DataFrame | list[str] | dict[str, str]]:
    """Build common-bin profiles, recurring-gene summaries, and support diagnostics."""
    common_bins: set[str] | None = None
    supported_types = []
    qualified_types = []
    skipped_types = {}
    for cell_type, metadata in metadata_by_type.items():
        support = metadata.get("age_trajectory", {})
        if cell_type not in results_by_type or support.get("status") != "complete":
            skipped_types[cell_type] = support.get("reason", "trajectory fit was not completed")
            continue
        retained = [str(value) for value in support.get("retained_bins", [])]
        if len(retained) < settings.minimum_bins:
            skipped_types[cell_type] = (
                f"only {len(retained)} age bins meet the sample-count threshold; "
                f"at least {settings.minimum_bins} are required"
            )
            continue
        if not retained:
            skipped_types[cell_type] = "no retained age bins are recorded"
            continue
        qualified_types.append(cell_type)
        common_bins = set(retained) if common_bins is None else common_bins.intersection(retained)
        supported_types.append(cell_type)
    bins = sorted(common_bins or [], key=lambda label: int(label.split("-", maxsplit=1)[0]))
    if len(bins) < settings.minimum_shared_bins:
        for cell_type in supported_types:
            skipped_types[cell_type] = (
                f"only {len(bins)} age bins are shared across fitted cell types; "
                f"at least {settings.minimum_shared_bins} are required"
            )
        supported_types = []
        bins = []

    profiles = []
    significance_records = []
    for cell_type, result in results_by_type.items():
        if cell_type not in qualified_types:
            continue
        result = result.copy()
        result["gene"] = result["gene"].astype(str)
        result["age_bin"] = result["age_bin"].astype(str)
        result["omnibus_padj"] = pd.to_numeric(result["omnibus_padj"], errors="coerce")
        for gene, gene_rows in result.groupby("gene", sort=False):
            unique_rows = gene_rows.drop_duplicates("age_bin").set_index("age_bin")
            omnibus_fdr = float(gene_rows["omnibus_padj"].dropna().iloc[0]) if (
                gene_rows["omnibus_padj"].notna().any()
            ) else np.nan
            significance_records.append({
                "gene": gene,
                "cell_type": cell_type,
                "omnibus_fdr": omnibus_fdr,
                "tested": bool(np.isfinite(omnibus_fdr)),
                "de_significant": bool(
                    np.isfinite(omnibus_fdr) and omnibus_fdr < settings.de_fdr_threshold
                ),
                "trajectory_cluster_significant": bool(
                    np.isfinite(omnibus_fdr)
                    and omnibus_fdr <= settings.cluster_fdr_threshold
                ),
            })
            if cell_type not in supported_types or not set(bins).issubset(unique_rows.index):
                continue
            if (
                not np.isfinite(omnibus_fdr)
                or omnibus_fdr > settings.cluster_fdr_threshold
            ):
                continue
            profile = pd.to_numeric(unique_rows.loc[bins, "log2FoldChange"], errors="coerce")
            if not np.isfinite(profile.to_numpy(dtype=float)).all():
                continue
            profile = profile.to_numpy(dtype=float)
            profile = profile - profile[0]
            profile_sd = float(np.std(profile, ddof=0))
            if not np.isfinite(profile_sd) or profile_sd <= 0:
                continue
            row = {"gene": gene, "cell_type": cell_type,
                   "omnibus_fdr": omnibus_fdr,
                   **dict(zip(bins, profile / profile_sd, strict=True))}
            profiles.append(row)

    trajectory_frame = pd.DataFrame(profiles)
    if not trajectory_frame.empty:
        trajectory_frame["trajectory_id"] = (
            trajectory_frame["gene"].astype(str) + " | "
            + trajectory_frame["cell_type"].astype(str)
        )
        trajectory_frame = trajectory_frame.set_index("trajectory_id", drop=False)
    significance = pd.DataFrame(significance_records)
    if significance.empty:
        recurrence = pd.DataFrame(columns=[
            "gene", "n_cell_types_tested", "n_cell_types_DE_significant",
            "significant_cell_types",
        ])
        pattern_summary = pd.DataFrame(columns=[
            "gene", "n_cell_types_cluster_significant", "dominant_cluster_fraction",
        ])
    else:
        recurrence_records = []
        for gene, rows in significance.groupby("gene", sort=False):
            recurrence_records.append({
                "gene": gene,
                "n_cell_types_tested": int(rows["tested"].sum()),
                "n_cell_types_DE_significant": int(rows["de_significant"].sum()),
                "significant_cell_types": "; ".join(sorted(
                    rows.loc[rows["de_significant"], "cell_type"].unique()
                )),
            })
        recurrence = pd.DataFrame(recurrence_records)
        recurrence = recurrence.sort_values(
            ["n_cell_types_DE_significant", "n_cell_types_tested", "gene"],
            ascending=[False, False, True], kind="stable",
        )
        pattern_summary = pd.DataFrame(columns=[
            "gene", "n_cell_types_cluster_significant", "dominant_cluster_fraction",
        ])
    return {
        "common_bins": bins,
        "supported_cell_types": supported_types,
        "eligible_cell_types": qualified_types,
        "skipped_cell_types": skipped_types,
        "trajectories": trajectory_frame,
        "significance": significance,
        "gene_recurrence": recurrence,
        "gene_pattern_summary": pattern_summary,
    }


def analyze_age_trajectories(
    differential_expression_dir: Path,
    output_dir: Path,
    settings: AgeTrajectorySettings,
) -> dict:
    """Write shared-bin cross-type clusters, recurrence, and support summaries."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for obsolete_name in (
        "cell_type_trajectory_clusters.csv", "cell_type_cluster_means.csv",
    ):
        (output_dir / obsolete_name).unlink(missing_ok=True)
    results_by_type, metadata_by_type = _read_cell_type_trajectory_results(
        differential_expression_dir
    )
    cell_type_support = []
    for cell_type, metadata in metadata_by_type.items():
        support = metadata.get("age_trajectory", {})
        fit = next(
            (model for model in metadata.get("models", [])
             if model.get("purpose") == "age_trajectory"),
            {},
        )
        retained_bins = [str(value) for value in support.get("retained_bins", [])]
        eligible_for_reporting = len(retained_bins) >= settings.minimum_bins
        cell_type_support.append({
            "cell_type": cell_type,
            "status": support.get("status", fit.get("status", "unknown")),
            "reason": support.get("reason", fit.get("reason")) if eligible_for_reporting else (
                f"only {len(retained_bins)} age bins meet the sample-count threshold; "
                f"at least {settings.minimum_bins} are required"
            ),
            "retained_bins": retained_bins,
            "age_bin_counts": support.get("age_bin_counts", {}),
        })
    cross = build_cross_cell_type_trajectories(
        results_by_type, metadata_by_type, settings=settings
    )
    common_bins = cross["common_bins"]
    cross_profiles = cross["trajectories"]
    if not cross_profiles.empty:
        cluster_profiles = cross_profiles.set_index("trajectory_id")[common_bins]
        cross_annotations, cross_means = cluster_trajectory_profiles(
            cluster_profiles,
            settings,
        )
        cross_profiles = cross_profiles.join(cross_annotations)
        cross_means = cross_means.reset_index()
    else:
        cross_annotations = pd.DataFrame()
        cross_means = pd.DataFrame(columns=[
            "trajectory_cluster", *common_bins, "n_trajectories",
        ])
        cross_profiles = pd.DataFrame(columns=[
            "trajectory_id", "gene", "cell_type", "omnibus_fdr", *common_bins,
            "trajectory_cluster", "umap_1", "umap_2",
        ])

    pattern_rows = []
    if not cross_profiles.empty:
        for gene, group in cross_profiles.groupby("gene", sort=False):
            if len(group) < 2:
                continue
            vectors = group[common_bins].to_numpy(dtype=float)
            distances = pdist(vectors, metric=settings.distance_metric)
            cluster_counts = group["trajectory_cluster"].value_counts()
            pattern_rows.append({
                "gene": gene,
                "n_cell_types_cluster_significant": int(group["cell_type"].nunique()),
                "dominant_cluster_fraction": float(cluster_counts.max() / cluster_counts.sum()),
                "mean_pairwise_shape_distance": float(np.mean(distances)),
                "same_cluster_pair_fraction": float(np.mean([
                    group.iloc[first]["trajectory_cluster"]
                    == group.iloc[second]["trajectory_cluster"]
                    for first in range(len(group)) for second in range(first + 1, len(group))
                ])),
                "cell_types": "; ".join(sorted(group["cell_type"].unique())),
            })
    pattern_summary = pd.DataFrame(pattern_rows)
    if not pattern_summary.empty:
        pattern_summary = pattern_summary.sort_values(
            ["n_cell_types_cluster_significant", "same_cluster_pair_fraction",
             "mean_pairwise_shape_distance", "gene"],
            ascending=[False, False, True, True], kind="stable",
        )
    else:
        pattern_summary = pd.DataFrame(columns=[
            "gene", "n_cell_types_cluster_significant", "dominant_cluster_fraction",
            "mean_pairwise_shape_distance", "same_cluster_pair_fraction", "cell_types",
        ])

    significance = cross["significance"]
    recurrence = cross["gene_recurrence"]
    outputs = {
        "cross_cell_type_trajectory_clusters.csv": cross_profiles,
        "cross_cell_type_cluster_means.csv": cross_means,
        "gene_cell_type_significance.csv": significance,
        "gene_recurrence.csv": recurrence,
        "gene_pattern_concordance.csv": pattern_summary,
    }
    for name, table in outputs.items():
        table.to_csv(output_dir / name, index=False)
    analysis_metadata = {
        **settings.to_mapping(),
        "minimum_retained_bins": settings.minimum_bins,
        "common_bins": common_bins,
        "cell_types_in_cross_cell_type_analysis": cross["supported_cell_types"],
        "cell_types_eligible_by_bin_support": cross["eligible_cell_types"],
        "cell_types_excluded_from_cross_cell_type_analysis": cross["skipped_cell_types"],
        "cell_type_support": cell_type_support,
        "n_cross_cell_type_trajectories": len(cross_profiles),
        "n_genes_tested_in_cross_cell_type_analysis": int(significance["gene"].nunique())
        if not significance.empty else 0,
    }
    (output_dir / "analysis_metadata.json").write_text(
        json.dumps(analysis_metadata, indent=2) + "\n"
    )
    return analysis_metadata

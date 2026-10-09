"""Clustering and cross-cell-type summaries for fitted age trajectories."""
from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from .config import AgeTrajectorySettings, ExpressionAtlasSettings
from .expression_atlas import analyze_expression_atlas
from .utils import slugify


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


def trajectory_merge_diagnostics(
    profiles: pd.DataFrame,
    settings: AgeTrajectorySettings,
    *,
    maximum_merges: int = 20,
) -> pd.DataFrame:
    """Return the final hierarchical merges for choosing a trajectory cut.

    Each row records the cost of changing from ``from_clusters`` to one fewer
    groups.  The final merges are the ones relevant to the small cluster-count
    cuts used in reports, and are retained in descending cluster-count order
    for direct plotting.
    """
    values = profiles.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    values = values.dropna(axis=0, how="any")
    columns = ["from_clusters", "to_clusters", "merge_height"]
    if len(values) < 2:
        return pd.DataFrame(columns=columns)
    if maximum_merges < 1:
        raise ValueError("maximum_merges must be positive")

    tree = linkage(
        values.to_numpy(dtype=float),
        method=settings.linkage_method,
        metric=settings.distance_metric,
    )
    n_profiles = len(values)
    all_merges = pd.DataFrame({
        "from_clusters": np.arange(n_profiles, 1, -1, dtype=int),
        "to_clusters": np.arange(n_profiles - 1, 0, -1, dtype=int),
        "merge_height": tree[:, 2],
    })
    return all_merges.tail(min(maximum_merges, len(all_merges))).reset_index(drop=True)


def cluster_sample_residuals(
    residuals: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    settings: AgeTrajectorySettings,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Cluster samples after centered PCA of their Pearson residual profiles."""
    values = residuals.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    values = values.dropna(axis=0, how="any")
    annotations = sample_metadata.reindex(values.index).copy()
    annotations.index.name = "sample_key"
    annotations["sample_key"] = annotations.index.astype(str)
    if len(values) < 2 or values.shape[1] < 1:
        annotations["residual_cluster"] = pd.Series(index=annotations.index, dtype="Int64")
        annotations["umap_1"] = np.nan
        annotations["umap_2"] = np.nan
        return annotations.reset_index(drop=True), pd.DataFrame(), 0

    n_components = min(
        settings.residual_pca_components, len(values) - 1, values.shape[1]
    )
    scores = PCA(
        n_components=n_components,
        svd_solver="randomized",
        random_state=settings.random_state,
    ).fit_transform(values.to_numpy(dtype=float))
    if len(scores) >= settings.minimum_umap_trajectories:
        import umap

        embedding = umap.UMAP(
            n_components=2,
            n_neighbors=min(settings.residual_umap_neighbors, len(scores) - 1),
            min_dist=settings.umap_min_dist,
            metric="euclidean",
            random_state=settings.random_state,
            transform_seed=settings.random_state,
            n_jobs=1,
        ).fit_transform(scores)
        annotations["umap_1"] = embedding[:, 0]
        annotations["umap_2"] = embedding[:, 1]
        candidates = range(2, min(settings.max_clusters, len(scores) // 3) + 1)
        best_labels = np.ones(len(scores), dtype=int)
        best_score = -np.inf
        for n_clusters in candidates:
            labels = KMeans(
                n_clusters=n_clusters, n_init=50, random_state=settings.random_state
            ).fit_predict(embedding)
            if np.bincount(labels).min() < 3:
                continue
            score = silhouette_score(embedding, labels)
            if score > best_score:
                best_labels = labels + 1
                best_score = score
        cluster_ids = best_labels
    else:
        annotations["umap_1"] = np.nan
        annotations["umap_2"] = np.nan
        cluster_ids = np.ones(len(scores), dtype=int)
    annotations["residual_cluster"] = cluster_ids
    pc_columns = [f"PC{index + 1}" for index in range(scores.shape[1])]
    pc_scores = pd.DataFrame(scores, index=values.index, columns=pc_columns)
    pc_scores.insert(0, "residual_cluster", cluster_ids)
    return annotations.reset_index(drop=True), pc_scores.reset_index(names="sample_key"), n_components


def summarize_residual_cluster_markers(
    residuals: pd.DataFrame, annotations: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, int | None]:
    """Compare each residual cluster with the largest cluster and track feature origins."""
    clusters = annotations.set_index("sample_key")["residual_cluster"].reindex(residuals.index)
    cluster_sizes = clusters.value_counts()
    if len(cluster_sizes) < 2:
        return pd.DataFrame(), pd.DataFrame(), None
    reference_cluster = int(cluster_sizes.index[0])
    reference_mean = residuals.loc[clusters == reference_cluster].mean(axis=0)
    marker_rows = []
    composition_rows = []
    for cluster, size in cluster_sizes.items():
        if int(cluster) == reference_cluster:
            continue
        difference = residuals.loc[clusters == cluster].mean(axis=0) - reference_mean
        for direction, ordered in (
            ("higher", difference.sort_values(ascending=False)),
            ("lower", difference.sort_values(ascending=True)),
        ):
            top = ordered.head(min(20, len(ordered)))
            for rank, (feature, value) in enumerate(top.items(), start=1):
                marker_rows.append({
                    "residual_cluster": int(cluster),
                    "reference_cluster": reference_cluster,
                    "cluster_size": int(size),
                    "direction": direction,
                    "rank": rank,
                    "feature": str(feature),
                    "mean_residual_difference_vs_reference": float(value),
                })
            top_features = ordered.head(min(500, len(ordered))).index.to_series()
            cell_types = top_features.str.split("::", n=1).str[0]
            for cell_type, count in cell_types.value_counts().items():
                composition_rows.append({
                    "residual_cluster": int(cluster),
                    "reference_cluster": reference_cluster,
                    "cluster_size": int(size),
                    "direction": direction,
                    "cell_type": cell_type,
                    "n_features_in_top_500": int(count),
                    "fraction_of_top_500_features": float(count / len(top_features)),
                })
    return pd.DataFrame(marker_rows), pd.DataFrame(composition_rows), reference_cluster


def _read_cell_type_residuals(
    differential_expression_dir: Path,
    settings: AgeTrajectorySettings,
    output_dir: Path,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, dict]]:
    for pattern in ("*_residual_clusters.csv", "*_residual_pc_scores.csv"):
        for stale_path in output_dir.glob(pattern):
            stale_path.unlink()
    for name in (
        "cross_cell_type_pearson_residuals.npz",
        "cross_cell_type_residual_sample_metadata.csv",
        "cross_cell_type_residual_clusters.csv",
        "cross_cell_type_residual_pc_scores.csv",
        "cross_cell_type_residual_cluster_markers.csv",
        "cross_cell_type_residual_cluster_feature_composition.csv",
    ):
        (output_dir / name).unlink(missing_ok=True)
    residuals_by_type: dict[str, pd.DataFrame] = {}
    metadata_by_type: dict[str, pd.DataFrame] = {}
    support_by_type: dict[str, dict] = {}
    for metadata_path in sorted(differential_expression_dir.glob("*/run_metadata.json")):
        metadata = json.loads(metadata_path.read_text())
        fit = next((model for model in metadata.get("models", [])
                    if model.get("purpose") == "age_trajectory"), None)
        if (
            fit is None or fit.get("status") != "complete"
            or not fit.get("pearson_residuals")
            or not fit.get("pearson_residual_sample_metadata")
        ):
            continue
        residual_path = metadata_path.parent / fit["pearson_residuals"]
        samples_path = metadata_path.parent / fit["pearson_residual_sample_metadata"]
        if not residual_path.is_file() or not samples_path.is_file():
            continue
        with np.load(residual_path, allow_pickle=False) as stored:
            residual_matrix = pd.DataFrame(
                stored["residuals"],
                index=stored["pseudobulk_id"].astype(str),
                columns=stored["genes"].astype(str),
            )
        sample_metadata = pd.read_csv(samples_path, dtype={"pseudobulk_id": str})
        sample_metadata["sample_key"] = (
            sample_metadata["study"].astype(str) + "|" + sample_metadata["sample"].astype(str)
        )
        residual_matrix.index = sample_metadata.set_index("pseudobulk_id").loc[
            residual_matrix.index, "sample_key"
        ].to_numpy()
        residual_matrix = residual_matrix.loc[~residual_matrix.index.duplicated(keep="first")]
        sample_metadata = sample_metadata.drop_duplicates("sample_key").set_index("sample_key")
        residuals_by_type[str(metadata["cell_type"])] = residual_matrix
        metadata_by_type[str(metadata["cell_type"])] = sample_metadata
        support_by_type[str(metadata["cell_type"])] = {
            "n_samples": len(residual_matrix),
            "n_genes": residual_matrix.shape[1],
            "n_pcs_requested": settings.residual_pca_components,
        }

    for cell_type, matrix in residuals_by_type.items():
        annotations, pc_scores, n_components = cluster_sample_residuals(
            matrix, metadata_by_type[cell_type], settings
        )
        annotations.to_csv(
            output_dir / f"{slugify(cell_type)}_residual_clusters.csv", index=False
        )
        pc_scores.to_csv(output_dir / f"{slugify(cell_type)}_residual_pc_scores.csv", index=False)
        support_by_type[cell_type]["n_pcs_used"] = n_components

    cross_types = sorted(residuals_by_type)
    all_samples = set().union(*(set(matrix.index) for matrix in residuals_by_type.values()))
    minimum_coverage = settings.cross_cell_type_residual_minimum_sample_coverage
    removed_cell_types = []

    def common_samples_for(types):
        return set.intersection(*(set(residuals_by_type[cell_type].index) for cell_type in types))

    while cross_types:
        common_samples = common_samples_for(cross_types)
        if not all_samples or len(common_samples) / len(all_samples) >= minimum_coverage:
            break
        if len(cross_types) == 1:
            break
        to_remove = min(
            cross_types,
            key=lambda cell_type: (len(residuals_by_type[cell_type]), cell_type.casefold()),
        )
        cross_types.remove(to_remove)
        remaining_common = common_samples_for(cross_types)
        removed_cell_types.append({
            "cell_type": to_remove,
            "n_samples_available": len(residuals_by_type[to_remove]),
            "n_samples_common_after_removal": len(remaining_common),
            "sample_coverage_after_removal": (
                len(remaining_common) / len(all_samples) if all_samples else 0.0
            ),
        })
    if cross_types:
        common_samples = sorted(common_samples_for(cross_types))
        cross_residuals = {cell_type: residuals_by_type[cell_type] for cell_type in cross_types}
        matrices = []
        for cell_type, matrix in cross_residuals.items():
            wide = matrix.loc[common_samples].copy()
            wide.columns = [f"{cell_type}::{gene}" for gene in wide.columns]
            matrices.append(wide)
        cross_matrix = pd.concat(matrices, axis=1)
        sample_metadata = metadata_by_type[cross_types[0]].reindex(common_samples)
        np.savez_compressed(
            output_dir / "cross_cell_type_pearson_residuals.npz",
            residuals=cross_matrix.to_numpy(dtype=np.float32),
            sample_key=np.asarray(cross_matrix.index.astype(str), dtype=np.str_),
            genes=np.asarray(cross_matrix.columns.astype(str), dtype=np.str_),
        )
        sample_metadata.to_csv(output_dir / "cross_cell_type_residual_sample_metadata.csv")
        annotations, pc_scores, n_components = cluster_sample_residuals(
            cross_matrix, sample_metadata, settings
        )
        annotations.to_csv(output_dir / "cross_cell_type_residual_clusters.csv", index=False)
        pc_scores.to_csv(output_dir / "cross_cell_type_residual_pc_scores.csv", index=False)
        markers, composition, reference_cluster = summarize_residual_cluster_markers(
            cross_matrix, annotations
        )
        markers.to_csv(
            output_dir / "cross_cell_type_residual_cluster_markers.csv", index=False
        )
        composition.to_csv(
            output_dir / "cross_cell_type_residual_cluster_feature_composition.csv", index=False
        )
        cross_support = {
            "n_cell_types": len(cross_residuals),
            "cell_types": cross_types,
            "removed_cell_types_for_sample_coverage": removed_cell_types,
            "n_samples_available_by_cell_type": {
                cell_type: len(matrix) for cell_type, matrix in cross_residuals.items()
            },
            "n_samples_common": len(common_samples),
            "n_samples_union": len(all_samples),
            "sample_coverage_fraction": (
                len(common_samples) / len(all_samples) if all_samples else 0.0
            ),
            "minimum_sample_coverage": minimum_coverage,
            "minimum_sample_coverage_met": (
                len(common_samples) / len(all_samples) >= minimum_coverage
                if all_samples else False
            ),
            "n_genes": cross_matrix.shape[1],
            "n_pcs_requested": settings.residual_pca_components,
            "n_pcs_used": n_components,
            "residual_cluster_method": "K-means on UMAP coordinates; silhouette-selected with at least 3 samples per cluster",
            "residual_marker_reference_cluster": reference_cluster,
        }
    else:
        cross_support = {
            "n_cell_types": 0, "cell_types": [], "n_samples_common": 0,
            "n_samples_union": len(all_samples), "n_genes": 0,
            "minimum_sample_coverage": minimum_coverage,
            "minimum_sample_coverage_met": False,
            "removed_cell_types_for_sample_coverage": removed_cell_types,
        }
    return residuals_by_type, metadata_by_type, {
        "cell_types": support_by_type,
        "cross_cell_type": cross_support,
    }


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
        retained = [
            str(value) for value in support.get("retained_bins", [])
            if int(str(value).split("-", maxsplit=1)[0])
            < settings.strict_age_cutoff_exclusive
        ]
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
            effect_values = pd.to_numeric(
                unique_rows["log2FoldChange"], errors="coerce"
            ).replace([np.inf, -np.inf], np.nan)
            if effect_values.notna().any():
                largest_effect_bin = str(effect_values.abs().idxmax())
                largest_effect = float(effect_values.loc[largest_effect_bin])
            else:
                largest_effect_bin = None
                largest_effect = np.nan
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
                "largest_absolute_log2_fold_change": largest_effect,
                "largest_effect_age_bin": largest_effect_bin,
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
                   "omnibus_fdr": omnibus_fdr, "trajectory_sd": profile_sd,
                   **dict(zip(bins, profile / profile_sd, strict=True)),
                   **{f"raw_{age_bin}": value for age_bin, value in zip(bins, profile, strict=True)}}
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
    *,
    pseudobulk_path: Path | None = None,
    expression_atlas_settings: ExpressionAtlasSettings | None = None,
    split_by: str = "aifi_l2_majority",
    l2_parent_l1: Mapping[str, str] | None = None,
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
            "samples_in_manually_excluded_bins": support.get(
                "samples_in_manually_excluded_bins", 0
            ),
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
    _, _, residual_support = _read_cell_type_residuals(
        differential_expression_dir, settings, output_dir,
    )
    if pseudobulk_path is not None and expression_atlas_settings is not None:
        expression_atlas = analyze_expression_atlas(
            pseudobulk_path, output_dir, expression_atlas_settings, split_by=split_by,
            l2_parent_l1=l2_parent_l1,
        )
    else:
        expression_atlas = {
            "status": "not_run",
            "reason": "merged pseudobulk input was not supplied to the trajectory report",
        }
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
        "pearson_residual_clustering": residual_support,
        "expression_atlas": expression_atlas,
    }
    (output_dir / "analysis_metadata.json").write_text(
        json.dumps(analysis_metadata, indent=2) + "\n"
    )
    return analysis_metadata

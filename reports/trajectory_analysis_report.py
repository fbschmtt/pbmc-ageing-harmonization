# %% [markdown]
# # Cross-cell-type age trajectories
#
# This report clusters standardized, omnibus-significant gene trajectories on
# shared age bins. It also summarizes how often each gene is DE-significant
# across cell types using each cell type's own omnibus FDR result.

# %%
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import Markdown, display

sns.set_theme(style="whitegrid")
analysis_dir = Path(os.environ["TRAJECTORY_ANALYSIS_DIR"])
metadata = json.loads((analysis_dir / "analysis_metadata.json").read_text())
cross_trajectories = pd.read_csv(analysis_dir / "cross_cell_type_trajectory_clusters.csv")
cluster_means = pd.read_csv(analysis_dir / "cross_cell_type_cluster_means.csv")
gene_recurrence = pd.read_csv(analysis_dir / "gene_recurrence.csv")
pattern_concordance = pd.read_csv(analysis_dir / "gene_pattern_concordance.csv")
cell_type_support = pd.DataFrame(metadata.get("cell_type_support", []))
common_bins = metadata.get("common_bins", [])
cluster_fdr = metadata["cluster_fdr_threshold"]
de_fdr = metadata["de_fdr_threshold"]
distance_metric = metadata["distance_metric"]
linkage_method = metadata["linkage_method"]
minimum_umap_trajectories = metadata["minimum_umap_trajectories"]
cross_report_top_n_genes = metadata["cross_report_top_n_genes"]
residual_support = metadata.get("pearson_residual_clustering", {})
cross_residual_support = residual_support.get("cross_cell_type", {})
expression_atlas_support = metadata.get("expression_atlas", {})
expression_matrix_path = analysis_dir / "expression_atlas_matrix.csv"
expression_clusters_path = analysis_dir / "expression_atlas_gene_clusters.csv"
expression_cluster_means_path = analysis_dir / "expression_atlas_cluster_means.csv"
expression_technology_contrast_path = analysis_dir / "expression_atlas_technology_contrast.csv"
expression_technology_summary_path = analysis_dir / "expression_atlas_technology_summary.csv"
expression_intronic_contrast_path = analysis_dir / "expression_atlas_intronic_contrast.csv"
expression_intronic_summary_path = analysis_dir / "expression_atlas_intronic_summary.csv"


def read_optional_csv(path, **kwargs):
    if not path.is_file() or not path.stat().st_size:
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


expression_matrix = read_optional_csv(expression_matrix_path, index_col="gene")
expression_clusters = read_optional_csv(expression_clusters_path)
expression_cluster_means = read_optional_csv(expression_cluster_means_path)
expression_technology_contrast = read_optional_csv(expression_technology_contrast_path)
expression_technology_summary = read_optional_csv(expression_technology_summary_path)
expression_intronic_contrast = read_optional_csv(expression_intronic_contrast_path)
expression_intronic_summary = read_optional_csv(expression_intronic_summary_path)

# %% [markdown]
# ## Support and scope

# %%
display(Markdown(
    f"The cross-cell-type trajectory comparison uses **{len(common_bins)} shared bins** "
    f"({', '.join(common_bins) if common_bins else 'none'}). "
    f"{len(metadata.get('cell_types_in_cross_cell_type_analysis', []))} cell types contribute; "
    f"{len(metadata.get('cell_types_excluded_from_cross_cell_type_analysis', {}))} are excluded. "
    f"Trajectory clustering uses omnibus FDR ≤ {cluster_fdr:g}; per-gene DE recurrence "
    f"counts use omnibus FDR < {de_fdr:g} within each cell type."
))
display(Markdown(
    f"Age-trajectory fits manually exclude the 90–100 bin and apply a strict age cutoff "
    f"of < {metadata['strict_age_cutoff_exclusive']} years. Each retained bin must contain "
    f"at least {metadata['minimum_samples_per_bin']} eligible pseudobulks. The support table "
    "lists the pre-filter sample count in the manually excluded bin for every cell type."
))
if not cell_type_support.empty:
    display(cell_type_support[[
        "cell_type", "status", "retained_bins",
        "samples_in_manually_excluded_bins",
    ]].fillna("—"))
excluded = metadata.get("cell_types_excluded_from_cross_cell_type_analysis", {})
if excluded:
    display(Markdown(
        f"{len(excluded)} cell types are omitted from trajectory summaries because "
        "they do not meet the age-bin support requirement or shared-bin threshold."
    ))

# %% [markdown]
# ## Study-balanced cell-type expression atlas
#
# This descriptive atlas is separate from differential expression. It uses only
# genes present in every study, sums counts within each study × cell type, drops
# groups below the configured depth threshold, computes `log2(CPM + pseudocount)`,
# and then averages those values equally across studies.

# %%
if expression_atlas_support.get("status") == "complete" and not expression_matrix.empty:
    display(Markdown(
        f"The atlas uses **{expression_atlas_support['n_shared_genes']:,} genes** shared by every "
        f"study. It retained **{expression_atlas_support['n_depth_qualified_groups']:,}** of "
        f"{expression_atlas_support['n_study_cell_type_groups']:,} study × cell-type groups with at least "
        f"{expression_atlas_support['minimum_study_cell_type_total_counts']:,} counts in that gene "
        f"intersection. Values are `log2(CPM + {expression_atlas_support['cpm_pseudocount']:g})`, "
        f"averaged equally over studies; each displayed cell type has at least "
        f"{expression_atlas_support['minimum_studies_per_cell_type']} retained studies. "
        "They are relative, study-balanced expression summaries rather than calibrated absolute RNA abundance."
    ))
    if not expression_clusters.empty:
        display(Markdown(
            f"Genes with CPM ≥ {expression_atlas_support['minimum_cpm_for_clustering']:g} in at least one "
            f"cell type were row-standardized and grouped into up to "
            f"{expression_atlas_support['max_clusters']} average-linkage correlation clusters. "
            f"Showing the {min(expression_atlas_support['report_top_n_genes'], len(expression_clusters)):,} "
            "most variable clustered genes."
        ))
        z_columns = [f"z_{column}" for column in expression_matrix.columns if f"z_{column}" in expression_clusters]
        heatmap_rows = expression_clusters.sort_values(
            ["expression_cluster", "profile_sd", "gene"], ascending=[True, False, True], kind="stable"
        ).head(expression_atlas_support["report_top_n_genes"])
        if z_columns and not heatmap_rows.empty:
            figure, axis = plt.subplots(
                figsize=(max(8, 0.95 * len(z_columns) + 3), max(5, 0.18 * len(heatmap_rows) + 2)),
            )
            sns.heatmap(
                heatmap_rows.set_index("gene")[z_columns].rename(columns=lambda name: name.removeprefix("z_")),
                cmap="vlag", center=0, yticklabels=True,
                cbar_kws={"label": "Expression profile (within-gene SD units)"}, ax=axis,
            )
            axis.set_xlabel("Cell type")
            axis.set_ylabel("Gene")
            axis.set_title("Cell-type expression-profile clusters")
            figure.tight_layout()
            plt.show()
        if not expression_cluster_means.empty:
            display(expression_cluster_means)
    display(Markdown(
        "Download the [study-balanced expression matrix](expression_atlas_matrix.csv), "
        "[study × cell-type depth support](expression_atlas_study_support.csv), and "
        "[gene-cluster assignments](expression_atlas_gene_clusters.csv)."
    ))
else:
    display(Markdown(
        "_The expression atlas was not produced: "
        + expression_atlas_support.get("reason", "no merged pseudobulk input was supplied")
        + "_"
    ))

# %% [markdown]
# ### 3′/5′ technology-associated contrast
#
# This sensitivity analysis calculates the difference between separately
# study-balanced 5′ and 3′ matrices. Technology is commonly study-confounded,
# so the contrast is descriptive and not an identified causal technology effect.

# %%
if not expression_technology_summary.empty and not expression_technology_contrast.empty:
    technology_support = expression_technology_contrast[[
        "cell_type", "n_studies_3_prime", "n_studies_5_prime",
        "meets_configured_study_replication",
    ]].drop_duplicates().sort_values("cell_type", kind="stable")
    display(Markdown(
        "This contrast is computed whenever at least one depth-qualified study exists in "
        "each technology family. The configured two-study threshold is a replication flag, "
        "not an estimability gate."
    ))
    display(technology_support)
    figure, axis = plt.subplots(figsize=(7.5, 4.8))
    sns.scatterplot(
        data=expression_technology_summary,
        x="mean_log2_cpm_across_technologies",
        y="median_log2_cpm_difference_5_prime_minus_3_prime",
        size="n_cell_types_with_technology_contrast",
        sizes=(8, 42), alpha=0.6, linewidth=0, legend=False, ax=axis,
    )
    axis.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Mean log2(CPM + pseudocount)")
    axis.set_ylabel("Median 5′ − 3′ log2(CPM + pseudocount)")
    axis.set_title("Technology-associated expression contrast")
    figure.tight_layout()
    plt.show()
    top_contrast_genes = expression_technology_summary.assign(
        absolute_difference=lambda table: table["median_log2_cpm_difference_5_prime_minus_3_prime"].abs()
    ).sort_values("absolute_difference", ascending=False, kind="stable").head(
        expression_atlas_support.get("report_top_n_genes", 200)
    )["gene"]
    contrast_heatmap = expression_technology_contrast.loc[
        expression_technology_contrast["gene"].isin(top_contrast_genes)
    ].pivot(index="gene", columns="cell_type", values="log2_cpm_difference_5_prime_minus_3_prime")
    if not contrast_heatmap.empty:
        contrast_heatmap = contrast_heatmap.reindex(top_contrast_genes.drop_duplicates())
        figure, axis = plt.subplots(
            figsize=(max(8, 0.95 * contrast_heatmap.shape[1] + 3), max(5, 0.18 * contrast_heatmap.shape[0] + 2)),
        )
        sns.heatmap(
            contrast_heatmap, cmap="vlag", center=0,
            cbar_kws={"label": "5′ − 3′ log2(CPM + pseudocount)"}, ax=axis,
        )
        axis.set_xlabel("Cell type")
        axis.set_ylabel("Gene")
        axis.set_title("Largest technology-associated contrasts")
        figure.tight_layout()
        plt.show()
    display(Markdown(
        "Download the [cell-type contrast matrix](expression_atlas_technology_contrast.csv) "
        "and [per-gene contrast/profile-concordance summary](expression_atlas_technology_summary.csv)."
    ))
else:
    technology_contrast_reason = expression_atlas_support.get("technology_contrast_reason")
    if not technology_contrast_reason:
        technology_contrast_reason = (
            "fewer than the configured number of depth-qualified studies supported one or both "
            "technology families for the same cell type"
        )
    display(Markdown(
        f"_No 3′/5′ contrast is available: {technology_contrast_reason}._"
    ))

# %% [markdown]
# ### Intronic-read inclusion contrast
#
# This is a separate, non-interaction sensitivity analysis: it contrasts studies
# whose alignment counted intronic reads with those that did not. It does not
# condition on or combine that label with the 3′/5′ technology contrast.

# %%
if not expression_intronic_summary.empty and not expression_intronic_contrast.empty:
    intronic_support = expression_intronic_contrast[[
        "cell_type", "n_studies_intronic", "n_studies_non_intronic",
        "meets_configured_study_replication",
    ]].drop_duplicates().sort_values("cell_type", kind="stable")
    display(Markdown(
        "This independent contrast is `intronic − non-intronic` on separately "
        "study-balanced log2(CPM + pseudocount) matrices. As for 3′/5′, a one-versus-one "
        "comparison is estimable but is flagged as unreplicated."
    ))
    display(intronic_support)
    figure, axis = plt.subplots(figsize=(7.5, 4.8))
    sns.scatterplot(
        data=expression_intronic_summary,
        x="mean_log2_cpm_across_intronic_status",
        y="median_log2_cpm_difference_intronic_minus_non_intronic",
        size="n_cell_types_with_intronic_contrast",
        sizes=(8, 42), alpha=0.6, linewidth=0, legend=False, ax=axis,
    )
    axis.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Mean log2(CPM + pseudocount)")
    axis.set_ylabel("Median intronic − non-intronic log2(CPM + pseudocount)")
    axis.set_title("Intronic-read inclusion expression contrast")
    figure.tight_layout()
    plt.show()
    top_intronic_genes = expression_intronic_summary.assign(
        absolute_difference=lambda table: table[
            "median_log2_cpm_difference_intronic_minus_non_intronic"
        ].abs()
    ).sort_values("absolute_difference", ascending=False, kind="stable").head(
        expression_atlas_support.get("report_top_n_genes", 200)
    )["gene"]
    intronic_heatmap = expression_intronic_contrast.loc[
        expression_intronic_contrast["gene"].isin(top_intronic_genes)
    ].pivot(
        index="gene", columns="cell_type",
        values="log2_cpm_difference_intronic_minus_non_intronic",
    )
    if not intronic_heatmap.empty:
        intronic_heatmap = intronic_heatmap.reindex(top_intronic_genes.drop_duplicates())
        figure, axis = plt.subplots(
            figsize=(max(8, 0.95 * intronic_heatmap.shape[1] + 3), max(5, 0.18 * intronic_heatmap.shape[0] + 2)),
        )
        sns.heatmap(
            intronic_heatmap, cmap="vlag", center=0,
            cbar_kws={"label": "Intronic − non-intronic log2(CPM + pseudocount)"}, ax=axis,
        )
        axis.set_xlabel("Cell type")
        axis.set_ylabel("Gene")
        axis.set_title("Largest intronic-read inclusion contrasts")
        figure.tight_layout()
        plt.show()
    display(Markdown(
        "Download the [cell-type intronic contrast matrix](expression_atlas_intronic_contrast.csv) "
        "and [per-gene intronic contrast summary](expression_atlas_intronic_summary.csv)."
    ))
else:
    intronic_contrast_reason = expression_atlas_support.get("intronic_contrast_reason")
    if not intronic_contrast_reason:
        intronic_contrast_reason = "no cell type had both depth-qualified intronic and non-intronic studies"
    display(Markdown(
        f"_No intronic-read inclusion contrast is available: {intronic_contrast_reason}._"
    ))

# %% [markdown]
# ## Sample clustering from bins-model Pearson residuals
#
# Per-cell-type matrices contain one sample per row and one gene per column.
# Values are Pearson residuals from the fitted age-bin model, including its
# estimable age-bin, sex, study/site, and log10(total_counts) terms. PCA centers
# the residual matrix and retains up to the configured number of components.
# The cross-cell-type matrix iteratively drops the cell type with the fewest
# fitted samples until the common sample set reaches the configured coverage.

# %%
residual_parameters = pd.DataFrame([
    ("Residual model", "Age bin + sex + log10(total counts) + study site when estimable"),
    ("Residual feature filter", "Finite fitted mean/dispersion and positive fitted variance; no HVG filter"),
    ("PCA", f"Centered, unscaled Pearson residuals; up to {metadata['residual_pca_components']} components"),
    ("Hierarchical clustering", f"{linkage_method} linkage; Euclidean distance in PC space; at most {metadata['max_clusters']} clusters"),
    ("Residual UMAP", f"2D Euclidean PC-space UMAP; {metadata['residual_umap_neighbors']} neighbors (capped at samples − 1); min_dist={metadata['umap_min_dist']}; seed={metadata['random_state']}"),
], columns=["Parameter", "Current value"])
display(residual_parameters)

# %%
def plot_residual_covariates(table, title, covariates, *, categorical_covariates):
    figure, axes = plt.subplots(2, 2, figsize=(16.8, 8.64))
    for axis, covariate in zip(axes.flat, covariates):
        if covariate not in table:
            axis.text(0.5, 0.5, "Not recorded", ha="center", va="center")
            axis.set_axis_off()
            continue
        values = table.dropna(subset=["umap_1", "umap_2", covariate]).copy()
        if values.empty:
            axis.text(0.5, 0.5, "No recorded values", ha="center", va="center")
            axis.set_axis_off()
            continue
        categorical = covariate in categorical_covariates
        sns.scatterplot(
            data=values, x="umap_1", y="umap_2", hue=covariate,
            palette="tab20" if categorical else "viridis",
            s=32, alpha=0.84, linewidth=0, legend="brief", ax=axis,
        )
        axis.set_title(f"{covariate.replace('_', ' ')} · {len(values)}/{len(table)} present")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        if categorical and axis.legend_ is not None:
            axis.legend_.set_title(covariate)
            axis.legend_.set_bbox_to_anchor((1.02, 1))
    for axis in axes.flat[len(covariates):]:
        axis.set_axis_off()
    figure.suptitle(title, y=1.01)
    figure.tight_layout()
    plt.show()

cross_residual_path = analysis_dir / "cross_cell_type_residual_clusters.csv"
if cross_residual_path.is_file() and cross_residual_support.get("n_cell_types", 0):
    cross_residual_clusters = pd.read_csv(cross_residual_path)
    included_types = cross_residual_support.get("cell_types", [])
    removed_types = cross_residual_support.get(
        "removed_cell_types_for_sample_coverage", []
    )
    removed_summary = "; ".join(
        f"{item['cell_type']} ({item['n_samples_available']} fitted samples; "
        f"intersection then {item['n_samples_common_after_removal']} samples, "
        f"{item['sample_coverage_after_removal']:.1%} coverage)"
        for item in removed_types
    )
    residual_scope = (
        f"### Cross-cell-type residual clusters\n\nThe merged residual matrix uses "
        f"{cross_residual_support['n_cell_types']} cell types and "
        f"{cross_residual_support['n_samples_common']} samples present in every included type "
        f"({cross_residual_support.get('sample_coverage_fraction', 0):.1%} of the "
        f"{cross_residual_support.get('n_samples_union', 0)} samples available across candidate types; "
        f"target {cross_residual_support.get('minimum_sample_coverage', 0):.0%}). "
        f"Included cell types: {', '.join(included_types)}. "
        + (
            "Types removed iteratively from lowest fitted-sample count: "
            + removed_summary
            + ". "
            if removed_types else "No cell types were removed. "
        )
    )
    if not cross_residual_support.get("minimum_sample_coverage_met", False):
        residual_scope += (
            "The coverage target was not reached after reducing to one cell type. "
        )
    residual_scope += (
        f"It has {cross_residual_support['n_genes']:,} cell-type × gene features; PCA used "
        f"{cross_residual_support.get('n_pcs_used', 0)} components. See the [cross-cell-type "
        "cluster assignments](cross_cell_type_residual_clusters.csv) and [wide Pearson "
        "residual matrix](cross_cell_type_pearson_residuals.npz)."
    )
    display(Markdown(residual_scope))
    residual_cluster_sizes = (
        cross_residual_clusters["residual_cluster"].value_counts(sort=False)
        .rename_axis("residual_cluster").reset_index(name="n_samples")
        .sort_values("n_samples", ascending=False, kind="stable")
    )
    residual_cluster_sizes["fraction_of_samples"] = (
        residual_cluster_sizes["n_samples"] / len(cross_residual_clusters)
    )
    display(Markdown(
        "Residual-cluster sizes are shown explicitly so very small outlier groups can be "
        "distinguished from broad structure. They remain descriptive unless a measured "
        "covariate or follow-up QC supports an explanation."
    ))
    display(residual_cluster_sizes)
    if cross_residual_support.get("n_samples_common", 0) >= minimum_umap_trajectories:
        figure, axis = plt.subplots(figsize=(6.6, 4.8))
        sns.scatterplot(
            data=cross_residual_clusters, x="umap_1", y="umap_2", hue="residual_cluster",
            palette="tab20", s=32, alpha=0.84, linewidth=0, ax=axis,
        )
        axis.set_title("Cross-cell-type residual UMAP: residual clusters")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        if axis.legend_ is not None:
            axis.legend_.set_title("Residual cluster")
            axis.legend_.set_bbox_to_anchor((1.02, 1))
        figure.tight_layout()
        plt.show()
        plot_residual_covariates(
            cross_residual_clusters,
            "Cross-cell-type residual UMAPs: model covariates",
            ["sex", "age", "log10_total_counts", "study"],
            categorical_covariates={"sex", "study"},
        )
        if any(covariate in cross_residual_clusters for covariate in ("bmi", "cmv")):
            plot_residual_covariates(
                cross_residual_clusters,
                "Cross-cell-type residual UMAPs: BMI and CMV coverage",
                [covariate for covariate in ("bmi", "cmv")
                 if covariate in cross_residual_clusters],
                categorical_covariates={"cmv"},
            )
    else:
        display(Markdown(
            f"_The selected set has {cross_residual_support['n_samples_common']} samples; "
            f"at least {minimum_umap_trajectories} are required to compute UMAP._"
        ))
else:
    display(Markdown(
        "_No residual-bearing cell types are available for cross-cell-type clustering._"
    ))

# %% [markdown]
# ## Cross-cell-type trajectory clusters
#
# Each point is one gene × cell-type trajectory. Profiles are re-standardized
# over the shared bins so shape comparisons use the same age intervals. The
# configured clustering groups trajectories by shape, not pathways.

# %%
if not cross_trajectories.empty and not cluster_means.empty:
    display(Markdown(
        "Each point is one gene × cell-type trajectory. Profiles are "
        f"re-standardized over the shared bins; {linkage_method} hierarchical "
        f"clustering and UMAP use {distance_metric} distance between bin values. "
        "Cluster labels identify broad shape groups, not pathways."
    ))
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    colors = sns.color_palette("husl", n_colors=len(cluster_means))
    color_map = {
        str(int(row.trajectory_cluster)): colors[index]
        for index, row in cluster_means.iterrows()
    }
    for _, row in cluster_means.iterrows():
        cluster = str(int(row.trajectory_cluster))
        axes[0].plot(
            common_bins, [row[age_bin] for age_bin in common_bins],
            marker="o", linewidth=2.2, color=color_map[cluster],
            label=f"Cluster {cluster} (n={int(row.n_trajectories)})",
        )
    axes[0].axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Age bin")
    axes[0].set_ylabel("Mean standardized trajectory (SD units)")
    axes[0].set_title("Mean trajectory by cluster")
    axes[0].tick_params(axis="x", rotation=35)
    axes[0].legend(fontsize=8)

    points = cross_trajectories.dropna(subset=["umap_1", "umap_2"]).copy()
    if len(points) >= minimum_umap_trajectories:
        points["trajectory_cluster"] = points["trajectory_cluster"].astype(str)
        sns.scatterplot(
            data=points, x="umap_1", y="umap_2", hue="trajectory_cluster",
            palette=color_map, s=22, alpha=0.75, linewidth=0, ax=axes[1],
        )
        axes[1].set_title(f"Trajectory UMAP ({distance_metric} bin distance)")
        axes[1].set_xlabel("UMAP 1")
        axes[1].set_ylabel("UMAP 2")
    else:
        axes[1].text(
            0.5, 0.5,
            f"At least {minimum_umap_trajectories} trajectories are needed for UMAP",
            ha="center", va="center", transform=axes[1].transAxes,
        )
        axes[1].set_axis_off()
    figure.tight_layout()
    plt.show()
    display(Markdown(
        f"{len(cross_trajectories):,} significant gene × cell-type trajectories "
        f"were clustered into {len(cluster_means)} groups."
    ))
    display(Markdown(
        "Download the [trajectory assignments](cross_cell_type_trajectory_clusters.csv) "
        "or [cluster mean profiles](cross_cell_type_cluster_means.csv)."
    ))
else:
    display(Markdown(
        "_No cross-cell-type trajectory clusters were produced. This can happen "
        "when too few age bins are shared or no gene passes the clustering threshold._"
    ))

# %% [markdown]
# ## Genes associated in multiple cell types
#
# Counts use per-cell-type omnibus FDR < 0.05 across completed age-bin fits
# that meet the minimum age-bin support requirement.
# The count is descriptive, does not apply a second correction across cell types,
# and depends on the number of cell types in which that gene was tested.

# %%
if not gene_recurrence.empty:
    recurrent = gene_recurrence.loc[
        gene_recurrence["n_cell_types_DE_significant"] > 0
    ].copy()
    if recurrent.empty:
        display(Markdown("_No genes passed the per-cell-type omnibus FDR threshold._"))
    else:
        recurrence_counts = recurrent["n_cell_types_DE_significant"].value_counts().sort_index()
        figure, axis = plt.subplots(figsize=(7, 4))
        sns.barplot(
            x=recurrence_counts.index.astype(str), y=recurrence_counts.values,
            color="#4c72b0", ax=axis,
        )
        axis.set_xlabel("Cell types with significant age-bin association")
        axis.set_ylabel("Number of genes")
        axis.set_title("Cross-cell-type recurrence of age-associated genes")
        figure.tight_layout()
        plt.show()
        display(Markdown(f"Showing at most {cross_report_top_n_genes} genes."))
        display(recurrent.head(cross_report_top_n_genes).reset_index(drop=True))
else:
    display(Markdown("_No gene-level cell-type significance results are available._"))

# %% [markdown]
# ## Do recurring genes share trajectory shapes?
#
# For genes with cluster-significant trajectories in at least two cell types,
# the table reports the share of cell-type pairs assigned to the same global
# shape cluster and their mean pairwise shape distance. These are descriptive
# measures of pattern agreement, not a formal cross-cell-type test.

# %%
if pattern_concordance.empty:
    display(Markdown(
        "_No genes have cluster-significant trajectories in more than one "
        "cell type on the shared bins._"
    ))
else:
    display(Markdown(
        f"Pairwise shape distances use the configured {distance_metric} metric. "
        f"Showing at most {cross_report_top_n_genes} genes."
    ))
    display(pattern_concordance.head(cross_report_top_n_genes).reset_index(drop=True))

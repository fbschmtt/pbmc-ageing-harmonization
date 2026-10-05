# %% [markdown]
# # Cross-cell-type age trajectories
#
# This report clusters standardized, omnibus-significant gene trajectories on
# shared age bins. It also summarizes how often each gene is DE-significant
# across cell types using each cell type's own omnibus FDR result.

# %%
import json
import os
import re
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
eligible_cell_types = metadata.get("cell_types_eligible_by_bin_support", [])
common_bins = metadata.get("common_bins", [])
cluster_fdr = metadata["cluster_fdr_threshold"]
de_fdr = metadata["de_fdr_threshold"]
distance_metric = metadata["distance_metric"]
linkage_method = metadata["linkage_method"]
minimum_umap_trajectories = metadata["minimum_umap_trajectories"]
cross_report_top_n_genes = metadata["cross_report_top_n_genes"]
residual_support = metadata.get("pearson_residual_clustering", {})
residual_types = residual_support.get("cell_types", {})
cross_residual_support = residual_support.get("cross_cell_type", {})

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
if not cell_type_support.empty:
    cell_type_support = cell_type_support.loc[
        cell_type_support["cell_type"].isin(eligible_cell_types)
    ]
    if not cell_type_support.empty:
        display(cell_type_support[["cell_type", "status", "retained_bins"]].fillna("—"))
excluded = metadata.get("cell_types_excluded_from_cross_cell_type_analysis", {})
if excluded:
    display(Markdown(
        f"{len(excluded)} cell types are omitted from trajectory summaries because "
        "they do not meet the age-bin support requirement or shared-bin threshold."
    ))

# %% [markdown]
# ## Sample clustering from bins-model Pearson residuals
#
# Per-cell-type matrices contain one sample per row and one gene per column.
# Values are Pearson residuals from the fitted age-bin model, including its
# estimable age-bin, sex, study/site, and log10(total_counts) terms. PCA centers
# the residual matrix and retains up to the configured number of components.
# The cross-cell-type matrix uses the cell types that contribute to the
# full-trajectory merge and samples present in every one of those types.

# %%
def plot_residual_covariates(table, title, figsize=(14, 3.6)):
    covariates = ["sex", "age", "log10_total_counts", "study"]
    figure, axes = plt.subplots(1, 4, figsize=figsize)
    for axis, covariate in zip(axes, covariates, strict=True):
        values = table.dropna(subset=["umap_1", "umap_2", covariate]).copy()
        if values.empty:
            axis.text(0.5, 0.5, "No UMAP coordinates", ha="center", va="center")
            axis.set_axis_off()
            continue
        categorical = covariate in {"sex", "study"}
        sns.scatterplot(
            data=values, x="umap_1", y="umap_2", hue=covariate,
            palette="tab20" if categorical else "viridis",
            hue_norm=None, s=27, alpha=0.82, linewidth=0,
            legend="brief", ax=axis,
        )
        axis.set_title(covariate.replace("_", " "))
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        if categorical and axis.legend_ is not None:
            axis.legend_.set_title(covariate)
            axis.legend_.set_bbox_to_anchor((1.02, 1))
    figure.suptitle(title, y=1.03)
    figure.tight_layout()
    plt.show()


if residual_types:
    display(Markdown(
        f"Pearson residual matrices were clustered for **{len(residual_types)} cell types**. "
        f"PCA used up to {metadata['residual_pca_components']} components; each type's "
        "actual component count is recorded in `analysis_metadata.json`."
    ))
    for cell_type, support in sorted(residual_types.items()):
        slug = re.sub(r"[^a-z0-9]+", "-", cell_type.lower()).strip("-") or "unnamed"
        cluster_path = analysis_dir / f"{slug}_residual_clusters.csv"
        if not cluster_path.is_file():
            continue
        residual_clusters = pd.read_csv(cluster_path)
        display(Markdown(
            f"### {cell_type}\n\n{support['n_samples']} samples and "
            f"{support['n_genes']:,} genes; PCA used {support.get('n_pcs_used', 0)} components. "
            f"Cluster assignments are available in [{slug} residual clusters]"
            f"({slug}_residual_clusters.csv)."
        ))
        plot_residual_covariates(residual_clusters, f"{cell_type}: residual sample UMAPs")
else:
    display(Markdown("_No completed age-bin fits with saved Pearson residuals are available._"))

cross_residual_path = analysis_dir / "cross_cell_type_residual_clusters.csv"
if cross_residual_path.is_file() and cross_residual_support.get("n_samples_common", 0) >= minimum_umap_trajectories:
    cross_residual_clusters = pd.read_csv(cross_residual_path)
    display(Markdown(
        f"### Cross-cell-type residual clusters\n\nThe merged residual matrix uses "
        f"{cross_residual_support['n_cell_types']} cell types that also contribute to the "
        "full-trajectory merge and the "
        f"{cross_residual_support['n_samples_common']} samples present in every included type. "
        f"It has {cross_residual_support['n_genes']:,} cell-type × gene features; PCA used "
        f"{cross_residual_support.get('n_pcs_used', 0)} components. See the [cross-cell-type "
        "cluster assignments](cross_cell_type_residual_clusters.csv) and [wide Pearson "
        "residual matrix](cross_cell_type_pearson_residuals.npz)."
    ))
    plot_residual_covariates(
        cross_residual_clusters, "Cross-cell-type residual sample UMAPs", figsize=(14, 3.8)
    )
elif cross_residual_support.get("n_cell_types", 0):
    display(Markdown(
        "_The cross-cell-type residual set has fewer samples than required to compute UMAP._"
    ))
elif cross_residual_support.get("missing_cell_types"):
    display(Markdown(
        "_The cross-cell-type residual view was skipped because residual matrices are missing "
        "for one or more cell types in the full-trajectory merge: "
        f"{', '.join(cross_residual_support['missing_cell_types'])}._"
    ))
else:
    display(Markdown(
        "_No residual-bearing cell types also met the full-trajectory merge support criteria._"
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

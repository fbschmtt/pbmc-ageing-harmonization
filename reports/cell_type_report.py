# %% [markdown]
# # Per-cell-type residual-variation report
#
# These are descriptive single-cell diagnostics. In particular, cell-level age
# correlations are not independent-sample inference and must be followed by a
# sample/subject-level analysis for any biological claim.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from IPython.display import Markdown, display

sns.set_theme(style="whitegrid")
from pbmc_pipeline.cell_type import analyse_cell_type
from pbmc_pipeline.config import read_json
from pbmc_pipeline.reporting import sample_cell_type_fractions, sample_cluster_fractions

primitive = sc.read_h5ad(Path(os.environ["CELL_TYPE_H5AD"]))
output_dir = Path(os.environ["CELL_TYPE_OUTPUT_DIR"])
pipeline = read_json(Path(os.environ["CELL_TYPE_CONFIG"]))
sc.settings.verbosity = 0
report = analyse_cell_type(Path(os.environ["CELL_TYPE_H5AD"]), output_dir, pipeline)
adata = sc.read_h5ad(output_dir / "analysis.h5ad")
display(Markdown(f"## {primitive.obs['aifi_l2_majority'].iloc[0]}\n\n{primitive.n_obs:,} cells × {primitive.n_vars:,} genes"))
display(pd.Series(report, name="value").to_frame())


def plot_umap(adata, *, color, title=None):
    """Render comparable low-resolution, equal-aspect UMAP panels."""
    sc.pl.umap(adata, color=color, size=max(1, 50_000 / adata.n_obs), show=False)
    figure = plt.gcf()
    figure.set_size_inches(8, 6)
    figure.set_dpi(90)
    for axis in figure.axes:
        if axis.has_data():
            axis.set_aspect("equal", adjustable="box")
    if title:
        figure.axes[0].set_title(title)
    figure.tight_layout()
    plt.show()

# %% [markdown]
# ## Contents
#
# - [Location on the global embedding](#global-location)
# - [Covariate coverage](#covariate-coverage)
# - [Sample-level fraction across age](#sample-fractions)
# - [Cluster composition across age](#cluster-composition)
# - [Cluster markers](#cluster-markers)
# - [PC--age correlations](#pc-age)

# %% [markdown]
# <a id="global-location"></a>
# ### Location on the global embedding

# %%
if "X_umap" in primitive.obsm:
    plot_umap(primitive, color="study", title="Former global UMAP coordinates, restricted to this cell type")

# %% [markdown]
# <a id="covariate-coverage"></a>
# ## Covariate coverage among retained samples
#
# These plots use one row per study × sample, rather than one row per cell, to
# make coverage and likely study-level confounding visible without overweighting
# large libraries. All currently retained cohorts are healthy by construction;
# the disease panel remains useful as an explicit guard against future changes.

# %%
sample_covariates = (
    primitive.obs.groupby(["study", "sample"], observed=True)
    .agg(
        age=("age", "first"),
        sex=("sex", "first"),
        technology=("technology", "first"),
        disease_status=("disease_status", "first"),
    )
    .reset_index()
)
figure, axes = plt.subplots(2, 2, figsize=(15, 10))
sns.stripplot(data=sample_covariates, x="study", y="age", hue="sex", dodge=True, ax=axes[0, 0])
axes[0, 0].set_title("Age coverage by study and sex")
axes[0, 0].tick_params(axis="x", rotation=45)
for axis, column, title in [
    (axes[0, 1], "sex", "Sex coverage by study"),
    (axes[1, 0], "technology", "Technology coverage by study"),
    (axes[1, 1], "disease_status", "Disease-status coverage by study"),
]:
    sns.countplot(data=sample_covariates, x=column, hue="study", ax=axis)
    axis.set_title(title)
    axis.tick_params(axis="x", rotation=45)
figure.tight_layout()
plt.show()

# %% [markdown]
# <a id="sample-fractions"></a>
# ## Sample-level fraction across age
#
# Each point is one study × sample. The left panel uses all retained PBMCs as
# denominator; the right panel uses the manually configured AIFI-L1 parent
# derived from every cell's merged L2 call, never from a separate L1 classifier.
# Both are descriptive rather than age-effect tests.

# %%
fractions_pbmc = sample_cell_type_fractions(
    adata,
    denominator="n_cells_in_sample",
    fraction_name="fraction_of_retained_pbmc",
)
fractions_l1 = sample_cell_type_fractions(
    adata,
    denominator="n_cells_in_sample_l1_parent",
    fraction_name="fraction_within_aifi_l1_parent",
)
figure, axes = plt.subplots(1, 2, figsize=(16, 6), sharex=True)
studies = sorted(fractions_pbmc["study"].unique())
palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
parent = adata.obs["aifi_l1_parent_for_l2"].iloc[0]
for axis, fractions, fraction_column, title, ylabel in [
    (axes[0], fractions_pbmc, "fraction_of_retained_pbmc", "Fraction of all retained PBMCs", "Fraction of retained PBMCs"),
    (axes[1], fractions_l1, "fraction_within_aifi_l1_parent", f"Fraction within {parent}", f"Fraction of {parent}"),
]:
    sns.scatterplot(data=fractions, x="age", y=fraction_column, hue="study", palette=palette, s=55, ax=axis)
    for study, subset in fractions.groupby("study", observed=True):
        if len(subset) >= 2 and subset["age"].nunique() >= 2:
            sns.regplot(data=subset, x="age", y=fraction_column, scatter=False, ci=None, color=palette[study], line_kws={"linestyle": "--", "alpha": 0.8}, ax=axis)
    axis.set_ylabel(ylabel)
    axis.set_title(title)
figure.tight_layout()
plt.show()

# %%
if report["status"] != "complete":
    display(Markdown("_Too few cells for a type-specific embedding and clustering._"))
else:
    colors = ["cluster", "study", "age"]
    if "aifi_l3_majority" in adata.obs:
        colors.append("aifi_l3_majority")
    for color in colors:
        plot_umap(adata, color=color)

# %% [markdown]
# <a id="cluster-composition"></a>
# ## Cluster composition

# %%
if report["status"] == "complete":
    composition = pd.read_csv(output_dir / "cluster_study_composition.tsv", sep="\t")
    display(composition.pivot(index="cluster", columns="study", values="n_cells").fillna(0).astype(int))
    cluster_fractions = sample_cluster_fractions(adata)
    clusters = sorted(cluster_fractions["cluster"].astype(str).unique())
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5 * n_columns, 4 * n_rows), squeeze=False, sharex=True, sharey=True)
    studies = sorted(cluster_fractions["study"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    for axis, cluster in zip(axes.flat, clusters):
        subset = cluster_fractions[cluster_fractions["cluster"].astype(str) == cluster]
        sns.scatterplot(data=subset, x="age", y="fraction_within_cell_type", hue="study", palette=palette, s=45, ax=axis, legend=False)
        for study, study_subset in subset.groupby("study", observed=True):
            if len(study_subset) >= 2 and study_subset["age"].nunique() >= 2:
                sns.regplot(data=study_subset, x="age", y="fraction_within_cell_type", scatter=False, ci=None, color=palette[study], line_kws={"linestyle": "--", "alpha": 0.8}, ax=axis)
        axis.set_title(f"Cluster {cluster}")
        axis.set_ylabel("Fraction within cell type")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    figure.tight_layout()
    plt.show()

# %% [markdown]
# <a id="cluster-markers"></a>
# ## Cluster markers

# %%
if report["status"] == "complete":
    markers = pd.read_csv(output_dir / "markers.tsv", sep="\t")
    top_markers = (
        markers.sort_values(["group", "scores"], ascending=[True, False])
        .groupby("group", observed=True, group_keys=False)
        .head(12)
    )
    clusters = sorted(top_markers["group"].astype(str).unique())
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5 * n_columns, 4 * n_rows), squeeze=False)
    for axis, cluster in zip(axes.flat, clusters):
        subset = top_markers[top_markers["group"].astype(str) == cluster].sort_values("scores", ascending=False)
        sns.barplot(data=subset, x="scores", y="names", color="#4c72b0", ax=axis)
        axis.set_title(f"Cluster {cluster}: top markers")
        axis.set_xlabel("Wilcoxon score")
        axis.set_ylabel("")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    figure.tight_layout()
    plt.show()

# %% [markdown]
# <a id="pc-age"></a>
# ## Exploratory cell-level PC–age correlations

# %%
if report["status"] == "complete":
    pc_age = pd.read_csv(output_dir / "pc_age_correlations.tsv", sep="\t")
    display(pc_age.reindex(pc_age["spearman_r"].abs().sort_values(ascending=False).index).head(10))
    plt.figure(figsize=(8, 4))
    sns.barplot(data=pc_age, x="pc", y="spearman_r", color="#4c72b0")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("Spearman correlation with age (cells; descriptive)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

# %% [markdown]
# ## TODO: pseudobulk age-associated genes on the local UMAP
#
# Consider loading the strongest age-associated genes from a matched,
# sample-level pseudobulk differential-expression analysis and plotting their
# expression on this type's local UMAP, alongside a compact table of its top
# hits. This would help connect robust sample-level evidence back to local cell
# states. We should first decide whether this belongs in this report or is out
# of scope: it introduces a cross-workflow dependency and UMAP colouring can
# easily be over-interpreted as validation of a pseudobulk association.

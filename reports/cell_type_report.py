# %% [markdown]
# # Per-cell-type residual-variation report
#
# These are descriptive single-cell diagnostics. In particular, cell-level age
# correlations are not independent-sample inference and must be followed by a
# sample/subject-level analysis for any biological claim.

# %%
import os
import warnings
from pathlib import Path

# The report runs in a notebook kernel where ipywidgets is intentionally not a
# runtime dependency; tqdm's optional-progress warning is not analytically
# relevant.
warnings.filterwarnings("ignore", message="IProgress not found.*")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from IPython.display import Markdown, display
from statsmodels.nonparametric.smoothers_lowess import lowess

sns.set_theme(style="whitegrid")
from pbmc_pipeline.config import read_json
from pbmc_pipeline.reporting import sample_cell_type_fractions, sample_cluster_fractions

primitive = sc.read_h5ad(Path(os.environ["CELL_TYPE_H5AD"]))
analysis_path = Path(os.environ["CELL_TYPE_ANALYSIS_H5AD"])
output_dir = analysis_path.parent
pipeline = read_json(Path(os.environ["CELL_TYPE_CONFIG"]))
sc.settings.verbosity = 0
report = read_json(output_dir / "report.json")
adata = sc.read_h5ad(analysis_path)
de_dir = os.environ.get("CELL_TYPE_DIFFERENTIAL_EXPRESSION_DIR")
de_run_metadata = {}
if de_dir:
    de_metadata_path = Path(de_dir) / "run_metadata.json"
    if de_metadata_path.is_file():
        de_run_metadata = read_json(de_metadata_path)
display(Markdown(f"## {primitive.obs['aifi_l2_majority'].iloc[0]}\n\n{primitive.n_obs:,} cells × {primitive.n_vars:,} genes"))
display(pd.Series(report, name="value").to_frame())
if de_run_metadata.get("analysis_mode") == "test_only":
    display(Markdown(
        "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
        + de_run_metadata["interpretation_warning"]
    ))


def plot_umap(adata, *, color, title=None, label_clusters=False):
    """Render comparable low-resolution, equal-aspect UMAP panels."""
    plot_kwargs = {"color": color, "size": max(2, 100_000 / adata.n_obs), "show": False}
    if label_clusters:
        plot_kwargs["legend_loc"] = "on data"
    sc.pl.umap(adata, **plot_kwargs)
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


def has_age_span(data, *, min_samples=2):
    return len(data) >= min_samples and data["age"].nunique() >= 2


def facet_grid_shape(n_panels, *, max_columns=4):
    """Keep study panels compact while growing to any number of studies."""
    n_columns = min(max_columns, n_panels)
    return (n_panels + n_columns - 1) // n_columns, n_columns


def linear_fit(data, *, y_column, x_grid):
    """Return an OLS line on a supplied grid."""
    x = data["age"].to_numpy(dtype=float)
    y = data[y_column].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(x)), x])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return coefficients[0] + coefficients[1] * x_grid


def lowess_fit(data, *, y_column, x_grid):
    """Return a LOWESS curve on a supplied grid."""
    x = data["age"].to_numpy(dtype=float)
    y = data[y_column].to_numpy(dtype=float)

    def predict(sample_x, sample_y):
        curve = lowess(sample_y, sample_x, frac=0.67, it=3, return_sorted=True)
        unique_x, first_indices = np.unique(curve[:, 0], return_index=True)
        return np.interp(x_grid, unique_x, curve[first_indices, 1])

    return predict(x, y)


def plot_study_linear_trends(
    data, *, fraction_column, title, ylabel, palette, show_fit_legend=True,
):
    """Show within-study trends with and without samples younger than 20."""
    if data.empty:
        display(Markdown(f"_No samples meet the denominator cutoff for {title}._"))
        return
    studies = sorted(data["study"].unique())
    n_rows, n_columns = facet_grid_shape(len(studies))
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(4.2 * n_columns, 3.5 * n_rows),
        squeeze=False, sharex=True, sharey=True,
    )
    for axis, study in zip(axes.flat, studies):
        subset = data.loc[data["study"] == study]
        sns.scatterplot(
            data=subset, x="age", y=fraction_column, color=palette[study],
            s=24, alpha=0.7, edgecolor="none", ax=axis,
        )
        for include_all_ages, color, linestyle in [
            (False, "#111111", "-"),
            (True, "#777777", "--"),
        ]:
            fit_data = subset if include_all_ages else subset.loc[subset["age"] >= 20]
            if has_age_span(fit_data):
                sns.regplot(
                    data=fit_data, x="age", y=fraction_column, scatter=False, ci=None,
                    color=color,
                    line_kws={"linewidth": 2.4, "alpha": 0.95, "linestyle": linestyle},
                    ax=axis,
                )
        if not has_age_span(subset):
            axis.text(0.5, 0.08, "No age span for trend", transform=axis.transAxes,
                      ha="center", va="bottom", fontsize=9)
        axis.set_title(f"{study} (n={len(subset)})")
        axis.set_ylabel(ylabel)
    for axis in axes.flat[len(studies):]:
        axis.remove()
    figure.suptitle(title, y=0.99)
    if show_fit_legend:
        fit_handles = [
            plt.Line2D([], [], color="#111111", linewidth=2.4,
                       label="Linear fit: age ≥20"),
            plt.Line2D([], [], color="#777777", linewidth=2.4, linestyle="--",
                       label="Linear fit: all samples"),
        ]
        figure.legend(handles=fit_handles, loc="upper center", ncol=2,
                      bbox_to_anchor=(0.5, 0.94))
        figure.tight_layout(rect=(0, 0, 1, 0.87))
    else:
        figure.tight_layout(rect=(0, 0, 1, 0.94))
    plt.show()


def plot_combined_study_fits(fraction_specs):
    """Show separate linear and LOWESS study fits with sample-share-weighted means."""
    studies = sorted(set().union(*(set(data["study"].astype(str)) for data, _, _ in fraction_specs)))
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    for fit_name, fit_function in [
        ("linear", linear_fit),
        ("LOWESS", lowess_fit),
    ]:
        figure, axes = plt.subplots(
            1, len(fraction_specs), figsize=(7 * len(fraction_specs), 5.5), squeeze=False,
        )
        for axis, (data, fraction_column, label) in zip(axes.flat, fraction_specs):
            plot_data = data.assign(study=data["study"].astype(str))
            eligible = {
                study: subset
                for study, subset in plot_data.groupby("study", sort=True)
                if has_age_span(subset)
            }
            total_n = sum(len(subset) for subset in eligible.values())
            legend_handles = []
            if not eligible:
                axis.text(0.5, 0.5, "Insufficient age variation for fit", transform=axis.transAxes,
                          ha="center", va="center")
            else:
                x_grid = np.linspace(
                    min(subset["age"].min() for subset in eligible.values()),
                    max(subset["age"].max() for subset in eligible.values()),
                    160,
                )
                study_fits = []
                for study, subset in eligible.items():
                    weight = len(subset) / total_n
                    local_grid = x_grid[
                        (x_grid >= subset["age"].min()) & (x_grid <= subset["age"].max())
                    ]
                    fitted = fit_function(
                        subset, y_column=fraction_column, x_grid=local_grid,
                    )
                    axis.plot(local_grid, fitted, color=palette[study], linewidth=2)
                    study_fits.append((study, subset, weight))
                    legend_handles.append(plt.Line2D(
                        [], [], color=palette[study], linewidth=2,
                        label=f"{study} (n={len(subset)}, {weight:.0%})",
                    ))

                weighted_sum = np.zeros(len(x_grid))
                weight_sum = np.zeros(len(x_grid))
                for _, subset, weight in study_fits:
                    supported = (x_grid >= subset["age"].min()) & (x_grid <= subset["age"].max())
                    fitted = fit_function(
                        subset, y_column=fraction_column, x_grid=x_grid[supported],
                    )
                    weighted_sum[supported] += weight * fitted
                    weight_sum[supported] += weight
                supported = weight_sum > 0
                axis.plot(
                    x_grid[supported], weighted_sum[supported] / weight_sum[supported],
                    color="#171717", linewidth=3.2, zorder=5,
                )
                legend_handles.append(plt.Line2D(
                    [], [], color="#171717", linewidth=3.2,
                    label="Sample-share-weighted mean",
                ))
            axis.set_title(label)
            axis.set_xlabel("Age (years)")
            axis.set_ylabel(label)
            if legend_handles:
                axis.legend(handles=legend_handles, title="Study (n samples, share)",
                            loc="best", fontsize=8, title_fontsize=8)
        figure.suptitle(f"Within-study {fit_name} fraction fits", y=0.99)
        figure.text(
            0.5, 0.015,
            "Study share = n / total eligible n; weights renormalize among studies spanning each age. Other schemes may be explored.",
            ha="center", fontsize=9,
        )
        figure.tight_layout(rect=(0, 0.05, 1, 0.93))
        plt.show()


def plot_cluster_fractions(data, *, palette):
    """Keep study identity in points while showing one uncluttered pooled trend."""
    if data.empty:
        display(Markdown("_No samples meet the denominator cutoff for cluster fractions._"))
        return
    plot_data = data.assign(study=data["study"].astype(str))
    clusters = sorted(plot_data["cluster"].astype(str).unique())
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(5 * n_columns, 4 * n_rows), squeeze=False,
        sharex=True, sharey=True,
    )
    studies = sorted(plot_data["study"].unique())
    for axis, cluster in zip(axes.flat, clusters):
        subset = plot_data.loc[plot_data["cluster"].astype(str) == cluster]
        sns.scatterplot(
            data=subset, x="age", y="fraction_within_cell_type", hue="study", palette=palette,
            s=20, alpha=0.65, edgecolor="none", ax=axis, legend=False,
        )
        if has_age_span(subset):
            sns.regplot(
                data=subset, x="age", y="fraction_within_cell_type", scatter=False, ci=None,
                color="#1f1f1f", line_kws={"linewidth": 2.6}, ax=axis,
            )
        axis.set_title(f"Cluster {cluster}")
        axis.set_ylabel("Fraction within cell type")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study, markersize=5)
        for study in studies
    ]
    figure.legend(handles=handles, title="Study", loc="upper center", ncol=min(4, len(studies)))
    figure.suptitle("Cluster fractions: pooled linear trends (study-unadjusted)", y=0.98)
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    plt.show()


def directional_cluster_genes(markers, *, cluster, n_genes=8):
    """Return the strongest positive and negative genes for one cluster contrast."""
    subset = markers.loc[markers["group"].astype(str) == str(cluster)].copy()
    positive = subset.loc[subset["logfoldchanges"] > 0].nlargest(n_genes, "logfoldchanges")
    negative = subset.loc[subset["logfoldchanges"] < 0].nsmallest(n_genes, "logfoldchanges")
    result = pd.concat([negative, positive]).copy()
    result["direction"] = result["logfoldchanges"].ge(0).map({
        True: "Higher in cluster", False: "Lower in cluster",
    })
    return result.sort_values("logfoldchanges")


def plot_pca_study_and_technology(adata):
    """Show the first four native local PC scores with study and chemistry encoded."""
    scores = adata.obsm["X_pca"]
    if scores.shape[1] < 2:
        display(Markdown("_Fewer than two local PCs are available to plot._"))
        return

    pairs = [(0, 1)]
    if scores.shape[1] >= 4:
        pairs.append((2, 3))
    score_data = adata.obs[["study", "technology"]].copy()
    score_data["study"] = score_data["study"].astype(str)
    score_data["technology"] = score_data["technology"].astype(str)
    for component in {component for pair in pairs for component in pair}:
        score_data[f"PC {component + 1}"] = scores[:, component]

    studies = sorted(score_data["study"].unique())
    technologies = sorted(score_data["technology"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    marker_cycle = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "h"]
    technology_markers = {
        technology: marker_cycle[index % len(marker_cycle)]
        for index, technology in enumerate(technologies)
    }
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])

    figure, axes = plt.subplots(1, len(pairs), figsize=(8 * len(pairs), 6.4), squeeze=False)
    for axis, (x_component, y_component) in zip(axes.flat, pairs):
        sns.scatterplot(
            data=score_data,
            x=f"PC {x_component + 1}",
            y=f"PC {y_component + 1}",
            hue="study",
            style="technology",
            palette=palette,
            markers=technology_markers,
            s=max(8, 120_000 / adata.n_obs),
            alpha=0.65,
            edgecolor="none",
            legend=False,
            rasterized=True,
            ax=axis,
        )
        axis.set_title(f"PC {x_component + 1} vs. PC {y_component + 1}")
        axis.set_xlabel(f"PC {x_component + 1} ({variance_ratio[x_component]:.1%} variance)")
        axis.set_ylabel(f"PC {y_component + 1} ({variance_ratio[y_component]:.1%} variance)")

    study_handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study, markersize=6)
        for study in studies
    ]
    technology_handles = [
        plt.Line2D(
            [], [], marker=technology_markers[technology], linestyle="", color="0.35",
            label=technology, markersize=6,
        )
        for technology in technologies
    ]
    study_legend = figure.legend(
        handles=study_handles, title="Study", loc="lower center",
        bbox_to_anchor=(0.5, 0.09), ncol=min(4, len(studies)),
    )
    figure.add_artist(study_legend)
    figure.legend(
        handles=technology_handles, title="Technology", loc="lower center",
        bbox_to_anchor=(0.5, 0.0), ncol=min(5, len(technologies)),
    )
    figure.suptitle("Native local PCA scores by study and technology", y=0.99)
    figure.tight_layout(rect=(0, 0.2, 1, 0.95))
    plt.show()

# %% [markdown]
# ## Contents
#
# - [Location on the global embedding](#global-location)
# - [Sample coverage](#sample-coverage)
# - [Sample-level fraction across age](#sample-fractions)
# - [Cluster composition across age](#cluster-composition)
# - [Cluster markers](#cluster-markers)
# - [PCA scores by study and technology](#pc-study-technology)
# - [PC--age correlations](#pc-age)

# %% [markdown]
# <a id="global-location"></a>
# ### Location on the global embedding

# %%
if "X_umap" in primitive.obsm:
    plot_umap(primitive, color="study", title="Former global UMAP coordinates, restricted to this cell type")

# %% [markdown]
# <a id="sample-coverage"></a>
# ## Sample coverage among retained samples
#
# One row represents one study × sample. The compact table retains the metadata
# context for the fraction plots below without repeating separate count plots
# for fields that are usually constant within a study.

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
study_coverage = (
    sample_covariates.groupby("study", observed=True)
    .agg(
        n_samples=("sample", "nunique"),
        age_min=("age", "min"),
        age_max=("age", "max"),
        sex=("sex", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        technology=("technology", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        disease_status=("disease_status", lambda values: ", ".join(sorted(values.astype(str).unique()))),
    )
    .reset_index()
)
display(study_coverage)
figure, axis = plt.subplots(figsize=(10, 4))
sns.stripplot(data=sample_covariates, x="study", y="age", hue="sex", dodge=True, ax=axis)
axis.set_title("Age coverage by study and sex")
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
# Samples with fewer than the configured denominator-cell cutoff are excluded.
# Study-specific trends and weighted summaries are descriptive rather than
# age-effect tests; the summaries do not adjust for study.
# Facets show within-study linear fits for all samples and for age ≥20. Separate
# combined figures show unshaded per-study linear or LOWESS fits, plus a mean
# curve weighted by each study's share of eligible samples. At each age, the
# mean is renormalized over studies whose observed age range covers it.

# %%
min_fraction_denominator = pipeline["cell_type_analysis"]["min_fraction_denominator_cells"]
fraction_specs = []
for denominator, fraction_column, title, ylabel in [
    ("n_cells_in_sample", "fraction_of_retained_pbmc", "Fraction of all retained PBMCs", "Fraction of retained PBMCs"),
    ("n_cells_in_sample_l1_parent", "fraction_within_aifi_l1_parent", None, None),
]:
    fraction_specs.append((
        sample_cell_type_fractions(
            adata, denominator=denominator, fraction_name=fraction_column,
            min_denominator_cells=min_fraction_denominator,
        ),
        fraction_column,
        title,
        ylabel,
    ))
parent = adata.obs["aifi_l1_parent_for_l2"].iloc[0]
fraction_specs[1] = (
    fraction_specs[1][0], fraction_specs[1][1], f"Fraction within {parent}", f"Fraction of {parent}"
)
all_fraction_samples = [
    len(sample_cell_type_fractions(adata, denominator=denominator, fraction_name=fraction_column))
    for denominator, fraction_column, _, _ in [
        ("n_cells_in_sample", "fraction_of_retained_pbmc", None, None),
        ("n_cells_in_sample_l1_parent", "fraction_within_aifi_l1_parent", None, None),
    ]
]
display(pd.DataFrame({
    "fraction": [spec[2] for spec in fraction_specs],
    "eligible_samples": [len(spec[0]) for spec in fraction_specs],
    "excluded_below_denominator_cutoff": [
        total - len(spec[0]) for total, spec in zip(all_fraction_samples, fraction_specs)
    ],
    "minimum_denominator_cells": min_fraction_denominator,
}))
studies = sorted(pd.concat([spec[0] for spec in fraction_specs])["study"].unique())
palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
for index, (fractions, fraction_column, title, ylabel) in enumerate(fraction_specs):
    plot_study_linear_trends(
        fractions, fraction_column=fraction_column, title=title, ylabel=ylabel, palette=palette,
        show_fit_legend=index == 0,
    )
plot_combined_study_fits([
    (fractions, fraction_column, ylabel)
    for fractions, fraction_column, _, ylabel in fraction_specs
])

# %%
if report["status"] != "complete":
    display(Markdown("_Too few cells for a type-specific embedding and clustering._"))
else:
    plot_umap(adata, color="cluster", title="Local UMAP: cluster IDs", label_clusters=True)
    colors = ["study", "age"]
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
    all_cluster_fractions = sample_cluster_fractions(adata)
    cluster_fractions = sample_cluster_fractions(
        adata, min_denominator_cells=min_fraction_denominator,
    )
    excluded_cluster_samples = (
        all_cluster_fractions[["study", "sample"]].drop_duplicates().shape[0]
        - cluster_fractions[["study", "sample"]].drop_duplicates().shape[0]
    )
    display(Markdown(
        f"_{excluded_cluster_samples} sample(s) excluded from cluster fractions because the "
        f"split cell-type denominator has fewer than {min_fraction_denominator} cells._"
    ))
    studies = sorted(cluster_fractions["study"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    plot_cluster_fractions(cluster_fractions, palette=palette)

# %% [markdown]
# <a id="cluster-markers"></a>
# ## Cluster markers
#
# Genes are ranked by a Wilcoxon comparison of each local cluster against all
# other local clusters. Signed log fold changes therefore show whether a gene
# is higher or lower in that cluster; these are descriptive contrasts.

# %%
if report["status"] == "complete":
    markers = pd.read_csv(output_dir / "markers.tsv", sep="\t")
    clusters = sorted(markers["group"].astype(str).unique())
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5 * n_columns, 4 * n_rows), squeeze=False)
    for axis, cluster in zip(axes.flat, clusters):
        subset = directional_cluster_genes(markers, cluster=cluster)
        sns.barplot(
            data=subset, x="logfoldchanges", y="names", hue="direction", dodge=False,
            palette={"Higher in cluster": "#c44e52", "Lower in cluster": "#4c72b0"}, ax=axis,
        )
        axis.axvline(0, color="black", linewidth=0.8)
        axis.set_title(f"Cluster {cluster}: differential genes")
        axis.set_xlabel("Log fold change versus all other clusters")
        axis.set_ylabel("")
        axis.legend(title="Direction", fontsize=8, title_fontsize=8, loc="best")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    figure.tight_layout()
    plt.show()

# %% [markdown]
# <a id="pc-study-technology"></a>
# ## PCA scores by study and technology
#
# These are the native, local PCA scores before Harmony adjustment. Colour uses
# the report's standard study palette and marker shape indicates the recorded
# single-cell technology, so technical structure in the first four PCs is
# visible directly. The Harmony-adjusted representation is used only for the
# local neighbour graph, UMAP, and clustering.

# %%
if report["status"] == "complete":
    plot_pca_study_and_technology(adata)

# %% [markdown]
# <a id="pc-age"></a>
# ## PCA gene loadings, variance explained, and exploratory PC–age correlations

# %%
if report["status"] == "complete":
    n_components_to_plot = min(10, adata.obsm["X_pca"].shape[1])
    sc.pl.pca_loadings(
        adata, components=list(range(1, n_components_to_plot + 1)), n_points=12, show=False,
    )
    figure = plt.gcf()
    n_columns = min(4, n_components_to_plot)
    n_rows = (n_components_to_plot + n_columns - 1) // n_columns
    figure.set_size_inches(14, 3.0 * n_rows)
    figure.subplots_adjust(left=0.06, right=0.98, bottom=0.08, top=0.94, hspace=0.75, wspace=0.35)
    plt.show()

# %%
if report["status"] == "complete":
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
    cumulative_variance = np.cumsum(variance_ratio)
    figure, axis = plt.subplots(figsize=(8, 4))
    axis.plot(np.arange(1, len(cumulative_variance) + 1), cumulative_variance, marker="o", markersize=3)
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Cumulative variance explained")
    axis.set_ylim(0, 1.02)
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    plt.show()

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

# %%
if de_dir:
    de_dir = Path(de_dir)
    per_study_files = sorted((de_dir / "per_study").glob("*.csv"))
    merged_file = de_dir / "merged.csv"
    merged_results = pd.read_csv(merged_file) if merged_file.is_file() else None
    per_study_results = []
    for result_path in per_study_files:
        result = pd.read_csv(result_path)
        per_study_results.append(result)
    has_per_study = bool(per_study_results)
    has_merged = merged_results is not None
    if has_per_study or has_merged or de_run_metadata:
        display(Markdown("## Pseudobulk differential expression"))
        if de_run_metadata.get("analysis_mode") == "test_only":
            display(Markdown(
                "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
                + de_run_metadata["interpretation_warning"]
            ))
        if not has_per_study and not has_merged:
            display(Markdown(
                "_No estimable PyDESeq2 result files were produced for this cell type._"
            ))
    if has_per_study or has_merged:
        de_alpha = pipeline["differential_expression"]["alpha"]
        per_study_all = pd.concat(per_study_results, ignore_index=True) if has_per_study else pd.DataFrame()
        if has_per_study:
            gene_universes = [set(result["gene"].astype(str)) for result in per_study_results]
            common_per_study_genes = set.intersection(*gene_universes)
            per_study_significant = per_study_all.loc[
                pd.to_numeric(per_study_all["padj"], errors="coerce") < de_alpha
            ].copy()
            per_study_hits = per_study_significant.loc[
                per_study_significant["gene"].astype(str).isin(common_per_study_genes)
            ].copy()
            recurrence = (
                per_study_hits.groupby("gene", as_index=False)
                .agg(
                    studies_associated=("study", "nunique"),
                    studies=("study", lambda values: ", ".join(sorted(set(values)))),
                    maximum_absolute_log2_fold_change=(
                        "log2FoldChange", lambda values: pd.to_numeric(values, errors="coerce").abs().max()
                    ),
                )
                .sort_values(
                    ["studies_associated", "maximum_absolute_log2_fold_change", "gene"],
                    ascending=[False, False, True],
                )
            )
        else:
            recurrence = pd.DataFrame()

        summary_rows = []
        if has_per_study:
            significant_by_study = (
                per_study_significant.groupby("study")["gene"].nunique().sort_index()
            )
            summary_rows.extend(
                {"model": f"Per study: {study}", "FDR-significant age-associated genes": int(count)}
                for study, count in significant_by_study.items()
            )
            summary_rows.append(
                {
                    "model": "Per-study union",
                    "FDR-significant age-associated genes": int(
                        per_study_significant["gene"].nunique()
                    ),
                }
            )
            summary_rows.append(
                {
                    "model": "Per-study common-universe union",
                    "FDR-significant age-associated genes": int(per_study_hits["gene"].nunique()),
                }
            )
        if has_merged:
            merged_padj = pd.to_numeric(merged_results["padj"], errors="coerce")
            summary_rows.append(
                {
                    "model": "Merged shared-slope fit",
                    "FDR-significant age-associated genes": int((merged_padj < de_alpha).sum()),
                }
            )
        display(Markdown(
            f"PyDESeq2 age effects; FDR threshold **{de_alpha:g}**. "
            "Gene counts use adjusted p-values in each fitted model. Study recurrence "
            "is restricted to genes tested in every available per-study result."
        ))
        display(pd.DataFrame(summary_rows))

        if has_per_study:
            figure, axes = plt.subplots(1, 2, figsize=(13, 5))
            if recurrence.empty:
                empty_message = (
                    "No genes were tested in every per-study model"
                    if not common_per_study_genes
                    else "No common-universe genes pass the FDR threshold"
                )
                axes[0].text(0.5, 0.5, empty_message,
                             transform=axes[0].transAxes, ha="center", va="center")
                axes[1].text(0.5, 0.5, empty_message,
                             transform=axes[1].transAxes, ha="center", va="center")
            else:
                top_genes = recurrence.head(20).sort_values("studies_associated")
                sns.barplot(
                    data=top_genes, x="studies_associated", y="gene", color="#4c72b0",
                    ax=axes[0],
                )
                axes[0].set_xlabel("Studies with FDR-significant age association")
                axes[0].set_ylabel("")
                axes[0].set_title("Genes recurring across studies (top 20)")
                recurrence_counts = recurrence["studies_associated"].value_counts().sort_index()
                sns.barplot(
                    x=recurrence_counts.index.astype(str), y=recurrence_counts.values,
                    color="#55a868", ax=axes[1],
                )
                axes[1].set_xlabel("Number of studies")
                axes[1].set_ylabel("Unique associated genes")
                axes[1].set_title(
                    f"Common-universe union: {recurrence['gene'].nunique():,} genes"
                )
            figure.tight_layout()
            plt.show()
            if not recurrence.empty:
                display(recurrence.head(20))

        if has_merged:
            volcano = merged_results.copy()
            volcano["padj"] = pd.to_numeric(volcano["padj"], errors="coerce")
            volcano["log2FoldChange"] = pd.to_numeric(
                volcano["log2FoldChange"], errors="coerce"
            )
            volcano = volcano.dropna(subset=["padj", "log2FoldChange"]).copy()
            if volcano.empty:
                display(Markdown("_The merged model has no finite adjusted p-values to plot._"))
            else:
                volcano["minus_log10_padj"] = -np.log10(
                    volcano["padj"].clip(lower=np.finfo(float).tiny)
                )
                volcano["association"] = np.where(
                    volcano["padj"] < de_alpha, "FDR significant", "Not significant"
                )
                figure, axis = plt.subplots(figsize=(8, 6))
                sns.scatterplot(
                    data=volcano, x="log2FoldChange", y="minus_log10_padj",
                    hue="association", hue_order=["Not significant", "FDR significant"],
                    palette={"Not significant": "#b8b8b8", "FDR significant": "#c44e52"},
                    s=16, alpha=0.7, linewidth=0, ax=axis,
                )
                axis.axhline(-np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1)
                axis.axvline(0, color="#555555", linewidth=0.8)
                axis.set_xlabel("Age effect (log2 fold change per year)")
                axis.set_ylabel("−log10(adjusted p-value)")
                axis.set_title("Merged shared-slope age association")
                labels = volcano.loc[volcano["padj"] < de_alpha].nsmallest(12, "padj")
                for _, row in labels.iterrows():
                    axis.annotate(
                        str(row["gene"]),
                        (row["log2FoldChange"], row["minus_log10_padj"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7,
                    )
                figure.tight_layout()
                plt.show()

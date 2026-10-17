# %%
import base64
import os
import re
import warnings
from html import escape
from io import BytesIO
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
from IPython.display import HTML, Markdown, display

from pbmc_pipeline.report_theme import render_report_header

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 110})
from pbmc_pipeline.cell_type import add_unintegrated_umap
from pbmc_pipeline.config import read_json
from pbmc_pipeline.reporting import sample_cluster_fractions

primitive = sc.read_h5ad(Path(os.environ["CELL_TYPE_H5AD"]))
analysis_path = Path(os.environ["CELL_TYPE_ANALYSIS_H5AD"])
output_dir = analysis_path.parent
pipeline = read_json(Path(os.environ["CELL_TYPE_CONFIG"]))
sc.settings.verbosity = 0
report = read_json(output_dir / "report.json")
adata = sc.read_h5ad(analysis_path)
cell_type_name = str(primitive.obs["aifi_l2_majority"].iloc[0])
status_label = str(report.get("status", "unknown")).replace("_", " ").title()
hero_samples = (
    primitive.obs.groupby(["study", "sample"], observed=True)
    .agg(age=("age", "first"))
    .reset_index()
)
hero_ages = pd.to_numeric(hero_samples["age"], errors="coerce").dropna()
hero_age_range = (
    f"{hero_ages.min():g}–{hero_ages.max():g} years" if not hero_ages.empty else "Not provided"
)
summary_metrics = [
    ("Cells", f"{primitive.n_obs:,}"),
    ("Studies", f"{hero_samples['study'].nunique():,}"),
    ("Retained samples", f"{len(hero_samples):,}"),
    ("Age support", hero_age_range),
    ("Local clusters", str(report.get("n_clusters", "—"))),
    ("Analysis status", status_label),
]
display(HTML(render_report_header(
    title=cell_type_name,
    eyebrow="PBMC ageing · per-cell-type report",
    subtitle=(
        "Descriptive sample-level and single-cell diagnostics. "
        "Cell-level associations are not independent-sample inference."
    ),
    metrics=summary_metrics,
)))
provenance = pd.Series(report, name="value").to_frame()
display(HTML(
    "<details class=\"report-provenance\"><summary>Analysis provenance</summary>"
    + provenance.to_html(escape=True)
    + "</details>"
))


def display_collapsible_table(data, summary, *, index=False, max_height=None):
    """Render a table inside a closed-by-default details panel."""
    table_html = data.to_html(index=index, escape=True, border=0)
    if max_height is not None:
        table_html = (
            f'<div style="max-height: {int(max_height)}px; overflow: auto;" '
            'tabindex="0">' + table_html + "</div>"
        )
    display(HTML(
        f'<details class="report-details"><summary>{escape(summary)}</summary>'
        f"{table_html}</details>"
    ))


def display_collapsible_figure(figure, summary, *, open_by_default=False):
    """Render a figure in a details panel, optionally open on load."""
    image_buffer = BytesIO()
    figure.savefig(image_buffer, format="png", bbox_inches="tight")
    image_data = base64.b64encode(image_buffer.getvalue()).decode("ascii")
    plt.close(figure)
    escaped_summary = escape(summary, quote=True)
    display(HTML(
        '<details class="report-details report-details-figure"'
        + (" open" if open_by_default else "")
        + ">"
        f"<summary>{escape(summary)}</summary>"
        f'<img alt="{escaped_summary}" src="data:image/png;base64,{image_data}">'
        "</details>"
    ))


def natural_sort_key(value):
    """Sort labels naturally so numbered clusters follow numeric order."""
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(value))
    )


def plot_umap(
    adata, *, color, title=None, label_clusters=False, marker_area_scale=1.0,
    basis="umap",
):
    """Render comparable low-resolution, equal-aspect UMAP panels."""
    # Scanpy's size controls marker area, so 4x area gives about 2x diameter.
    point_size = max(2, 100_000 / adata.n_obs) * marker_area_scale
    plot_kwargs = {"color": color, "size": point_size, "show": False}
    if label_clusters:
        plot_kwargs["legend_loc"] = "on data"
    sc.pl.embedding(adata, basis=basis, **plot_kwargs)
    figure = plt.gcf()
    figure.set_size_inches(9.5, 6.75)
    figure.set_dpi(110)
    for axis in figure.axes:
        if axis.has_data():
            axis.set_aspect("equal", adjustable="box")
    if title:
        figure.axes[0].set_title(title)
    figure.tight_layout()
    plt.show()


def plot_umap_split_by_study(adata):
    """Show study subsets in their shared local UMAP coordinate system."""
    if "study" not in adata.obs or "X_umap" not in adata.obsm:
        display(Markdown("_Study-split UMAP requires study labels and UMAP coordinates._"))
        return
    coordinates = np.asarray(adata.obsm["X_umap"])
    studies = adata.obs["study"].astype(str)
    study_order = sorted(studies.unique())
    n_columns = min(4, len(study_order))
    n_rows = int(np.ceil(len(study_order) / n_columns))
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(4.3 * n_columns, 3.8 * n_rows), squeeze=False
    )
    point_size = max(2, 100_000 / max(adata.n_obs, 1)) * 4.0
    x_padding = max(np.ptp(coordinates[:, 0]) * 0.03, 0.1)
    y_padding = max(np.ptp(coordinates[:, 1]) * 0.03, 0.1)
    for axis, study in zip(axes.flat, study_order):
        selected = studies.to_numpy() == study
        axis.scatter(
            coordinates[selected, 0], coordinates[selected, 1],
            s=point_size, color="#4c72b0", alpha=0.55, linewidths=0,
            rasterized=True,
        )
        axis.set_title(f"{study} (n={int(selected.sum()):,})")
        axis.set_xlim(
            float(coordinates[:, 0].min() - x_padding),
            float(coordinates[:, 0].max() + x_padding),
        )
        axis.set_ylim(
            float(coordinates[:, 1].min() - y_padding),
            float(coordinates[:, 1].max() + y_padding),
        )
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
    for axis in axes.flat[len(study_order):]:
        axis.set_axis_off()
    figure.suptitle("Local UMAP split by study")
    figure.tight_layout()
    display_collapsible_figure(figure, "Show local UMAP split by study")


def has_age_span(data, *, min_samples=2):
    return len(data) >= min_samples and data["age"].nunique() >= 2


def facet_grid_shape(n_panels, *, max_columns=4):
    """Keep study panels compact while growing to any number of studies."""
    n_columns = min(max_columns, n_panels)
    return (n_panels + n_columns - 1) // n_columns, n_columns


def plot_cluster_fractions(data, *, palette):
    """Keep study identity in points while showing one uncluttered pooled trend."""
    if data.empty:
        display(Markdown("_No samples meet the denominator cutoff for cluster fractions._"))
        return
    plot_data = data.assign(study=data["study"].astype(str))
    clusters = sorted(plot_data["cluster"].astype(str).unique(), key=natural_sort_key)
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(5.5 * n_columns, 4.4 * n_rows), squeeze=False,
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
    display_collapsible_figure(figure, "Show cluster-composition curves")


def directional_cluster_genes(markers, *, cluster, n_genes=8):
    """Return the strongest positive and negative genes for one cluster contrast."""
    subset = markers.loc[markers["group"].astype(str) == str(cluster)].copy()
    positive = subset.loc[subset["logfoldchanges"] > 0].nlargest(n_genes, "logfoldchanges")
    negative = subset.loc[subset["logfoldchanges"] < 0].nsmallest(n_genes, "logfoldchanges")
    result = pd.concat([positive, negative]).copy()
    result["direction"] = result["logfoldchanges"].ge(0).map({
        True: "Higher in cluster", False: "Lower in cluster",
    })
    return result.reset_index(drop=True)


def plot_pca_study_technology_and_introns(adata):
    """Show native PC scores with detailed and collapsed protocol covariates."""
    scores = adata.obsm["X_pca"]
    if scores.shape[1] < 2:
        display(Markdown("_Fewer than two local PCs are available to plot._"))
        return

    pairs = [(0, 1)]
    if scores.shape[1] >= 4:
        pairs.append((2, 3))
    score_data = adata.obs[["study", "technology", "include_intronic"]].copy()
    score_data["study"] = score_data["study"].astype(str)
    score_data["technology"] = score_data["technology"].astype(str)
    score_data["include_intronic"] = score_data["include_intronic"].astype(str)

    def technology_end(value):
        normalized = value.lower().replace("′", "'").replace(" ", "")
        if "3'" in normalized:
            return "3′"
        if "5'" in normalized:
            return "5′"
        return "other"

    intronic_labels = {
        "yes": "Yes",
        "no": "No",
        "not_provided": "Not provided",
    }
    score_data["technology_end"] = score_data["technology"].map(technology_end)
    score_data["include_intronic_label"] = score_data["include_intronic"].map(
        lambda value: intronic_labels.get(value.lower(), value)
    )
    for component in {component for pair in pairs for component in pair}:
        score_data[f"PC {component + 1}"] = scores[:, component]

    studies = sorted(score_data["study"].unique())
    technologies = sorted(score_data["technology"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    intronic_levels = [
        label for label in ["Yes", "No", "Not provided"]
        if label in score_data["include_intronic_label"].unique()
    ]
    intronic_palette = {
        level: {"Yes": "#0072B2", "No": "#D55E00", "Not provided": "#7F7F7F"}[level]
        for level in intronic_levels
    }
    marker_cycle = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "h"]
    technology_markers = {
        technology: marker_cycle[index % len(marker_cycle)]
        for index, technology in enumerate(technologies)
    }
    technology_end_order = [
        end for end in ["3′", "5′", "other"] if end in score_data["technology_end"].unique()
    ]
    technology_end_markers = {"3′": "o", "5′": "s", "other": "X"}
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])

    figure, axes = plt.subplots(2, len(pairs), figsize=(9 * len(pairs), 12), squeeze=False)
    for axis, (x_component, y_component) in zip(axes[0], pairs):
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
        axis.set_title(f"PC {x_component + 1} vs. PC {y_component + 1}: study and technology")
        axis.set_xlabel(f"PC {x_component + 1} ({variance_ratio[x_component]:.1%} variance)")
        axis.set_ylabel(f"PC {y_component + 1} ({variance_ratio[y_component]:.1%} variance)")
    for axis, (x_component, y_component) in zip(axes[1], pairs):
        sns.scatterplot(
            data=score_data,
            x=f"PC {x_component + 1}",
            y=f"PC {y_component + 1}",
            hue="include_intronic_label",
            hue_order=intronic_levels,
            style="technology_end",
            style_order=technology_end_order,
            palette=intronic_palette,
            markers=technology_end_markers,
            s=max(8, 120_000 / adata.n_obs),
            alpha=0.65,
            edgecolor="none",
            legend=False,
            rasterized=True,
            ax=axis,
        )
        axis.set_title(f"PC {x_component + 1} vs. PC {y_component + 1}: intronic reads and 3′/5′")
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
    intronic_handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=intronic_palette[level], label=level, markersize=6)
        for level in intronic_levels
    ]
    technology_end_handles = [
        plt.Line2D(
            [], [], marker=technology_end_markers[end], linestyle="", color="0.35",
            label=end, markersize=6,
        )
        for end in technology_end_order
    ]
    study_legend = figure.legend(
        handles=study_handles, title="Study", loc="upper left",
        bbox_to_anchor=(0.81, 0.94),
    )
    figure.add_artist(study_legend)
    technology_legend = figure.legend(
        handles=technology_handles, title="Technology", loc="upper left",
        bbox_to_anchor=(0.81, 0.67),
    )
    figure.add_artist(technology_legend)
    intronic_legend = figure.legend(
        handles=intronic_handles, title="Intronic reads used\nin alignment", loc="upper left",
        bbox_to_anchor=(0.81, 0.39),
    )
    figure.add_artist(intronic_legend)
    figure.legend(
        handles=technology_end_handles, title="Technology end", loc="upper left",
        bbox_to_anchor=(0.81, 0.17),
    )
    figure.suptitle("Native local PCA scores by study and protocol covariates", y=0.99)
    figure.tight_layout(rect=(0, 0, 0.8, 0.96))
    display_collapsible_figure(figure, "Show PCA scores by study and protocol")


def plot_pca_scores_on_umap(adata):
    """Show the first four native PC scores over the Harmony-derived UMAP."""
    if "X_umap" not in adata.obsm or adata.obsm["X_pca"].shape[1] < 4:
        display(Markdown("_PC scores on UMAP require UMAP coordinates and at least four PCs._"))
        return
    coordinates = np.asarray(adata.obsm["X_umap"])
    scores = np.asarray(adata.obsm["X_pca"])
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
    figure, axes = plt.subplots(2, 2, figsize=(12, 10))
    for component, axis in enumerate(axes.flat[:4]):
        points = axis.scatter(
            coordinates[:, 0], coordinates[:, 1], c=scores[:, component],
            cmap="viridis", s=max(5, 80_000 / adata.n_obs), alpha=0.8,
            linewidths=0, rasterized=True,
        )
        axis.set_title(f"PC {component + 1} ({variance_ratio[component]:.1%} variance)")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        figure.colorbar(points, ax=axis, label=f"PC {component + 1} score")
    figure.suptitle("Native PCA scores over the Harmony-derived UMAP")
    figure.tight_layout()
    display_collapsible_figure(figure, "Show PC 1–4 scores on UMAP")

# %% [markdown]
# <a id="sample-coverage"></a>
# ## Sample Coverage Among Retained Samples
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
        include_intronic=("include_intronic", "first"),
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
        include_intronic=("include_intronic", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        disease_status=("disease_status", lambda values: ", ".join(sorted(values.astype(str).unique()))),
    )
    .reset_index()
)
display_collapsible_table(study_coverage, "Show study coverage and protocol details")
figure, axis = plt.subplots(figsize=(10, 4))
sns.stripplot(data=sample_covariates, x="study", y="age", hue="sex", dodge=True, ax=axis)
axis.set_title("Age coverage by study and sex")
axis.tick_params(axis="x", rotation=45)
figure.tight_layout()
plt.show()

# %% [markdown]
# <a id="local-embedding"></a>
# ## Local Embedding and Clusters
#
# This is a type-specific exploratory embedding. PCA is calculated from the
# local highly variable genes after normalization, log transformation, and
# scaling; configured V(D)J genes are excluded. When more than one study is
# represented, Harmony corrects these PCs by study. The configured local
# neighbor graph, UMAP, and Leiden clusters are calculated from the
# Harmony-adjusted PCs; native PCs are retained below for diagnostic plots.
# Integration can also change biological structure, so study colouring and
# sample-level cluster composition should be reviewed before interpreting a
# cluster as a biological state.

# %%
if report["status"] != "complete":
    display(Markdown("_Too few cells for a type-specific embedding and clustering._"))
else:
    unintegrated_basis = add_unintegrated_umap(adata, pipeline)
    plot_umap(
        adata, color="study", basis=unintegrated_basis.removeprefix("X_"),
        title="Local UMAP without study integration",
        marker_area_scale=4.0,
    )
    plot_umap(adata, color="cluster", title="Local UMAP: cluster IDs", label_clusters=True)
    plot_umap_split_by_study(adata)
    colors = ["study", "age"]
    if "aifi_l3_majority" in adata.obs:
        colors.append("aifi_l3_majority")
    for color in colors:
        plot_umap(adata, color=color, marker_area_scale=4.0)

# %% [markdown]
# <a id="cluster-markers"></a>
# ### Cluster Markers
#
# Genes are ranked by a Wilcoxon comparison of each local cluster against all
# other local clusters. Signed log fold changes therefore show whether a gene
# is higher or lower in that cluster; these are descriptive contrasts.

# %%
if report["status"] == "complete":
    markers = pd.read_csv(output_dir / "markers.tsv", sep="\t")
    clusters = sorted(markers["group"].astype(str).unique(), key=natural_sort_key)
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5.5 * n_columns, 4.4 * n_rows), squeeze=False)
    for axis, cluster in zip(axes.flat, clusters):
        subset = directional_cluster_genes(markers, cluster=cluster)
        sns.barplot(
            data=subset, x="logfoldchanges", y="names", hue="direction", dodge=False,
            order=subset["names"].astype(str).tolist(),
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
    display_collapsible_figure(
        figure, f"Show marker genes across {len(clusters)} local clusters", open_by_default=True,
    )

# %% [markdown]
# <a id="cluster-composition"></a>
# ### Cluster Composition
#
# The plots use one study × sample fraction per point, so they describe
# sample-level cluster composition rather than treating cells as independent
# observations. Curves are minimized on load.

# %%
if report["status"] == "complete":
    min_fraction_denominator = pipeline["cell_type_analysis"][
        "min_cluster_fraction_denominator_cells"
    ]
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
    display(Markdown(
        "_Each point is a sample-level cluster fraction. Differences among coloured study series "
        "are descriptive and should be read alongside the coverage and protocol context above._"
    ))

# %% [markdown]
# <a id="pca-diagnostics"></a>
# ### PCA Diagnostics
#
# These are the native, local PCA scores before Harmony adjustment. Colour uses
# the report's standard study palette and marker shape indicates the recorded
# single-cell technology in the top row. The lower row collapses technology to
# 3′ or 5′ marker shapes and maps intronic read inclusion to colour. The
# Harmony-adjusted representation is used only for the local neighbour graph,
# UMAP, and clustering.

# %%
if report["status"] == "complete":
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
    figure, axis = plt.subplots(figsize=(9, 4.5))
    n_components_to_plot = min(20, len(variance_ratio))
    components = np.arange(1, n_components_to_plot + 1)
    axis.bar(components, variance_ratio[:n_components_to_plot] * 100, color="#4c72b0")
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Variance explained (%)")
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    plt.show()

# %%
if report["status"] == "complete":
    plot_pca_study_technology_and_introns(adata)
    plot_pca_scores_on_umap(adata)

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
    display_collapsible_figure(figure, f"Show gene loadings for the first {n_components_to_plot} PCs")

# %%
if report["status"] == "complete":
    pc_age = pd.read_csv(output_dir / "pc_age_correlations.tsv", sep="\t")
    figure, axis = plt.subplots(figsize=(9, 4.5))
    sns.barplot(data=pc_age, x="pc", y="spearman_r", color="#4c72b0", ax=axis)
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_ylabel("Spearman correlation with age (cells; descriptive)")
    axis.tick_params(axis="x", rotation=45)
    figure.tight_layout()
    display_collapsible_figure(figure, "Show PC–age correlations")

# %%

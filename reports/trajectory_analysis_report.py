# %% [markdown]
# # Cross-cell-type age trajectories
#
# This report clusters standardized, omnibus-significant gene trajectories on
# shared age bins. It also summarizes how often each gene is DE-significant
# across cell types using each cell type's own omnibus FDR result.

# %%
import base64
import json
import os
import warnings
from dataclasses import fields
from html import escape
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display
from upsetplot import UpSet, from_indicators

from pbmc_pipeline.config import AgeTrajectorySettings
from pbmc_pipeline.reporting import plot_trajectory_merge_diagnostic
from pbmc_pipeline.trajectory_analysis import (
    summarize_residual_clustering,
    trajectory_merge_diagnostics,
)

sns.set_theme(style="whitegrid")
analysis_dir = Path(os.environ["TRAJECTORY_ANALYSIS_DIR"])
metadata = json.loads((analysis_dir / "analysis_metadata.json").read_text())
age_trajectory_settings = AgeTrajectorySettings.from_mapping(
    {
        field.name: metadata[field.name]
        for field in fields(AgeTrajectorySettings) if field.init
    }
)
cross_trajectories = pd.read_csv(analysis_dir / "cross_cell_type_trajectory_clusters.csv")
cluster_means = pd.read_csv(analysis_dir / "cross_cell_type_cluster_means.csv")
gene_recurrence = pd.read_csv(analysis_dir / "gene_recurrence.csv")
gene_cell_type_significance = pd.read_csv(analysis_dir / "gene_cell_type_significance.csv")
pattern_concordance = pd.read_csv(analysis_dir / "gene_pattern_concordance.csv")
cell_type_support = pd.DataFrame(metadata.get("cell_type_support", []))
common_bins = metadata.get("common_bins", [])
cluster_fdr = metadata["cluster_fdr_threshold"]
de_fdr = metadata["de_fdr_threshold"]
distance_metric = metadata["distance_metric"]
linkage_method = metadata["linkage_method"]
minimum_umap_trajectories = metadata["minimum_umap_trajectories"]
cross_report_top_n_genes = metadata["cross_report_top_n_genes"]
max_clusters = metadata["max_clusters"]
residual_support = metadata.get("pearson_residual_clustering", {})
cross_residual_support = residual_support.get("cross_cell_type", {})
expression_atlas_support = metadata.get("expression_atlas", {})
atlas_cpm_pseudocount = expression_atlas_support.get("cpm_pseudocount", 1.0)
atlas_log_cpm_label = f"log2(CPM + {atlas_cpm_pseudocount:g})"
expression_matrix_path = analysis_dir / "expression_atlas_matrix.csv"
expression_gene_umap_path = analysis_dir / "expression_atlas_gene_umap.csv"
expression_clusters_path = analysis_dir / "expression_atlas_gene_clusters.csv"
expression_cluster_means_path = analysis_dir / "expression_atlas_cluster_means.csv"
expression_technology_contrast_path = analysis_dir / "expression_atlas_technology_contrast.csv"
expression_technology_summary_path = analysis_dir / "expression_atlas_technology_summary.csv"
expression_intronic_contrast_path = analysis_dir / "expression_atlas_intronic_contrast.csv"
expression_intronic_summary_path = analysis_dir / "expression_atlas_intronic_summary.csv"
expression_gene_availability_path = analysis_dir / "expression_atlas_gene_availability.csv"
expression_gene_intersection_accounting_path = (
    analysis_dir / "expression_atlas_gene_intersection_accounting.csv"
)


def read_optional_csv(path, **kwargs):
    if not path.is_file() or not path.read_text().strip():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


expression_matrix = read_optional_csv(expression_matrix_path, index_col="gene")
expression_gene_umap = read_optional_csv(expression_gene_umap_path)
expression_clusters = read_optional_csv(expression_clusters_path)
expression_cluster_means = read_optional_csv(expression_cluster_means_path)
expression_technology_contrast = read_optional_csv(expression_technology_contrast_path)
expression_technology_summary = read_optional_csv(expression_technology_summary_path)
expression_intronic_contrast = read_optional_csv(expression_intronic_contrast_path)
expression_intronic_summary = read_optional_csv(expression_intronic_summary_path)
expression_gene_availability = read_optional_csv(
    expression_gene_availability_path, index_col="gene"
).astype(bool)
expression_gene_intersection_accounting = read_optional_csv(
    expression_gene_intersection_accounting_path
)

display(HTML("""
<style>
:root {
  --report-ink: #182230;
  --report-muted: #5f6b7a;
  --report-line: #d8e0e8;
  --report-paper: #ffffff;
  --report-page: #f5f7fa;
  --report-accent: #176b87;
  --report-accent-soft: #e6f3f7;
}
body.jp-Notebook { background: var(--report-page); color: var(--report-ink); }
main { max-width: 1600px; margin: 0 auto; padding: 2rem clamp(1.25rem, 3vw, 2.5rem) 4rem; }
.jp-Cell-outputWrapper, .jp-OutputArea { width: 100%; }
.jp-RenderedImage img { display: block; max-width: 100%; height: auto; margin: 0 auto; }
.report-toc-shell { position: fixed; z-index: 1000; top: 1rem; bottom: 1rem; left: 1rem; width: 260px; }
.report-toc-toggle { position: absolute; width: 1px; height: 1px; margin: -1px; overflow: hidden; clip: rect(0, 0, 0, 0); white-space: nowrap; clip-path: inset(50%); }
.report-toc-toggle-label { display: none; }
.report-toc { height: 100%; overflow-y: auto; padding: .9rem .8rem 1.1rem; border: 1px solid var(--report-line); border-radius: 14px; background: var(--report-paper); box-shadow: 0 10px 28px rgba(24, 34, 48, .12); }
.report-toc-heading { margin: 0 0 .65rem; padding: .25rem .45rem .65rem; border-bottom: 1px solid var(--report-line); color: var(--report-muted); font-size: .72rem; font-weight: 750; letter-spacing: .08em; text-transform: uppercase; }
.report-toc-list, .report-toc-list ul { margin: 0; padding: 0; list-style: none; }
.report-toc-list > li { margin: .12rem 0 .45rem; }
.report-toc-list ul { margin: .2rem 0 .35rem .65rem; padding-left: .65rem; border-left: 1px solid var(--report-line); }
.report-toc-list ul li { margin: .12rem 0; }
.report-toc a, .report-toc a:visited { display: block; padding: .25rem .45rem; border-radius: 6px; color: #344454; font-size: .88rem; font-weight: 600; line-height: 1.35; text-decoration: none; }
.report-toc-list ul a { color: var(--report-muted); font-size: .82rem; font-weight: 500; }
.report-toc a:hover, .report-toc a:focus-visible { background: var(--report-accent-soft); color: var(--report-accent); outline: none; }
html { scroll-behavior: smooth; scroll-padding-top: 1.25rem; }
@media print {
  body.jp-Notebook { background: #fff; }
  main { max-width: none; padding: 0; }
  .report-toc-shell { display: none; }
}
@media (min-width: 1200px) {
  main { width: calc(100% - 320px); max-width: 1500px; margin: 0 1.5rem 0 295px; }
}
@media (max-width: 1199px) {
  main { padding-top: 5rem; }
  .report-toc-shell { top: .5rem; right: .5rem; bottom: auto; left: .5rem; width: auto; }
  .report-toc-toggle-label { display: flex; align-items: center; justify-content: space-between; min-height: 2.8rem; padding: .6rem .9rem; border: 1px solid var(--report-line); border-radius: 10px; background: var(--report-paper); box-shadow: 0 6px 18px rgba(24, 34, 48, .12); color: var(--report-accent); cursor: pointer; font-weight: 700; }
  .report-toc-toggle-label::after { content: "＋"; font-size: 1.2rem; }
  .report-toc-toggle:checked + .report-toc-toggle-label::after { content: "−"; }
  .report-toc { display: none; height: auto; max-height: min(72vh, 680px); margin-top: .35rem; border-radius: 10px; }
  .report-toc-toggle:checked ~ .report-toc { display: block; }
  .report-toc-list { columns: 2; column-gap: 1rem; }
  .report-toc-list > li { break-inside: avoid; }
}
@media (max-width: 560px) { .report-toc-list { columns: 1; } }
</style>
<div class="report-toc-shell">
  <input class="report-toc-toggle" type="checkbox" id="cross-report-toc-toggle">
  <label class="report-toc-toggle-label" for="cross-report-toc-toggle">Contents</label>
  <nav class="report-toc" aria-label="Cross-cell-type report table of contents">
    <p class="report-toc-heading">On this page</p>
    <ul class="report-toc-list">
      <li><a href="#Support-and-scope">Support and scope</a></li>
      <li><a href="#Study-balanced-cell-type-expression-atlas">Expression atlas</a>
        <ul>
          <li><a href="#Gene-intersection-across-studies">Gene intersection</a></li>
          <li><a href="#AIFI-L1-mean-scaled-expression">AIFI L1 scaled expression</a></li>
          <li><a href="#3%E2%80%B2/5%E2%80%B2-technology-associated-contrast">3′/5′ contrast</a></li>
          <li><a href="#Intronic-read-inclusion-contrast">Intronic contrast</a></li>
        </ul>
      </li>
      <li><a href="#Sample-clustering-from-bins-model-Pearson-residuals">Residual clustering</a>
        <ul><li><a href="#Cluster-focused-residual-characterization">Cluster characterization</a></li></ul>
      </li>
      <li><a href="#Age-bin-DE-genes-recurring-in-multiple-cell-types">Recurring age-bin DE genes</a></li>
      <li><a href="#Cross-cell-type-trajectory-clusters">Trajectory clusters</a></li>
      <li><a href="#Do-recurring-genes-share-trajectory-shapes?">Trajectory-shape concordance</a></li>
    </ul>
  </nav>
</div>
"""))

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
    f"Age-trajectory fits exclude the 90–100 bin and apply a strict age cutoff "
    f"of < {metadata['strict_age_cutoff_exclusive']} years; each retained bin has at least "
    f"{metadata['minimum_samples_per_bin']} eligible pseudobulks."
))
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

# %% [markdown]
# ### Gene intersection across studies

# %%
display(Markdown(
    "**Gene universes in this workflow**\n\n"
    "| Context | Gene universe | Purpose |\n"
    "|---|---|---|\n"
    "| Merged pseudobulk | Outer union across studies | Retains observed genes for per-study models. |\n"
    "| Per-study DE | Genes available in that study | Fits each study's observed counts. |\n"
    "| Combined DE | Intersection across the studies in that fit | Excludes study-absent outer-join zeros. |\n"
    "| This expression atlas | Intersection across all studies | Makes study-balanced expression comparable. |"
))
if expression_gene_availability.empty or expression_gene_availability.shape[1] < 2:
    display(Markdown(
        "_Gene-availability overlap requires an outer pseudobulk merge with at least two studies._"
    ))
else:
    n_shared_genes = int(expression_gene_availability.all(axis=1).sum())
    n_union_genes = len(expression_gene_availability)
    study_gene_counts = expression_gene_availability.sum(axis=0).rename("available genes")
    display(Markdown(
        f"The pseudobulk merge retains the outer union of **{n_union_genes:,} genes** so "
        "per-study models can use each study's observed genes. The study-balanced expression "
        f"atlas below deliberately uses the **{n_shared_genes:,} genes** available in every "
        "study. This is the relevant gene-intersection diagnostic; it is shown here rather "
        "than in merge QC because the merge itself does not discard study-specific genes."
    ))
    display(study_gene_counts.to_frame())
    if not expression_gene_intersection_accounting.empty:
        accounting = expression_gene_intersection_accounting.copy()
        total_input_counts = int(accounting["input_counts"].sum())
        total_discarded_counts = int(accounting["discarded_counts"].sum())
        total_discarded_fraction = (
            total_discarded_counts / total_input_counts if total_input_counts else 0.0
        )
        display(Markdown(
            f"Restricting to the shared-gene intersection would exclude **{total_discarded_counts:,} "
            f"of {total_input_counts:,} pseudobulk counts ({100 * total_discarded_fraction:.3f}%)**. "
            "These are hypothetical exclusions for the study-balanced atlas; the outer pseudobulk "
            "merge retains those counts for per-study models."
        ))
        plot_rows = accounting.sort_values("discarded_fraction", ascending=False)
        figure, axis = plt.subplots(figsize=(9, max(3, 0.4 * len(plot_rows))))
        percentages = 100 * plot_rows["discarded_fraction"]
        bars = axis.barh(plot_rows["study"], percentages, color="#c44e52")
        axis.invert_yaxis()
        axis.set_xlabel("Counts outside shared-gene intersection (%)")
        axis.set_ylabel("Study")
        axis.set_title("Pseudobulk counts outside the all-study gene intersection")
        axis.set_xlim(0, max(0.5, float(percentages.max()) * 1.15))
        axis.bar_label(bars, labels=[f"{value:.2f}%" for value in percentages], padding=3)
        figure.tight_layout()
        plt.show()
        accounting["counts outside intersection (%)"] = (
            100 * accounting["discarded_fraction"]
        ).map(lambda value: f"{value:.3f}%")
        accounting["top excluded gene"] = accounting.apply(
            lambda row: (
                f"{row['top_discarded_gene']} ({int(row['top_discarded_gene_counts']):,})"
                if pd.notna(row["top_discarded_gene"])
                else "—"
            ),
            axis=1,
        )
        display(accounting[[
            "study", "input_gene_symbols", "retained_gene_symbols", "discarded_gene_symbols",
            "input_counts", "discarded_counts", "counts outside intersection (%)",
            "top excluded gene",
        ]].rename(columns={
            "input_gene_symbols": "input genes",
            "retained_gene_symbols": "shared genes",
            "discarded_gene_symbols": "genes outside intersection",
            "input_counts": "input counts",
            "discarded_counts": "counts outside intersection",
        }))
    # UpSetPlot 0.9 uses chained assignment internally; suppress only its
    # known pandas-3.0 compatibility warning while rendering this plot.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"A value is trying to be set on a copy.*",
            category=FutureWarning,
            module=r"upsetplot\.plotting",
        )
        axes = UpSet(
            from_indicators(expression_gene_availability), subset_size="count", sort_by="cardinality"
        ).plot()
    figure = axes["intersections"].figure
    figure.suptitle("Gene availability across studies", y=1.02)
    figure.subplots_adjust(top=0.9)
    plt.show()
    display(Markdown(
        "Download the [gene-availability matrix](expression_atlas_gene_availability.csv) "
        "and [gene-intersection count accounting](expression_atlas_gene_intersection_accounting.csv)."
    ))

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
    if not expression_gene_umap.empty:
        expression_umap_plot = expression_gene_umap.merge(
            expression_matrix.mean(axis=1).rename("mean_log2_cpm"), left_on="gene", right_index=True
        )
        if not expression_technology_summary.empty:
            expression_umap_plot = expression_umap_plot.merge(
                expression_technology_summary[[
                    "gene", "median_log2_cpm_difference_5_prime_minus_3_prime"
                ]], on="gene", how="left"
            )
        if not expression_intronic_summary.empty:
            expression_umap_plot = expression_umap_plot.merge(
                expression_intronic_summary[[
                    "gene", "median_log2_cpm_difference_intronic_minus_non_intronic"
                ]], on="gene", how="left"
            )
        color_columns = [
            ("scaling_sd_log2_cpm", "Cross-L2 SD used for scaling"),
            ("mean_log2_cpm", f"Mean {atlas_log_cpm_label}"),
            (
                "median_log2_cpm_difference_5_prime_minus_3_prime",
                f"Cross-type median 5′ − 3′ difference\nin {atlas_log_cpm_label}",
            ),
            (
                "median_log2_cpm_difference_intronic_minus_non_intronic",
                f"Cross-type median intronic − non-intronic difference\nin {atlas_log_cpm_label}",
            ),
        ]
        figure, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")
        for axis, (column, title) in zip(axes.ravel(), color_columns):
            if column not in expression_umap_plot or expression_umap_plot[column].notna().sum() == 0:
                axis.text(0.5, 0.5, "Not available", ha="center", va="center")
                axis.set_axis_off()
                continue
            points = expression_umap_plot.dropna(subset=[column])
            scatter = axis.scatter(
                points["umap_1"], points["umap_2"], c=points[column], s=8,
                cmap="viridis" if column in {"mean_log2_cpm", "scaling_sd_log2_cpm"} else "vlag",
                vmin=None if column in {"mean_log2_cpm", "scaling_sd_log2_cpm"} else -np.nanmax(np.abs(points[column])),
                vmax=None if column in {"mean_log2_cpm", "scaling_sd_log2_cpm"} else np.nanmax(np.abs(points[column])),
                linewidths=0, alpha=0.8,
            )
            figure.colorbar(scatter, ax=axis, label=title)
            axis.set_title(title)
            axis.set_xlabel("Gene-expression UMAP 1")
            axis.set_ylabel("Gene-expression UMAP 2")
        figure.suptitle(
            f"Gene expression UMAP: per-gene standardized only for PCA ({expression_atlas_support.get('n_gene_expression_umap_pcs', 0)} PCs) and UMAP; colors use unscaled values",
            wrap=True,
        )
        plt.show()
    display(Markdown("### AIFI L1 mean scaled expression"))
    l1_parent_by_l2 = expression_atlas_support.get("l1_parent_by_l2", {})
    l1_groups = {}
    for cell_type in expression_matrix.columns:
        l1_parent = l1_parent_by_l2.get(str(cell_type))
        if l1_parent:
            l1_groups.setdefault(str(l1_parent), []).append(str(cell_type))
    if l1_groups and not expression_gene_umap.empty:
        l1_expression = expression_matrix.reindex(expression_gene_umap["gene"]).astype(float)
        scaling_sd = expression_gene_umap.set_index("gene").get("scaling_sd_log2_cpm")
        if scaling_sd is None:
            scaling_sd = l1_expression.std(axis=1, ddof=0)
        scaling_sd = scaling_sd.reindex(l1_expression.index).replace(0, np.nan)
        scaled_expression = l1_expression.sub(l1_expression.mean(axis=1), axis=0).div(
            scaling_sd, axis=0
        )
        l1_values = {
            l1_parent: scaled_expression[cell_types].mean(axis=1)
            for l1_parent, cell_types in sorted(l1_groups.items())
        }
        l1_limit = np.nanmax(np.abs(np.concatenate([
            values.to_numpy(dtype=float) for values in l1_values.values()
        ])))
        l1_limit = l1_limit if np.isfinite(l1_limit) and l1_limit > 0 else 1.0
        n_l1_columns = min(3, len(l1_values))
        n_l1_rows = int(np.ceil(len(l1_values) / n_l1_columns))
        figure, axes = plt.subplots(
            n_l1_rows, n_l1_columns,
            figsize=(5.4 * n_l1_columns, 4.6 * n_l1_rows),
            layout="constrained",
        )
        axes = np.atleast_1d(axes).ravel()
        scatter = None
        for axis, (l1_parent, values) in zip(axes, l1_values.items()):
            scatter = axis.scatter(
                expression_gene_umap["umap_1"], expression_gene_umap["umap_2"],
                c=values.to_numpy(dtype=float), s=8, cmap="vlag",
                vmin=-l1_limit, vmax=l1_limit, linewidths=0, alpha=0.8,
            )
            axis.set_title(f"{l1_parent} ({len(l1_groups[l1_parent])} L2 types)")
            axis.set_xlabel("Gene-expression UMAP 1")
            axis.set_ylabel("Gene-expression UMAP 2")
        for axis in axes[len(l1_values):]:
            axis.set_axis_off()
        if scatter is not None:
            figure.colorbar(
                scatter, ax=axes[:len(l1_values)],
                label="Mean scaled expression across member L2 types",
            )
        plt.show()
    elif l1_parent_by_l2:
        display(Markdown("_No gene-expression UMAP coordinates are available for L1 summaries._"))
    else:
        display(Markdown(
            "_No AIFI L2-to-L1 mapping was recorded with this atlas; L1 summaries are unavailable._"
        ))
    display(Markdown(
        "Download the [study-balanced expression matrix](expression_atlas_matrix.csv) and "
        "[gene-expression UMAP coordinates](expression_atlas_gene_umap.csv)."
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
def plot_labeled_expression_effect(table, *, x, y, title, y_label):
    plotted = table.dropna(subset=[x, y]).copy()
    labels = plotted.assign(
        absolute_effect=lambda values: values[y].abs()
    ).sort_values("absolute_effect", ascending=False, kind="stable").head(
        20
    )
    figure, axis = plt.subplots(figsize=(9, 6))
    axis.scatter(plotted[x], plotted[y], s=12, alpha=0.45, color="#4c72b0", linewidths=0)
    x_margin = max(0.05, 0.06 * (plotted[x].max() - plotted[x].min()))
    y_margin = max(0.03, 0.08 * (plotted[y].max() - plotted[y].min()))
    axis.set_xlim(plotted[x].min() - x_margin, plotted[x].max() + x_margin)
    axis.set_ylim(plotted[y].min() - y_margin, plotted[y].max() + y_margin)
    for _, row in labels.iterrows():
        axis.text(row[x], row[y], str(row["gene"]), fontsize=7, alpha=0.85)
    axis.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    axis.set_xlabel(f"Mean {atlas_log_cpm_label}")
    axis.set_ylabel(y_label)
    axis.set_title(title + f" ({len(plotted):,} genes; top {len(labels)} labeled)")
    figure.tight_layout()
    plt.show()


# %%
if not expression_technology_summary.empty and not expression_technology_contrast.empty:
    display(Markdown(
        "This contrast is computed whenever at least one depth-qualified study exists in "
        "each technology family. The configured two-study threshold is a replication flag, "
        "not an estimability gate."
    ))
    plot_labeled_expression_effect(
        expression_technology_summary,
        x="mean_log2_cpm_across_technologies",
        y="median_log2_cpm_difference_5_prime_minus_3_prime",
        title="Technology-associated expression contrast",
        y_label=f"Cross-type median 5′ − 3′ difference\nin {atlas_log_cpm_label}",
    )
    top_contrast_genes = expression_technology_summary.assign(
        absolute_difference=lambda table: table["median_log2_cpm_difference_5_prime_minus_3_prime"].abs()
    ).sort_values("absolute_difference", ascending=False, kind="stable").head(20)["gene"]
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
            cbar_kws={"label": f"5′ − 3′ difference in {atlas_log_cpm_label}"}, ax=axis,
        )
        axis.set_xlabel("Cell type")
        axis.set_ylabel("Gene")
        axis.set_title("Top genes by cross-type median 5′ − 3′ contrast")
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
    display(Markdown(
        "This independent contrast is `intronic − non-intronic` on separately "
        f"study-balanced {atlas_log_cpm_label} matrices. As for 3′/5′, a one-versus-one "
        "comparison is estimable but is flagged as unreplicated."
    ))
    plot_labeled_expression_effect(
        expression_intronic_summary,
        x="mean_log2_cpm_across_intronic_status",
        y="median_log2_cpm_difference_intronic_minus_non_intronic",
        title="Intronic-read inclusion expression contrast",
        y_label=(
            f"Cross-type median intronic − non-intronic difference\n"
            f"in {atlas_log_cpm_label}"
        ),
    )
    top_intronic_genes = expression_intronic_summary.assign(
        absolute_difference=lambda table: table[
            "median_log2_cpm_difference_intronic_minus_non_intronic"
        ].abs()
    ).sort_values("absolute_difference", ascending=False, kind="stable").head(20)["gene"]
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
            cbar_kws={
                "label": f"Intronic − non-intronic difference in {atlas_log_cpm_label}"
            }, ax=axis,
        )
        axis.set_xlabel("Cell type")
        axis.set_ylabel("Gene")
        axis.set_title("Top genes by cross-type median intronic-read contrast")
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
display(Markdown(
    "Residuals use the age-bin model with sex, log10(total counts), and estimable study-site terms; "
    f"PCA is centered and unscaled (up to {metadata['residual_pca_components']} PCs). "
    f"Residual UMAP uses {metadata['residual_umap_neighbors']} neighbors (capped at samples − 1), "
    f"min_dist={metadata['umap_min_dist']}, Euclidean PC distance, and seed={metadata['random_state']}. "
    f"Primary residual clusters use Leiden on the graph of the first up to "
    f"{metadata['residual_pc_hdbscan_components']} PCs, with "
    f"{metadata['residual_leiden_neighbors']} neighbors. Per-type and combined analyses "
    f"select the lowest tested resolution reaching {metadata['residual_leiden_min_clusters']} "
    "clusters, or the tested resolution with the most clusters if none reach that target. "
    f"Trajectory UMAP separately uses {metadata['umap_neighbors']} neighbors; gene-expression UMAP uses "
    f"{expression_atlas_support.get('umap_neighbors', 'not run')} neighbors after up to "
    f"{expression_atlas_support.get('umap_pca_components', 'not run')} per-gene-standardized PCs."
))

# %%
def display_collapsible_figure(figure, summary, *, open_by_default=False):
    image_buffer = BytesIO()
    figure.savefig(image_buffer, format="png", bbox_inches="tight")
    plt.close(figure)
    image_data = base64.b64encode(image_buffer.getvalue()).decode("ascii")
    display(HTML(
        '<details class="report-details report-details-figure"'
        + (" open" if open_by_default else "")
        + ">"
        f"<summary>{escape(summary)}</summary>"
        f'<img alt="{escape(summary, quote=True)}" '
        f'src="data:image/png;base64,{image_data}"></details>'
    ))


def plot_residual_covariates(
    table, title, covariates, *, categorical_covariates, open_by_default=False
):
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
    display_collapsible_figure(figure, title, open_by_default=open_by_default)


def plot_residual_cluster_alternatives(table, settings):
    summary = summarize_residual_clustering(table, settings)
    shown = summary[[
        "method", "n_clusters", "n_assigned", "n_noise", "noise_fraction", "cluster_sizes",
    ]].copy()
    shown["noise_fraction"] = shown["noise_fraction"].map(
        lambda value: f"{value:.1%}" if pd.notna(value) else "—"
    )
    display(HTML(
        '<details class="report-details"><summary>Residual clustering summary</summary>'
        + shown.to_html(index=False, escape=True, border=0)
        + "</details>"
    ))
    plot_summary = summary.loc[
        summary["label_column"].isin([
            "residual_cluster_kmeans",
            "residual_cluster_hdbscan_umap",
            "residual_cluster_hdbscan_pcs",
        ])
    ]
    coordinates = table.dropna(subset=["umap_1", "umap_2"])
    n_rows = max(1, (len(plot_summary) + 1) // 2)
    figure, axes = plt.subplots(n_rows, 2, figsize=(14, 5.5 * n_rows), squeeze=False)
    for axis, row in zip(axes.flat, plot_summary.to_dict("records")):
        column = row["label_column"]
        values = coordinates.dropna(subset=[column]) if column in coordinates else pd.DataFrame()
        if values.empty:
            axis.text(0.5, 0.5, "UMAP or cluster labels unavailable", ha="center", va="center")
            axis.set_axis_off()
            continue
        labels = pd.to_numeric(values[column], errors="coerce")
        noise = labels.eq(-1)
        if noise.any():
            axis.scatter(
                values.loc[noise, "umap_1"], values.loc[noise, "umap_2"],
                color="#8a8a8a", s=32, alpha=0.75, linewidths=0,
            )
        clusters = sorted(labels.loc[~noise].dropna().astype(int).unique())
        for index, cluster in enumerate(clusters):
            selected = labels.eq(cluster)
            axis.scatter(
                values.loc[selected, "umap_1"], values.loc[selected, "umap_2"],
                color=plt.get_cmap("tab20")(index % 20), s=32, alpha=0.82, linewidths=0,
            )
        axis.set_title(
            f"{row['method']}\n{row['n_clusters']} clusters; {row['n_noise']} noise"
        )
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
    for axis in axes.flat[len(plot_summary):]:
        axis.set_axis_off()
    figure.tight_layout()
    display_collapsible_figure(figure, "Residual cluster method comparison")

cross_residual_path = analysis_dir / "cross_cell_type_residual_clusters.csv"
cross_residual_pc_scores_path = analysis_dir / "cross_cell_type_residual_pc_scores.csv"
cross_residual_markers_path = analysis_dir / "cross_cell_type_residual_cluster_markers.csv"
cross_residual_composition_path = analysis_dir / "cross_cell_type_residual_cluster_feature_composition.csv"
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
        "residual matrix](cross_cell_type_pearson_residuals.npz). "
        f"Primary labels use Leiden at resolution "
        f"{cross_residual_support.get('residual_leiden_selected_resolution', '—')}, selected "
        f"to reach the target of "
        f"{cross_residual_support.get('residual_leiden_min_clusters', '—')} clusters."
    )
    display(Markdown(residual_scope))
    display(Markdown(
        "Residual cluster sizes are encoded in the UMAP legend. Noise samples are excluded from "
        "marker contrasts, which compare each assigned cluster's "
        "mean residual against the largest cluster; feature-origin fractions summarize the top/bottom "
        "500 residual differences by cell type."
    ))
    residual_markers = read_optional_csv(cross_residual_markers_path)
    residual_composition = read_optional_csv(cross_residual_composition_path)
    if cross_residual_support.get("n_samples_common", 0) >= minimum_umap_trajectories:
        residual_cluster_counts = cross_residual_clusters["residual_cluster"].value_counts()
        cross_residual_clusters["residual_cluster_label"] = cross_residual_clusters[
            "residual_cluster"
        ].map(lambda cluster: (
            f"Noise (n={residual_cluster_counts[cluster]})" if int(cluster) == -1
            else f"Cluster {cluster} (n={residual_cluster_counts[cluster]})"
        ))
        figure, axis = plt.subplots(figsize=(6.6, 4.8))
        sns.scatterplot(
            data=cross_residual_clusters, x="umap_1", y="umap_2", hue="residual_cluster_label",
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
                open_by_default=True,
            )
        cross_residual_pc_scores = read_optional_csv(cross_residual_pc_scores_path)
        if {"sample_key", "residual_cluster", "PC1", "PC2"}.issubset(
            cross_residual_pc_scores.columns
        ):
            plot_data = cross_residual_clusters.merge(
                cross_residual_pc_scores[
                    ["sample_key", "residual_cluster", "PC1", "PC2"]
                ],
                on=["sample_key", "residual_cluster"],
                how="inner",
                validate="one_to_one",
            ).dropna(subset=["PC1", "PC2"])
            if not plot_data.empty:
                figure, axes = plt.subplots(1, 3, figsize=(18, 5.2))
                sns.scatterplot(
                    data=plot_data, x="PC1", y="PC2", hue="residual_cluster",
                    palette="tab20", s=34, alpha=0.84, linewidth=0, ax=axes[0],
                )
                axes[0].set_title("Residual PCA: PC1 vs. PC2")
                axes[0].legend(title="Residual cluster", bbox_to_anchor=(1.02, 1))
                for axis, component in zip(axes[1:], ("PC1", "PC2")):
                    values = plot_data.dropna(subset=["umap_1", "umap_2", component])
                    points = axis.scatter(
                        values["umap_1"], values["umap_2"], c=values[component],
                        cmap="viridis", s=34, alpha=0.84, linewidths=0,
                    )
                    axis.set_title(f"Residual UMAP colored by {component}")
                    axis.set_xlabel("UMAP 1")
                    axis.set_ylabel("UMAP 2")
                    figure.colorbar(points, ax=axis, label=f"{component} score")
                figure.tight_layout()
                display_collapsible_figure(figure, "Residual PCA and PC scores over UMAP")
        else:
            display(Markdown("_Cross-cell-type residual PCA scores are unavailable._"))
        display(Markdown("### Alternative clusterings"))
        display(Markdown(
            "Leiden on the PC graph is the primary clustering. It uses the lowest tested "
            f"resolution reaching {age_trajectory_settings.residual_leiden_min_clusters} "
            "clusters, or the tested resolution with the most clusters if none reach that "
            f"target. The summary table includes the tested resolutions "
            f"{age_trajectory_settings.residual_leiden_resolutions}; plots compare K-means "
            "and one configured HDBSCAN result on UMAP and PCs. Marker contrasts use at "
            f"most the largest {age_trajectory_settings.residual_marker_max_clusters} clusters."
        ))
        plot_residual_cluster_alternatives(
            cross_residual_clusters, age_trajectory_settings
        )
        display(Markdown("### Cluster-focused residual characterization"))
        if not residual_markers.empty:
            figure, axis = plt.subplots(figsize=(6.6, 4.8))
            sns.scatterplot(
                data=cross_residual_clusters, x="umap_1", y="umap_2", hue="residual_cluster_label",
                palette="tab20", s=32, alpha=0.84, linewidth=0, ax=axis,
            )
            axis.set_title("Residual clusters characterized below")
            axis.set_xlabel("UMAP 1")
            axis.set_ylabel("UMAP 2")
            if axis.legend_ is not None:
                axis.legend_.set_title("Residual cluster")
                axis.legend_.set_bbox_to_anchor((1.02, 1))
            figure.tight_layout()
            plt.show()
            marker_clusters = sorted(residual_markers["residual_cluster"].unique())
            figure, axes = plt.subplots(
                len(marker_clusters), 2, figsize=(13, 3.4 * len(marker_clusters))
            )
            axes = np.atleast_2d(axes)
            for row_index, cluster in enumerate(marker_clusters):
                for column_index, direction in enumerate(("lower", "higher")):
                    axis = axes[row_index, column_index]
                    markers = residual_markers.loc[
                        (residual_markers["residual_cluster"] == cluster)
                        & (residual_markers["direction"] == direction)
                    ].sort_values("rank").head(10)
                    axis.barh(
                        markers["feature"].iloc[::-1],
                        markers["mean_residual_difference_vs_reference"].iloc[::-1],
                        color="#4c72b0",
                    )
                    axis.set_title(f"Cluster {cluster}: {direction} vs. largest cluster")
                    axis.set_xlabel("Mean Pearson-residual difference")
            figure.tight_layout()
            plt.show()
        if not residual_composition.empty:
            composition = residual_composition.pivot_table(
                index=["residual_cluster", "direction"], columns="cell_type",
                values="fraction_of_top_500_features", fill_value=0,
            )
            figure, axis = plt.subplots(figsize=(9, 4.5))
            composition.plot(kind="bar", stacked=True, ax=axis, colormap="tab20")
            axis.set_xlabel("Residual cluster and feature direction")
            axis.set_ylabel("Fraction of top 500 features")
            axis.set_title("Cell-type contribution to residual-cluster marker features")
            axis.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left")
            figure.tight_layout()
            plt.show()
        elif cross_residual_clusters.loc[
            cross_residual_clusters["residual_cluster"] >= 0, "residual_cluster"
        ].nunique() < 2:
            display(Markdown(
                "_Leiden found fewer than two clusters; there is therefore no "
                "between-cluster marker or feature-origin contrast to report._"
            ))
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
# ## Age-bin DE genes recurring in multiple cell types
#
# Counts use per-cell-type omnibus FDR < 0.05 across completed age-bin fits.

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
        top_recurrent_genes = recurrent.head(12)["gene"]
        effects = gene_cell_type_significance.loc[
            gene_cell_type_significance["gene"].isin(top_recurrent_genes)
            & gene_cell_type_significance["de_significant"]
        ].copy()
        if not effects.empty:
            effect_matrix = effects.pivot(
                index="gene", columns="cell_type", values="largest_absolute_log2_fold_change"
            ).reindex(top_recurrent_genes)
            figure, axis = plt.subplots(
                figsize=(max(8, 1.25 * effect_matrix.shape[1] + 3), max(4, 0.45 * len(effect_matrix) + 2)),
            )
            sns.heatmap(
                effect_matrix, cmap="vlag", center=0, annot=True, fmt=".2g",
                cbar_kws={"label": "Largest absolute age-bin log2 fold change"}, ax=axis,
            )
            axis.set_xlabel("Cell type (only omnibus-significant effects shown)")
            axis.set_ylabel("Gene")
            axis.set_title("Top recurring genes: strongest age-bin coefficient by cell type")
            figure.tight_layout()
            plt.show()
else:
    display(Markdown("_No gene-level cell-type significance results are available._"))

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
    trajectory_settings = AgeTrajectorySettings.from_mapping(
        {
            field.name: metadata[field.name]
            for field in fields(AgeTrajectorySettings) if field.init
        }
    )
    merge_diagnostics = trajectory_merge_diagnostics(
        cross_trajectories.set_index("trajectory_id")[common_bins],
        trajectory_settings,
        maximum_merges=20,
    )
    figure, axes = plt.subplots(1, 3, figsize=(18, 5.2))
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
    plot_trajectory_merge_diagnostic(
        axes[2], merge_diagnostics, max_clusters=max_clusters,
    )
    figure.tight_layout()
    plt.show()
    display(Markdown(
        "The merge-cost panel shows up to the final 20 hierarchical merges. "
        "A sharp increase at **K → K−1** supports retaining K trajectory groups; "
        "the red marker shows the configured cluster cap when that transition is available."
    ))
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
        f"Pairwise shape distances use the configured {distance_metric} metric. The plots below "
        "restore each cell type's trajectory amplitude by multiplying standardized values by its "
        "stored trajectory SD."
    ))
    top_pattern_genes = pattern_concordance.head(12)["gene"].tolist()
    figure = plt.figure(figsize=(21, 11), layout="constrained")
    grid = figure.add_gridspec(3, 5, width_ratios=[1, 1, 1, 1, 0.72])
    axes = np.array([
        figure.add_subplot(grid[row, column])
        for row in range(3) for column in range(4)
    ]).reshape(3, 4)
    legend_axis = figure.add_subplot(grid[:, 4])
    legend_axis.set_axis_off()
    legend_handles = {}
    raw_columns = [f"raw_{age_bin}" for age_bin in common_bins]
    for axis, gene in zip(axes.flat, top_pattern_genes):
        trajectories = cross_trajectories.loc[cross_trajectories["gene"] == gene]
        for _, row in trajectories.iterrows():
            if all(column in trajectories for column in raw_columns):
                values = [row[column] for column in raw_columns]
            else:
                values = [row[age_bin] * row["trajectory_sd"] for age_bin in common_bins]
            line, = axis.plot(
                common_bins, values, marker="o", linewidth=1.6, label=row["cell_type"]
            )
            legend_handles.setdefault(str(row["cell_type"]), line)
        axis.axhline(0, color="#555555", linestyle="--", linewidth=0.6)
        axis.set_title(str(gene), fontsize=9)
        axis.tick_params(axis="x", rotation=35, labelsize=7)
        axis.tick_params(axis="y", labelsize=7)
    for axis in axes.flat[len(top_pattern_genes):]:
        axis.set_axis_off()
    if legend_handles:
        legend_axis.legend(
            legend_handles.values(), legend_handles.keys(), title="Cell type",
            loc="upper left", frameon=True,
        )
    plt.show()

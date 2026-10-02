# %% [markdown]
# # Cross-study merge QC report
#
# This single document always reports the pseudobulk merge and appends global
# single-cell merge diagnostics when that optional branch was run.

# %%
import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from IPython.display import Markdown, display
from matplotlib.ticker import MaxNLocator
from scipy import sparse
from upsetplot import UpSet, from_indicators

from pbmc_pipeline.reporting import gene_presence_indicators, pseudobulk_celltype_fractions

sns.set_theme(style="whitegrid")


def load_merge(prefix: str):
    path = Path(os.environ[f"QC_{prefix}_H5AD"])
    report = json.loads(Path(os.environ[f"QC_{prefix}_REPORT"]).read_text())
    return sc.read_h5ad(path), report, path


def show_summary(adata, report, path, heading):
    display(Markdown(f"## {heading}\n\nInput: `{path}`"))
    display(pd.Series({key: report.get(key) for key in ["kind", "status", "timestamp_utc", "n_observations", "n_genes"]}, name="value").to_frame())
    print(f"{adata.n_obs:,} observations × {adata.n_vars:,} genes")


def show_study_coverage(adata):
    study_sizes = adata.obs["study"].astype(str).value_counts().sort_values()
    _, axes = plt.subplots(1, 2, figsize=(13, 4))
    study_sizes.plot.barh(ax=axes[0], color="#4c72b0")
    axes[0].set_title("Observations by study")
    if "n_cells" in adata.obs:
        sns.boxplot(data=adata.obs, x="study", y="n_cells", ax=axes[1])
        axes[1].set_yscale("log")
        axes[1].set_title("Cells represented by each pseudobulk")
    else:
        adata.obs.groupby("study", observed=True)["sample"].nunique().sort_values().plot.barh(ax=axes[1], color="#55a868")
        axes[1].set_title("Distinct samples by study")
    plt.tight_layout()
    plt.show()
    display(pd.DataFrame({"observations": study_sizes, "samples": adata.obs.groupby("study", observed=True)["sample"].nunique(), "subjects": adata.obs.groupby("study", observed=True)["subject"].nunique()}))


def show_technical_covariates(adata):
    covariates = ["technology", "include_intronic", "frozen"]
    available = [column for column in covariates if column in adata.obs]
    if not available:
        return

    metadata = adata.obs[["study", "sample", *available]].copy()
    metadata["study"] = metadata["study"].astype(str)
    fig, axes = plt.subplots(1, len(available), figsize=(5 * len(available), 4), squeeze=False)
    for axis, column in zip(axes.flat, available):
        values = metadata[["study", column]].copy()
        values[column] = values[column].astype("string").fillna("not_provided").astype(str)
        study_counts = (
            values.drop_duplicates()
            .groupby(column, observed=True)["study"]
            .nunique()
            .sort_values(ascending=True)
        )
        axis.barh(study_counts.index, study_counts.values, color="#4c72b0")
        axis.set_title(column.replace("_", " "))
        axis.set_xlabel("Studies")
        axis.set_ylabel("Value")
        axis.set_xlim(0, max(1, int(study_counts.max())))
        axis.xaxis.set_major_locator(MaxNLocator(integer=True))

    display(Markdown("### Technical covariates across studies"))
    display(Markdown(
        "Bars count distinct studies represented by each value. A study with "
        "multiple values for a covariate is counted under each value; missing "
        "values are shown as `not_provided`."
    ))
    fig.tight_layout()
    plt.show()


def show_gene_join_accounting(report, heading):
    accounting = report.get("gene_join_accounting")
    if not accounting:
        display(Markdown(f"### {heading}\n\n_Count accounting is unavailable in this merge report._"))
        return

    join = accounting["gene_join"]
    join_description = "union" if join == "outer" else "intersection"
    total = int(accounting["input_counts"])
    discarded = int(accounting["discarded_counts"])
    percent = 100 * float(accounting["discarded_fraction"])
    display(Markdown(
        f"### {heading}\n\nThe {join_description} retains "
        f"{accounting['joined_gene_symbols']:,} gene symbols. Across all studies, "
        f"{discarded:,} of {total:,} input counts ({percent:.3f}%) are excluded by "
        "the gene join. The table reports each study's counts and the excluded gene "
        "symbol with the largest count contribution."
    ))
    rows = pd.DataFrame(accounting["studies"])
    if rows.empty:
        return
    rows["counts discarded (%)"] = (100 * rows["discarded_fraction"]).map(lambda value: f"{value:.3f}")
    rows["top excluded gene"] = rows.apply(
        lambda row: (
            f"{row['top_discarded_gene']} ({int(row['top_discarded_gene_counts']):,})"
            if row["top_discarded_gene"] is not None
            else "—"
        ),
        axis=1,
    )
    if report.get("kind") == "single_cell_merge":
        plot_rows = rows.sort_values("discarded_fraction", ascending=False)
        percentages = 100 * plot_rows["discarded_fraction"]
        fig, axis = plt.subplots(figsize=(9, max(3, 0.4 * len(plot_rows))))
        bars = axis.barh(plot_rows["study"], percentages, color="#c44e52")
        axis.invert_yaxis()
        axis.set_xlabel("Counts discarded (%)")
        axis.set_ylabel("Study")
        axis.set_title("Counts removed by the single-cell gene intersection")
        axis.set_xlim(0, max(0.5, float(percentages.max()) * 1.15))
        axis.bar_label(bars, labels=[f"{value:.2f}%" for value in percentages], padding=3)
        fig.tight_layout()
        plt.show()

    display(rows[[
        "study", "input_gene_symbols", "retained_gene_symbols", "discarded_gene_symbols",
        "input_counts", "discarded_counts", "counts discarded (%)", "top excluded gene",
    ]].rename(columns={
        "input_gene_symbols": "input genes",
        "retained_gene_symbols": "retained genes",
        "discarded_gene_symbols": "discarded genes",
        "input_counts": "input counts",
        "discarded_counts": "discarded counts",
    }))


# %% [markdown]
# ## Pseudobulk merge

# %%
pseudobulk, pseudobulk_report, pseudobulk_path = load_merge("PSEUDOBULK")
show_summary(pseudobulk, pseudobulk_report, pseudobulk_path, "Pseudobulk merge")
show_study_coverage(pseudobulk)
show_gene_join_accounting(pseudobulk_report, "Pseudobulk gene join")
show_technical_covariates(pseudobulk)

# %% [markdown]
# ### AIFI L2 composition by cell

# %%
composition = pseudobulk_celltype_fractions(pseudobulk.obs)
plt.figure(figsize=(max(10, 0.55 * composition.shape[1]), max(4, 0.55 * composition.shape[0])))
sns.heatmap(composition, cmap="viridis", vmin=0)
plt.title("AIFI L2 fraction of cells within each study")
plt.xlabel("AIFI L2")
plt.ylabel("Study")
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Gene-presence overlap

# %%
presence = gene_presence_indicators(pseudobulk.var)
if presence.empty or presence.shape[1] < 2:
    display(Markdown("_Gene-presence overlap requires an outer gene join across at least two studies._"))
else:
    # UpSetPlot 0.9 has a rendering error with ``show_counts`` under current
    # Matplotlib, so intersection labels are intentionally omitted.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module=r"upsetplot\.plotting")
        axes = UpSet(from_indicators(presence), subset_size="count", sort_by="cardinality").plot()
    figure = axes["intersections"].figure
    figure.suptitle("Gene presence across studies", y=1.02)
    figure.subplots_adjust(top=0.9)
    plt.show()

# %% [markdown]
# ## Global single-cell merge

# %%
single_cell_available = "QC_SINGLE_CELL_H5AD" in os.environ
if not single_cell_available:
    display(Markdown("_The global single-cell merge was not requested for this run._"))
else:
    single_cell, single_cell_report, single_cell_path = load_merge("SINGLE_CELL")
    show_summary(single_cell, single_cell_report, single_cell_path, "Global single-cell merge")
    show_study_coverage(single_cell)
    show_gene_join_accounting(single_cell_report, "Single-cell gene join")

    display(Markdown("### AIFI L2 label concordance"))
    label_columns = [
        ("Per-study L2 (downstream ground truth)", "aifi_l2_majority"),
        ("Experimental Harmony graph", "experimental_aifi_l2_majority"),
        ("Experimental unintegrated PCA graph", "experimental_aifi_l2_unintegrated_majority"),
    ]
    available_labels = [item for item in label_columns if item[1] in single_cell.obs]
    if len(available_labels) < 2:
        display(Markdown("_Fewer than two AIFI L2 label sets were available for comparison._"))
    else:
        comparisons = [
            (available_labels[left], available_labels[right])
            for left in range(len(available_labels))
            for right in range(left + 1, len(available_labels))
        ]
        fig, axes = plt.subplots(1, len(comparisons), figsize=(7 * len(comparisons), 6), squeeze=False)
        for axis, ((left_name, left_column), (right_name, right_column)) in zip(axes.flat, comparisons):
            concordance = pd.crosstab(
                single_cell.obs[left_column].astype(str),
                single_cell.obs[right_column].astype(str),
                normalize="index",
            )
            sns.heatmap(concordance, cmap="viridis", vmin=0, annot=concordance.shape[0] <= 15, fmt=".2f", ax=axis)
            axis.set_title(f"{left_name} vs {right_name}")
            axis.set_xlabel(right_name)
            axis.set_ylabel(left_name)
        fig.tight_layout()
        plt.show()
    display(Markdown("### Embedding and depth"))
    values = single_cell.X
    totals = np.asarray(values.sum(axis=1)).ravel()
    detected = np.diff(values.tocsr().indptr) if sparse.issparse(values) else np.count_nonzero(values, axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    sns.boxplot(data=pd.DataFrame({"study": single_cell.obs["study"].astype(str), "log10 library size": np.log10(totals + 1)}), x="study", y="log10 library size", ax=axes[0])
    sns.boxplot(data=pd.DataFrame({"study": single_cell.obs["study"].astype(str), "detected genes": detected}), x="study", y="detected genes", ax=axes[1])
    for axis in axes:
        axis.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.show()
    umap_colors = ["study", "aifi_l2_majority"]
    if "experimental_aifi_l2_majority" in single_cell.obs:
        umap_colors.append("experimental_aifi_l2_majority")
    if "experimental_aifi_l2_unintegrated_majority" in single_cell.obs:
        umap_colors.append("experimental_aifi_l2_unintegrated_majority")
    for color in umap_colors:
        sc.pl.umap(
            single_cell,
            color=color,
            size=max(2, 120000 / single_cell.n_obs),
            legend_loc="right margin",
            show=False,
        )
        figure = plt.gcf()
        figure.set_size_inches(12, 7)
        figure.tight_layout()
        plt.show()

    display(Markdown("### Single-cell QC UMAPs"))
    single_cell.obs["log_UMIs_per_cell"] = np.log1p(single_cell.obs["UMIs_per_cell"])
    for color, kwargs in [
        ("log_UMIs_per_cell", {}),
        ("percent_mito", {"vmin": 0, "vmax": 15}),
    ]:
        sc.pl.umap(
            single_cell,
            color=color,
            size=max(2, 120000 / single_cell.n_obs),
            show=False,
            **kwargs,
        )
        figure = plt.gcf()
        figure.set_size_inches(10, 7)
        figure.tight_layout()
        plt.show()

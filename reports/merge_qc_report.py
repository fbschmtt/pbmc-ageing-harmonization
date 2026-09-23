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
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
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


# %% [markdown]
# ## Pseudobulk merge

# %%
pseudobulk, pseudobulk_report, pseudobulk_path = load_merge("PSEUDOBULK")
show_summary(pseudobulk, pseudobulk_report, pseudobulk_path, "Pseudobulk merge")
show_study_coverage(pseudobulk)

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

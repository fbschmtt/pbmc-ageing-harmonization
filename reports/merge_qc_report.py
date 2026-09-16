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

from pbmc_pipeline.reporting import aifi_l2_concordance, gene_presence_indicators

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
# ### AIFI L2 composition

# %%
composition = pd.crosstab(pseudobulk.obs["study"].astype(str), pseudobulk.obs["aifi_l2_majority"].astype(str), normalize="index")
plt.figure(figsize=(max(10, 0.55 * composition.shape[1]), max(4, 0.55 * composition.shape[0])))
sns.heatmap(composition, cmap="viridis", vmin=0)
plt.title("AIFI L2 fraction within each study")
plt.xlabel("AIFI L2")
plt.ylabel("Study")
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Gene-presence overlap

# %%
presence = gene_presence_indicators(pseudobulk.var)
if presence.empty:
    display(Markdown("_No outer gene join was used, so this merge has no synthetic-zero gene-presence annotations._"))
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

    display(Markdown("### Merged versus individual-study AIFI L2 labels"))
    concordance = aifi_l2_concordance(single_cell.obs)
    if concordance.empty:
        display(Markdown("_Individual-study AIFI L2 labels were unavailable for comparison._"))
    else:
        plt.figure(figsize=(max(8, 0.6 * concordance.shape[1]), max(6, 0.45 * concordance.shape[0])))
        sns.heatmap(concordance, cmap="viridis", vmin=0, annot=concordance.shape[0] <= 15, fmt=".2f")
        plt.title("Merged AIFI L2 prediction by original study-level AIFI L2 label")
        plt.xlabel("Merged AIFI L2 prediction")
        plt.ylabel("Individual-study AIFI L2 label")
        plt.tight_layout()
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
    if "aifi_l2_study_majority" in single_cell.obs:
        umap_colors.insert(1, "aifi_l2_study_majority")
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

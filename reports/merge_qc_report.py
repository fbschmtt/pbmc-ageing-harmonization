# %%
import json
import os
from html import escape
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from IPython.display import HTML, Markdown, display
from matplotlib.ticker import MaxNLocator
from scipy import sparse

from pbmc_pipeline.report_theme import render_report_header
from pbmc_pipeline.reporting import pseudobulk_celltype_fractions

sns.set_theme(style="whitegrid")


def load_merge(prefix: str):
    path = Path(os.environ[f"QC_{prefix}_H5AD"])
    report = json.loads(Path(os.environ[f"QC_{prefix}_REPORT"]).read_text())
    return sc.read_h5ad(path), report, path


def show_summary(adata, report, path, heading):
    metadata = pd.Series(
        {
            key: report.get(key)
            for key in ["kind", "status", "timestamp_utc", "n_observations", "n_genes"]
        },
        name="value",
    ).to_frame()
    display(HTML(
        f'<details class="report-details"><summary>{escape(heading)} Metadata</summary>'
        f'<p><strong>Input:</strong> {escape(str(path))}</p>'
        + metadata.to_html(escape=True, border=0)
        + "</details>"
    ))


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

    display(Markdown("### Technical Covariates Across Studies"))
    display(Markdown(
        "Bars count distinct studies represented by each value. A study with "
        "multiple values for a covariate is counted under each value; missing "
        "values are shown as `not_provided`."
    ))
    fig.tight_layout()
    plt.show()


pseudobulk, pseudobulk_report, pseudobulk_path = load_merge("PSEUDOBULK")
single_cell_available = "QC_SINGLE_CELL_H5AD" in os.environ
if single_cell_available:
    single_cell, single_cell_report, single_cell_path = load_merge("SINGLE_CELL")

report_metrics = [
    ("Studies", f"{pseudobulk.obs['study'].nunique():,}"),
    ("Pseudobulk profiles", f"{pseudobulk.n_obs:,}"),
    ("Pseudobulk genes", f"{pseudobulk.n_vars:,}"),
    (
        "Single-cell cells",
        f"{single_cell.n_obs:,}" if single_cell_available else "Not included",
    ),
]
display(HTML(render_report_header(
    title="Cross-Study Merge QC",
    eyebrow="PBMC ageing · merge report",
    subtitle=(
        "Study coverage and technical metadata for the pseudobulk merge, "
        "with global single-cell diagnostics when available."
    ),
    metrics=report_metrics,
)))


# %% [markdown]
# ## Pseudobulk Merge

# %%
show_summary(pseudobulk, pseudobulk_report, pseudobulk_path, "Pseudobulk Merge")
display(Markdown(
    "The pseudobulk merge retains the **outer union** of study genes so per-study models "
    "can use observed counts. The cross-cell-type expression report records the corresponding "
    "shared-gene intersection, count-loss accounting, and availability overlap."
))
show_study_coverage(pseudobulk)
show_technical_covariates(pseudobulk)

# %% [markdown]
# ### AIFI L2 Composition by Cell

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
# ## Global Single-Cell Merge

# %%
if not single_cell_available:
    display(Markdown("_The global single-cell merge was not requested for this run._"))
else:
    show_summary(single_cell, single_cell_report, single_cell_path, "Global Single-Cell Merge")
    show_study_coverage(single_cell)

    display(Markdown("### Depth"))
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
    display(Markdown(
        "_The core merge intentionally stores raw counts, metadata, and QC only. "
        "Run `make run-integration-benchmark` for global UMAP and diagnostic "
        "CellTypist comparisons._"
    ))

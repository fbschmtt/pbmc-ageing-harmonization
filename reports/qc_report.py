# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display
from scipy import sparse

from pbmc_pipeline.report_theme import render_report_header

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message=r"Importing read_.* from `anndata` is deprecated",
        category=FutureWarning,
        module=r"anndata\.utils",
    )
    import scanpy as sc

sns.set_theme(style="whitegrid")
sc.settings.verbosity = 1
input_path = Path(os.environ["QC_INPUT_H5AD"])
run_report_path = Path(os.environ["QC_RUN_REPORT"])
study = os.environ["QC_STUDY"]
run_report = json.loads(run_report_path.read_text())
display(HTML(render_report_header(
    title=study,
    eyebrow="PBMC ageing · per-study QC",
    subtitle="Harmonized single-cell metadata, labels, and expression diagnostics.",
    metrics=[
        ("Status", str(run_report.get("status", "unknown"))),
        ("Cells", f"{run_report.get('n_cells', 0):,}"),
        ("Genes", f"{run_report.get('n_genes', 0):,}"),
        ("Excluded input cells", f"{run_report.get('excluded_input_cells', 0):,}"),
    ],
)))

# %% [markdown]
# ## Run provenance and dimensions

# %%
summary_keys = ["status", "timestamp_utc", "input_cells", "input_genes", "n_cells", "n_genes", "aggregated_duplicate_features", "configuration_sha256"]
display(pd.DataFrame({"value": [run_report.get(key) for key in summary_keys]}, index=summary_keys))
display(Markdown("**Warnings:** " + ("; ".join(run_report.get("warnings", [])) or "none")))
display(Markdown("**Recorded assumptions/provenance:**"))
display(pd.Series(run_report.get("provenance", {}), name="value").to_frame())
adata_backed = sc.read_h5ad(input_path, backed="r")
obs = adata_backed.obs.copy()
print(f"Stored object: {adata_backed.n_obs:,} cells × {adata_backed.n_vars:,} genes")

# %% [markdown]
# ## Metadata completeness

# %%
missing = obs.isna().mean()
not_provided = pd.Series({column: obs[column].astype(str).eq("not_provided").mean() for column in obs.columns})
completeness = pd.DataFrame({"missing": missing, "not_provided": not_provided}).sort_values(["missing", "not_provided"], ascending=False)
ax = completeness.plot.barh(stacked=True, figsize=(9, max(5, 0.28 * len(completeness))), color=["#c44e52", "#dd8452"])
ax.set_xlabel("Fraction of cells")
ax.set_title("Missing and explicitly unavailable metadata")
plt.tight_layout(); plt.show()
display(completeness.style.format("{:.1%}"))

# %% [markdown]
# ## Cohort, samples, and batches

# %%
subject_columns = [column for column in ["subject", "age", "sex", "bmi", "ethnicity", "cmv", "smoking_status"] if column in obs]
subjects = obs[subject_columns].drop_duplicates("subject") if "subject" in subject_columns else pd.DataFrame()
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
if not subjects.empty and subjects["age"].notna().any():
    sns.histplot(subjects, x="age", hue="sex" if "sex" in subjects else None, multiple="stack", bins=20, ax=axes[0])
    axes[0].set_title(f"Age distribution ({len(subjects):,} subjects)")
sample_sizes = obs["sample"].astype(str).value_counts().head(30).sort_values()
sample_sizes.plot.barh(ax=axes[1], color="#4c72b0")
axes[1].set_title("Cells per sample (30 largest)"); axes[1].set_xlabel("Cells")
plt.tight_layout(); plt.show()
print(f"Subjects: {obs['subject'].nunique():,}; samples: {obs['sample'].nunique():,}; batches: {obs['batch_single_cell'].nunique():,}")

# %% [markdown]
# ## Cell-type composition

# %%
label_columns = [column for column in ["study_celltype", "aifi_l1_majority", "aifi_l2_majority", "aifi_l3_majority"] if column in obs]
fig, axes = plt.subplots(len(label_columns), 1, figsize=(11, 4.5 * len(label_columns)))
axes = np.atleast_1d(axes)
for ax, column in zip(axes, label_columns):
    counts = obs[column].astype(str).value_counts().head(25).sort_values()
    counts.plot.barh(ax=ax, color="#55a868")
    ax.set_title(f"{column} (25 most frequent)", pad=16); ax.set_xlabel("Cells")
plt.tight_layout(h_pad=3.0); plt.show()

# %% [markdown]
# ## Expression and embedding checks
#
# Plots below use all cells in the harmonized study.

# %%
adata = adata_backed.to_memory()
adata_backed.file.close()
values = adata.X
total_counts = np.asarray(values.sum(axis=1)).ravel()
detected_genes = np.asarray((values > 0).sum(axis=1)).ravel() if sparse.issparse(values) else (values > 0).sum(axis=1)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
sns.histplot(np.log10(total_counts + 1), bins=40, ax=axes[0], color="#4c72b0")
axes[0].set_xlabel("log10(total counts + 1)")
sns.histplot(detected_genes, bins=40, ax=axes[1], color="#55a868")
axes[1].set_xlabel("Detected genes")
fig.suptitle(f"Count-depth checks ({adata.n_obs:,} cells)")
plt.tight_layout(); plt.show()

# %%
umap_columns = [column for column in ["study_celltype", "aifi_l1_majority", "aifi_l2_majority"] if column in adata.obs]
if "X_umap" in adata.obsm and umap_columns:
    sc.pl.umap(adata, color=umap_columns, ncols=2, size=max(2, 120000 / adata.n_obs), show=False)
    plt.show()
else:
    display(Markdown("_No stored UMAP coordinates were available._"))

# %% [markdown]
# ## Canonical PBMC marker expression

# %%
marker_panels = {
    "T cells": ["CD3D", "CD3E", "CD4", "CD8A", "CCR7"],
    "NK cells": ["NKG7", "GNLY", "KLRD1"],
    "B cells": ["MS4A1", "CD79A", "CD37"],
    "Monocytes/DC": ["LYZ", "S100A8", "FCGR3A", "FCER1A"],
    "Platelets": ["PPBP", "PF4"],
}
present = {panel: [gene for gene in genes if gene in adata.var_names] for panel, genes in marker_panels.items()}
present = {panel: genes for panel, genes in present.items() if genes}
groupby = "aifi_l2_majority" if "aifi_l2_majority" in adata.obs else "study_celltype"
if present:
    expression = adata.copy()
    sc.pp.normalize_total(expression, target_sum=1e4)
    sc.pp.log1p(expression)
    sc.pl.dotplot(expression, present, groupby=groupby, standard_scale="var", show=False)
    plt.show()
else:
    display(Markdown("_None of the configured canonical markers were present._"))

# %% [markdown]
# ## Study labels versus AIFI predictions

# %%
if {"study_celltype", "aifi_l2_majority"} <= set(obs.columns):
    top_study = obs["study_celltype"].astype(str).value_counts().head(20).index
    top_aifi = obs["aifi_l2_majority"].astype(str).value_counts().head(20).index
    selected = obs["study_celltype"].astype(str).isin(top_study) & obs["aifi_l2_majority"].astype(str).isin(top_aifi)
    agreement = pd.crosstab(obs.loc[selected, "study_celltype"], obs.loc[selected, "aifi_l2_majority"], normalize="index")
    plt.figure(figsize=(max(8, 0.5 * agreement.shape[1]), max(6, 0.4 * agreement.shape[0])))
    sns.heatmap(agreement, cmap="viridis", vmin=0, vmax=1)
    plt.title("AIFI L2 composition within study-provided labels")
    plt.xlabel("AIFI L2"); plt.ylabel("Study label"); plt.tight_layout(); plt.show()

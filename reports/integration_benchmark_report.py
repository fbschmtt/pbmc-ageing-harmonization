# %% [markdown]
# # Global integration benchmark
#
# This report compares diagnostic CellTypist majority-voting labels produced
# from the configured global graph variants. The retained per-study AIFI-L2
# call is the downstream label; benchmark labels are diagnostic only.

# %%
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from IPython.display import Markdown, display

from pbmc_pipeline.reporting import aifi_l2_concordance

sns.set_theme(style="whitegrid")
artifact_path = Path(os.environ["INTEGRATION_BENCHMARK_H5AD"])
run_report_path = Path(os.environ["INTEGRATION_BENCHMARK_REPORT"])
adata = sc.read_h5ad(artifact_path)
run_report = json.loads(run_report_path.read_text())

# %% [markdown]
# ## Benchmark scope

# %%
display(Markdown(
    f"**Input:** `{run_report['input']}`  \\n+**Cells:** {adata.n_obs:,}  \\n+**Input genes:** {run_report['n_genes_input']:,}  \\n+**Harmony graph reused for UMAP:** {run_report['harmony_graph_reused_for_umap']}"
))
display(pd.DataFrame(run_report["methods"])[["name", "label_column", "basis", "integration"]].fillna("—"))
display(Markdown(
    "_The benchmark artifact stores labels and UMAP coordinates only. The raw "
    "matrix remains in the core merged H5AD._"
))

# %% [markdown]
# ## AIFI L2 label concordance
#
# Each row is a per-study AIFI-L2 label and each column is a diagnostic
# benchmark prediction. Values are row-normalized cell fractions.

# %%
reference_column = "aifi_l2_study_majority"
methods = run_report["methods"]
figure, axes = plt.subplots(1, len(methods), figsize=(7 * len(methods), 6), squeeze=False)
for axis, method in zip(axes.flat, methods):
    prediction_column = method["label_column"]
    concordance = aifi_l2_concordance(
        adata.obs[[reference_column, prediction_column]].rename(
            columns={reference_column: "aifi_l2_majority"}
        ),
        prediction_column,
    )
    sns.heatmap(
        concordance, cmap="viridis", vmin=0, annot=concordance.shape[0] <= 15,
        fmt=".2f", ax=axis,
    )
    axis.set_title(method["name"])
    axis.set_xlabel("Benchmark AIFI L2")
    axis.set_ylabel("Per-study AIFI L2")
figure.tight_layout()
plt.show()

# %% [markdown]
# ## Global UMAP
#
# The UMAP uses the Harmony graph. The unintegrated-PCA labels are projected
# onto the same coordinates for comparison; they do not define this embedding.

# %%
for color, title in [
    ("study", "Study"),
    (reference_column, "Per-study AIFI L2"),
    *[(method["label_column"], f"Benchmark AIFI L2: {method['name']}") for method in methods],
]:
    sc.pl.umap(
        adata, color=color, title=title, legend_loc="right margin",
        size=max(2, 120000 / adata.n_obs), show=False,
    )
    figure = plt.gcf()
    figure.set_size_inches(12, 7)
    figure.tight_layout()
    plt.show()

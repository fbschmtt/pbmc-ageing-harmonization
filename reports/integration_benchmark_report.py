# %%
import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display

from pbmc_pipeline.cell_type_colors import set_aifi_l2_scanpy_palette
from pbmc_pipeline.report_theme import render_report_header
from pbmc_pipeline.reporting import aifi_l2_concordance

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message=r"Importing read_.* from `anndata` is deprecated",
        category=FutureWarning,
        module=r"anndata\.utils",
    )
    import scanpy as sc

sns.set_theme(style="whitegrid")
artifact_path = Path(os.environ["INTEGRATION_BENCHMARK_H5AD"])
run_report_path = Path(os.environ["INTEGRATION_BENCHMARK_REPORT"])
adata = sc.read_h5ad(artifact_path)
run_report = json.loads(run_report_path.read_text())
sampling = run_report.get("sampling", {
    "n_cells": adata.n_obs,
    "source_n_cells": adata.n_obs,
    "fraction": 1.0,
    "strategy": "none",
    "seed": "not recorded",
})
method_labels = {
    "unintegrated_pca": "unintegrated PCA",
    "harmony": "Harmony",
    "scanorama": "Scanorama",
    "bbknn": "BBKNN",
}
active_method_labels = [
    method_labels.get(method["name"], method["name"].replace("_", " "))
    for method in run_report["methods"]
]

# %%
display(HTML(render_report_header(
    title="Global integration benchmark",
    eyebrow="PBMC ageing · integration and type calling",
    subtitle=(
        f"Compare {', '.join(active_method_labels)} embeddings, "
        "then assess CellTypist calls made from each method's own neighbor graph."
    ),
    metrics=[
        ("Cells", f"{adata.n_obs:,}"),
        ("Studies", f"{adata.obs['study'].nunique():,}"),
        ("Input genes", f"{run_report['n_genes_input']:,}"),
        ("Methods", str(len(run_report["methods"]))),
    ],
)))
display(Markdown(
    f"This benchmark used **{sampling['n_cells']:,} of "
    f"{sampling['source_n_cells']:,} cells** "
    f"({sampling['fraction']:.1%}; "
    f"{sampling['strategy'].replace('_', ' ')}, "
    f"seed {sampling['seed']}). All methods were evaluated on "
    "the same sample; benchmark embeddings and metrics do not cover cells "
    "outside it."
))

# %% [markdown]
# ## 1. Embedding integration
#
# Compare the study-mixed geometry before and after batch integration. Each
# embedding is paired with the graph used to produce it.

# %%
reference_column = "aifi_l2_study_majority"
methods = run_report["methods"]
method_order = {
    "unintegrated_pca": 0,
    "harmony": 2,
    "scanorama": 3,
    "bbknn": 4,
}
embedding_methods = sorted(
    methods,
    key=lambda method: (method_order.get(method["name"], 99), method["name"]),
)
configured_method_names = [method["name"] for method in embedding_methods]
scanorama_note = (
    " Scanorama is temporarily disabled in the active configuration because "
    "its correction step can allocate large dense temporary arrays."
    if "scanorama" not in configured_method_names else ""
)
display(Markdown(
    "Each method's graph is used for both its UMAP and CellTypist majority "
    "voting. Sparse graph overlap is calculated during aggregation; the final "
    "artifact retains labels and UMAPs, but omits expression, PCA, and graphs. "
    "**scVI was not included** because the "
    "expected input is around six million cells; model training at that scale "
    "would be highly resource intensive and may exceed available GPU memory. "
    "Configured methods in this run: "
    + ", ".join(name.replace("_", " ") for name in configured_method_names)
    + "."
    + scanorama_note
))
failed_methods = run_report.get("failed_methods", [])
if failed_methods:
    display(Markdown(
        "**Partial benchmark:** these configured methods did not produce an artifact: "
        + "; ".join(
            f"{method['name']} ({method.get('status', 'failed')}): "
            f"{method.get('error', 'no details')}"
            for method in failed_methods
        )
        + "."
    ))
display(pd.DataFrame([
    {
        "method": method["name"],
        "embedding": method["embedding_key"],
        "graph representation": method.get("neighbors_use_rep", method.get("basis")),
        "integration": method.get("integration"),
    }
    for method in embedding_methods
]))

for method in embedding_methods:
    heading = {
        "unintegrated_pca": "Unintegrated",
        "harmony": "Harmony integrated",
        "scanorama": "Scanorama integrated",
        "bbknn": "BBKNN integrated neighbors",
    }.get(method["name"], method["name"].replace("_", " ").title())
    display(Markdown(f"### {heading}"))
    key = method["embedding_key"]
    set_aifi_l2_scanpy_palette(adata, reference_column)
    for color, title in [("study", "Study"), (reference_column, "Per-study AIFI L2")]:
        sc.pl.embedding(
            adata, basis=key, color=color, title=f"{method['name']}: {title}",
            legend_loc="right margin", size=max(2, 120000 / adata.n_obs), show=False,
        )
        figure = plt.gcf()
        figure.set_size_inches(12, 7)
        figure.tight_layout()
        plt.show()

# %% [markdown]
# ## 2. Type calling
#
# CellTypist majority-voting predictions are generated on each method's own
# graph. Sparse edge overlap is calculated during aggregation and recorded
# below; the published H5AD omits the graph matrices.
# Concordance rows are per-study AIFI-L2 labels and columns are benchmark calls;
# values are row-normalized cell fractions.

# %%
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

# %%
graph_comparisons = run_report.get("graph_comparisons", [])
if graph_comparisons:
    display(Markdown("### Neighbor graph overlap against unintegrated PCA"))
    display(pd.DataFrame(graph_comparisons))

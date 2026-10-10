# %%
import json
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display

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

# %%
display(HTML(render_report_header(
    title="Global integration benchmark",
    eyebrow="PBMC ageing · integration and type calling",
    subtitle=(
        "Compare unintegrated PCA and scVI, Harmony, Scanorama, and BBKNN embeddings, "
        "then assess CellTypist calls made from each method's own neighbor graph."
    ),
    metrics=[
        ("Cells", f"{adata.n_obs:,}"),
        ("Studies", f"{adata.obs['study'].nunique():,}"),
        ("Input genes", f"{run_report['n_genes_input']:,}"),
        ("Methods", str(len(run_report["methods"]))),
    ],
)))

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
    "scvi_unintegrated": 1,
    "harmony": 2,
    "scanorama": 3,
    "bbknn": 4,
    "scvi": 5,
}
embedding_methods = sorted(
    methods,
    key=lambda method: (method_order.get(method["name"], 99), method["name"]),
)
display(Markdown(
    "Each method's graph is used for both its UMAP and CellTypist majority "
    "voting. The artifact retains both embeddings and neighbor graphs, but "
    "omits the merged expression matrix. scVI produces a latent representation; "
    "the benchmark builds neighbors and UMAP from that representation. The "
    "unintegrated scVI model omits the study covariate, while both scVI runs "
    "use the same study-aware HVG selection."
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
        "scvi_unintegrated": "Unintegrated scVI latent",
        "scvi": "scVI integrated latent",
    }.get(method["name"], method["name"].replace("_", " ").title())
    display(Markdown(f"### {heading}"))
    key = method["embedding_key"]
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
# graph. The stored connectivities also let us compare graph membership.
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
def edge_set(connectivities):
    graph = connectivities.tocsr(copy=True)
    graph.setdiag(0)
    graph.eliminate_zeros()
    rows, cols = graph.nonzero()
    return set(zip(rows[rows < cols].tolist(), cols[rows < cols].tolist()))

baseline = next((method for method in methods if method["name"] == "unintegrated_pca"), None)
if baseline is not None:
    baseline_edges = edge_set(adata.obsp[baseline["graph_keys"]["connectivities"]])
    rows = []
    for method in methods:
        if method is baseline:
            continue
        edges = edge_set(adata.obsp[method["graph_keys"]["connectivities"]])
        union = baseline_edges | edges
        rows.append({
            "graphs": f"unintegrated vs {method['name']}",
            "shared_undirected_edges": len(baseline_edges & edges),
            "union_undirected_edges": len(union),
            "edge_jaccard": len(baseline_edges & edges) / len(union) if union else 1.0,
        })
    display(pd.DataFrame(rows))

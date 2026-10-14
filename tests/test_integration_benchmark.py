import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.integration_benchmark import _thin_benchmark_artifact


def test_thin_benchmark_artifact_keeps_labels_embeddings_and_graphs_without_counts():
    source = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [0, 3]]),
        obs=pd.DataFrame(
            {
                "aifi_l2_majority": ["T", "B"],
                "benchmark_harmony_aifi_l2_majority": ["T", "B"],
            },
            index=["cell_1", "cell_2"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    source.obsm["X_umap_harmony"] = np.array([[1.0, 2.0], [3.0, 4.0]])
    source.obsm["X_pca_harmony"] = np.array([[1.0], [2.0]])
    source.obsp["harmony_connectivities"] = sparse.csr_matrix([[0, 1], [1, 0]])
    source.obsp["harmony_distances"] = sparse.csr_matrix([[0, 2], [2, 0]])

    result = _thin_benchmark_artifact(source, methods=[{
        "name": "harmony",
        "embedding_key": "X_umap_harmony",
        "representation_key": "X_pca_harmony",
        "graph_keys": {
            "connectivities": "harmony_connectivities",
            "distances": "harmony_distances",
        },
    }])

    assert result.shape == (2, 0)
    assert result.obs.columns.tolist() == source.obs.columns.tolist()
    assert np.array_equal(result.obsm["X_umap_harmony"], source.obsm["X_umap_harmony"])
    assert np.array_equal(result.obsm["X_pca_harmony"], source.obsm["X_pca_harmony"])
    assert set(result.obsp) == {"harmony_connectivities", "harmony_distances"}
    assert (result.obsp["harmony_connectivities"] != source.obsp["harmony_connectivities"]).nnz == 0
    assert result.uns["integration_benchmark"]["source_matrix_retained"] is False

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.integration_benchmark import _thin_benchmark_artifact


def test_thin_benchmark_artifact_keeps_labels_and_umap_without_counts():
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
    source.obsm["X_umap"] = np.array([[1.0, 2.0], [3.0, 4.0]])

    result = _thin_benchmark_artifact(source, methods=[{"name": "harmony"}])

    assert result.shape == (2, 0)
    assert result.obs.columns.tolist() == source.obs.columns.tolist()
    assert np.array_equal(result.obsm["X_umap"], source.obsm["X_umap"])
    assert result.uns["integration_benchmark"]["source_matrix_retained"] is False

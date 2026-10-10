import sys
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.integration_benchmark import (
    _compute_scvi_representation,
    _thin_benchmark_artifact,
)


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


def test_scvi_without_batch_covariate_still_uses_study_aware_hvgs(monkeypatch):
    source = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [0, 3]], dtype=np.int32),
        obs=pd.DataFrame({"study": ["one", "two"]}, index=["cell_1", "cell_2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    hvg_arguments = {}
    setup_arguments = {}

    def select_hvgs(adata, **kwargs):
        hvg_arguments.update(kwargs)
        adata.var["highly_variable"] = [True, True]

    class FakeSCVI:
        @staticmethod
        def setup_anndata(adata, **kwargs):
            setup_arguments.update(kwargs)

        def __init__(self, adata, n_latent):
            self.n_latent = n_latent

        def train(self, **kwargs):
            pass

        def get_latent_representation(self):
            return np.ones((2, self.n_latent))

    fake_scanpy = SimpleNamespace(
        experimental=SimpleNamespace(pp=SimpleNamespace(highly_variable_genes=select_hvgs))
    )
    fake_scvi = SimpleNamespace(
        settings=SimpleNamespace(seed=None),
        model=SimpleNamespace(SCVI=FakeSCVI),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    latent, details = _compute_scvi_representation(
        source,
        {
            "name": "scvi_unintegrated",
            "hvg_batch_key": "study",
            "n_top_genes": 2,
            "hvg_chunksize": 10,
            "n_latent": 3,
            "max_epochs": 2,
            "batch_size": 2,
        },
        random_seed=13,
        exclude_vdj_regex=r"TR[ABDG]",
    )

    assert latent.shape == (2, 3)
    assert hvg_arguments["batch_key"] == "study"
    assert setup_arguments == {}
    assert details["batch_key"] is None
    assert details["hvg_batch_key"] == "study"
    assert details["integration"] == "scVI latent without batch correction"

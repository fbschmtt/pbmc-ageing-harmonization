import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.cell_type import _marker_table, analyse_cell_type, split_cell_types

PIPELINE = {
    "processing": {"output_compression": "gzip"},
    "cell_type_analysis": {
        "split_by": "aifi_l2_majority",
        "l2_parent_l1": {"Memory B cell": "B cell", "Naive CD4 T cell": "T cell"},
    },
}


def test_split_cell_types_writes_raw_count_partitions_with_global_umap(tmp_path):
    source = ad.AnnData(
        X=sparse.csr_matrix([[1, 0], [2, 3], [4, 5]]),
        obs=pd.DataFrame(
            {
                "aifi_l2_majority": ["Memory B cell", "Memory B cell", "Naive CD4 T cell"],
                "study": ["one"] * 3,
                "sample": ["one", "one", "two"],
            },
            index=["a", "b", "c"],
        ),
        var=pd.DataFrame(index=["A", "B"]),
    )
    source.obsm["X_umap"] = np.ones((3, 2))
    source_path = tmp_path / "merged.h5ad"
    source.write_h5ad(source_path)

    document = split_cell_types(source_path, tmp_path / "cells", tmp_path / "cell_types.json", PIPELINE)

    assert [(item["cell_type"], item["slug"], item["n_cells"]) for item in document["cell_types"]] == [
        ("Memory B cell", "memory-b-cell", 2), ("Naive CD4 T cell", "naive-cd4-t-cell", 1)
    ]
    split = ad.read_h5ad(tmp_path / "cells" / "memory-b-cell.h5ad")
    assert split.obs_names.tolist() == ["a", "b"]
    assert set(split.obsm) == {"X_umap"}
    assert split.X.toarray().tolist() == [[1, 0], [2, 3]]
    assert split.obs["n_cells_in_sample"].tolist() == [2, 2]
    assert split.obs["n_cells_in_sample_l1_parent"].tolist() == [2, 2]
    assert json.loads((tmp_path / "cell_types.json").read_text())["split_by"] == "aifi_l2_majority"


def test_analysis_marks_small_cell_type_without_attempting_an_embedding(tmp_path):
    primitive = ad.AnnData(
        X=sparse.csr_matrix([[1, 0], [2, 3], [4, 5]]),
        obs=pd.DataFrame(
            {
                "study": ["one"] * 3,
                "sample": ["one"] * 3,
                "age": [20, 20, 20],
                "n_cells_in_sample": [3] * 3,
                "aifi_l2_majority": ["Naive CD4 T cell"] * 3,
                "aifi_l1_parent_for_l2": ["T cell"] * 3,
            },
            index=["a", "b", "c"],
        ),
        var=pd.DataFrame(index=["A", "B"]),
    )
    input_path = tmp_path / "small.h5ad"
    primitive.write_h5ad(input_path)
    pipeline = {
        "processing": {"output_compression": "gzip"},
        "cell_type_analysis": {"min_cells": 4},
    }

    report = analyse_cell_type(input_path, tmp_path / "analysis", pipeline)

    assert report["status"] == "insufficient_cells"
    assert ad.read_h5ad(tmp_path / "analysis" / "analysis.h5ad").n_obs == 3


def test_marker_table_adds_group_when_scanpy_returns_a_single_cluster_without_one():
    markers = pd.DataFrame({"names": ["A", "B"], "scores": [2.0, 1.0]})

    result = _marker_table(markers, pd.Series(["0", "0"]))

    assert result["group"].tolist() == ["0", "0"]


def test_analysis_uses_harmony_for_local_neighbors_but_keeps_native_pcs(tmp_path, monkeypatch):
    import scanpy as sc
    import scanpy.external as sce

    rng = np.random.default_rng(42)
    counts = rng.poisson(2, size=(30, 12)) + 1
    counts[:15, :6] += 8
    counts[15:, 6:] += 8
    primitive = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=pd.DataFrame(
            {
                "study": ["one"] * 15 + ["two"] * 15,
                "sample": ["a"] * 15 + ["b"] * 15,
                "age": [30.0] * 15 + [60.0] * 15,
                "n_cells_in_sample": [15] * 30,
                "aifi_l2_majority": ["Naive CD4 T cell"] * 30,
                "aifi_l1_parent_for_l2": ["T cell"] * 30,
            },
            index=[f"cell-{index}" for index in range(30)],
        ),
        var=pd.DataFrame(index=[f"gene-{index}" for index in range(12)]),
    )
    input_path = tmp_path / "type.h5ad"
    primitive.write_h5ad(input_path)
    observed = {}

    def fake_harmony(adata, key, *, basis, adjusted_basis, **kwargs):
        observed.update({"key": key, "basis": basis, "adjusted_basis": adjusted_basis})
        adata.obsm[adjusted_basis] = adata.obsm[basis].copy()

    def all_genes_are_hvg(adata, **kwargs):
        adata.var["highly_variable"] = True

    def two_clusters(adata, *, key_added, **kwargs):
        adata.obs[key_added] = ["0"] * 15 + ["1"] * 15

    monkeypatch.setattr(sce.pp, "harmony_integrate", fake_harmony)
    monkeypatch.setattr(sc.pp, "highly_variable_genes", all_genes_are_hvg)
    monkeypatch.setattr(sc.tl, "leiden", two_clusters)
    pipeline = {
        "processing": {
            "output_compression": "gzip", "target_sum": 10_000,
            "hvg": {"min_mean": 0.0125, "max_mean": 3.0, "min_disp": 0.5},
            "exclude_vdj_regex": "^$", "neighbors": {"n_neighbors": 10, "n_pcs": 5},
            "random_seed": 42,
        },
        "cell_type_analysis": {
            "min_cells": 5, "clustering": {"resolution": 0.5},
            "integration": {"method": "harmony", "batch_key": "study", "adjusted_basis": "X_pca_harmony"},
        },
    }

    report = analyse_cell_type(input_path, tmp_path / "analysis", pipeline)
    analysed = ad.read_h5ad(tmp_path / "analysis" / "analysis.h5ad")

    assert report["status"] == "complete"
    assert observed == {"key": "study", "basis": "X_pca", "adjusted_basis": "X_pca_harmony"}
    assert {"X_pca", "X_pca_harmony", "X_umap"}.issubset(analysed.obsm)
    assert analysed.uns["neighbors"]["params"]["use_rep"] == "X_pca_harmony"
    markers = pd.read_csv(tmp_path / "analysis" / "markers.tsv", sep="\t")
    assert markers["logfoldchanges"].notna().all()

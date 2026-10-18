import json
import types

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pbmc_pipeline.harmonize import _compute_embedding, _predict_celltypist
from pbmc_pipeline.merge import (
    _match_csr_index_dtypes,
    _shared_gene_names,
    _stringify_mixed_object_columns,
    merge_pseudobulks,
    pseudobulk_study,
)

PIPELINE = {
    "processing": {"output_compression": "gzip"},
    "merge": {"pseudobulk": {"groupby": ["sample", "aifi_l2_majority"]}},
}
OBS_SCHEMA = {
    "pseudobulk_metadata": ["study", "study_doi", "subject", "sex", "age"],
}


def _study(path, study, genes, rows):
    obs = pd.DataFrame(rows, index=[f"{study}_{i}" for i in range(len(rows))])
    obs["study"] = study
    obs["study_doi"] = f"10.0000/{study}"
    ad.AnnData(X=sparse.csr_matrix(np.arange(1, len(rows) * len(genes) + 1).reshape(len(rows), len(genes))), obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(path)


def test_mixed_object_metadata_is_stringified_without_losing_missingness():
    obs = pd.DataFrame({
        "age": pd.Series([35.0, "Heterogeneous During Bulk", np.nan], dtype=object),
    })

    repairs = _stringify_mixed_object_columns(obs)

    assert obs["age"].iloc[:2].tolist() == ["35.0", "Heterogeneous During Bulk"]
    assert pd.isna(obs["age"].iloc[2])
    assert repairs == [
        {"column": "age", "source_types": ["number", "string"], "converted_to": "string"},
    ]


def test_pseudobulk_is_sample_and_l2_specific_then_outer_merged(tmp_path, caplog):
    first = tmp_path / "one.h5ad"
    second = tmp_path / "two.h5ad"
    rows = [
        {
            "sample": "s1", "subject": "d1", "sex": "female", "aifi_l2_majority": "T",
            "age": 35.0,
            "n_genes_by_counts": 2, "percent_mito": 1.0,
            "aifi_l1_majority": "lymphocyte", "aifi_l3_majority": "naive T",
            "study_celltype": "T cell", "unmapped_source_obs": "debug only",
        },
        {
            "sample": "s1", "subject": "d1", "sex": "male", "aifi_l2_majority": "T",
            "age": 45.0,
            "n_genes_by_counts": 3, "percent_mito": 2.0,
            "aifi_l1_majority": "lymphocyte", "aifi_l3_majority": "memory T",
            "study_celltype": "T cell", "unmapped_source_obs": "debug only",
        },
        {
            "sample": "s1", "subject": "d1", "sex": "female", "aifi_l2_majority": "B",
            "age": 35.0,
            "n_genes_by_counts": 4, "percent_mito": 3.0,
            "aifi_l1_majority": "lymphocyte", "aifi_l3_majority": "naive B",
            "study_celltype": "B cell", "unmapped_source_obs": "debug only",
        },
    ]
    _study(first, "one", ["A", "B"], rows)
    _study(second, "two", ["B", "C"], rows[:1])
    one_pb, one_report = tmp_path / "one.pb.h5ad", tmp_path / "one.pb.json"
    two_pb, two_report = tmp_path / "two.pb.h5ad", tmp_path / "two.pb.json"
    report = pseudobulk_study(first, one_pb, one_report, PIPELINE, OBS_SCHEMA)
    pseudobulk_study(second, two_pb, two_report, PIPELINE, OBS_SCHEMA)
    assert report["n_pseudobulks"] == 2
    result = ad.read_h5ad(one_pb)
    assert result.obs.loc["one::s1::T", "n_cells"] == 2
    assert result.obs.loc["one::s1::T", "total_counts"] == 10
    assert result.obs.loc["one::s1::T", "subject"] == "d1"
    assert result.obs.loc["one::s1::T", "sex"] == "Heterogeneous During Bulk"
    assert "n_genes_by_counts" not in result.obs
    assert "percent_mito" not in result.obs
    assert "aifi_l1_majority" not in result.obs
    assert "aifi_l3_majority" not in result.obs
    assert "study_celltype" not in result.obs
    assert "unmapped_source_obs" not in result.obs
    assert isinstance(result.obs["study_doi"].dtype, pd.CategoricalDtype)
    assert result.obs["study_doi"].cat.categories.tolist() == ["10.0000/one"]
    assert result.uns["pipeline_provenance"]["artifact_kind"] == "pseudobulk"
    merged_path, merged_report = tmp_path / "merged.h5ad", tmp_path / "merged.json"
    merge_pseudobulks([one_pb, two_pb], merged_path, merged_report, PIPELINE)
    merged = ad.read_h5ad(merged_path)
    assert set(merged.var_names) == {"A", "B", "C"}
    assert merged.n_obs == 3
    assert bool(merged.var.loc["A", "has_synthetic_zeros"])
    assert merged.var.loc["A", "synthetic_zero_filled_studies"] == "two"
    assert not bool(merged.var.loc["B", "has_synthetic_zeros"])
    assert merged.uns["pipeline_provenance"]["gene_join"] == "outer"
    assert merged.obs.loc["one::s1::T", "age"] == "Heterogeneous During Bulk"
    assert merged.obs.loc["two::s1::T", "age"] == "35.0"
    expected_normalizations = [
        {"column": "age", "source_types": ["number", "string"], "converted_to": "string"},
    ]
    assert "mixed value types in obs column age" in caplog.text
    assert isinstance(merged.obs["study_doi"].dtype, pd.CategoricalDtype)
    assert merged.obs["study_doi"].cat.categories.tolist() == ["10.0000/one", "10.0000/two"]
    report = json.loads(merged_report.read_text())
    assert report["obs_type_normalizations"] == expected_normalizations
    assert report["observations_by_study"] == {"one": 2, "two": 1}
    assert report["aifi_l2_cells_by_study"] == [
        {"study": "one", "aifi_l2_majority": "B", "n_cells": 1},
        {"study": "one", "aifi_l2_majority": "T", "n_cells": 2},
        {"study": "two", "aifi_l2_majority": "T", "n_cells": 1},
    ]


def test_single_cell_embedding_uses_shared_genes_then_harmony_neighbors(monkeypatch):
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.random.default_rng(42).poisson(2, size=(80, 20))),
        obs=pd.DataFrame(
            {"study": ["one"] * 40 + ["two"] * 40},
            index=[f"cell_{index}" for index in range(80)],
        ),
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(20)]),
    )
    shared_genes = pd.Index([f"gene_{index}" for index in range(12)])
    seen: dict[str, object] = {}

    def fake_harmony(embedding, key, *, basis, adjusted_basis, random_state):
        seen.update({"genes": embedding.var_names.copy(), "key": key, "basis": basis})
        embedding.obsm[adjusted_basis] = embedding.obsm[basis].copy()
        return {"method": "harmony", "implementation": "harmonypy"}

    monkeypatch.setattr("pbmc_pipeline.harmonize._harmony_integrate", fake_harmony)
    details = _compute_embedding(
        adata,
        {
            "processing": {
                "hvg": {"min_mean": 0.0125, "min_disp": 0.5},
                "exclude_vdj_regex": "^VDJ",
                "neighbors": {"n_neighbors": 10, "n_pcs": 5},
                "random_seed": 42,
            }
        },
        "merged_single_cell",
        batch_key="study",
        embedding_genes=shared_genes,
    )

    assert set(seen["genes"]).issubset(set(shared_genes))
    assert seen["key"] == "study"
    assert seen["basis"] == "X_pca"
    assert details["input_genes"] == len(shared_genes)
    assert details["neighbors_use_rep"] == "X_pca_harmony"
    assert details["neighbors_n_pcs"] == min(5, adata.obsm["X_pca_harmony"].shape[1])
    assert "X_pca_harmony" in adata.obsm
    assert adata.uns["neighbors"]["params"]["use_rep"] == "X_pca_harmony"
    assert adata.uns["neighbors"]["params"]["n_pcs"] == details["neighbors_n_pcs"]


def test_shared_gene_names_preserves_first_study_order():
    first = ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["C", "A", "B"]))
    second = ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["B", "C"]))
    assert _shared_gene_names([first, second]).tolist() == ["C", "B"]


def test_matching_csr_index_dtypes_unblocks_scanpy_normalization():
    import scanpy as sc

    matrix = sparse.csr_matrix([[3, 1], [2, 4]], dtype=np.int64)
    matrix.indices = matrix.indices.astype(np.int32)
    matrix.indptr = matrix.indptr.astype(np.int64)
    assert matrix.indices.dtype != matrix.indptr.dtype
    adata = ad.AnnData(X=_match_csr_index_dtypes(matrix))

    assert adata.X.indices.dtype == adata.X.indptr.dtype == np.dtype(np.int32)
    sc.pp.normalize_total(adata, target_sum=10_000)
    assert np.allclose(np.asarray(adata.X.sum(axis=1)).ravel(), 10_000)


def test_matching_csr_index_dtypes_leaves_compatible_matrix_in_place():
    matrix = sparse.csr_matrix([[1, 0], [0, 1]], dtype=np.float32)

    assert _match_csr_index_dtypes(matrix) is matrix
    assert matrix.indices.dtype == matrix.indptr.dtype


def test_experimental_prediction_prefix_does_not_replace_per_study_l2(monkeypatch, tmp_path):
    adata = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"aifi_l2_majority": ["per-study T", "per-study B"]}, index=["a", "b"]),
    )
    adata.uns["neighbors"] = {"params": {"use_rep": "X_pca_harmony"}}

    fake_celltypist = types.SimpleNamespace(
        models=types.SimpleNamespace(Model=types.SimpleNamespace(load=lambda path: path)),
        annotate=lambda *args, **kwargs: types.SimpleNamespace(
            predicted_labels=pd.DataFrame({"majority_voting": ["experimental T", "experimental B"]}, index=adata.obs_names)
        ),
    )
    monkeypatch.setitem(__import__("sys").modules, "celltypist", fake_celltypist)

    _predict_celltypist(
        adata, {"levels": ["l2"]}, {"models": {"aifi_l2": "l2.pkl"}}, tmp_path,
        "merged", label_prefix="experimental_",
    )

    assert adata.obs["aifi_l2_majority"].tolist() == ["per-study T", "per-study B"]
    assert adata.obs["experimental_aifi_l2_majority"].tolist() == ["experimental T", "experimental B"]

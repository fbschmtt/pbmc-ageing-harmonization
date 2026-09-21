import anndata as ad
import numpy as np
import pandas as pd

from pbmc_pipeline.harmonize import (
    _normalize_missing_metadata,
    _normalize_nullable_strings_for_h5ad,
    _repair_features,
)
from pbmc_pipeline.metadata import build_homogeneous_obs


def test_missing_optional_source_uses_fallback():
    source = pd.DataFrame({"sampleName": ["S1", "S2"]}, index=["c1", "c2"])
    rules = {
        "sample": {"source": "sampleName"},
        "sex": {"source": "Gender", "fallback_constant": "not_provided", "lower": True},
        "bmi": {"source": "BMI (kg/m2)", "fallback_constant": None, "dtype": "float"},
    }
    result = build_homogeneous_obs(source, "wang25", rules)
    assert result["sex"].tolist() == ["not_provided", "not_provided"]
    assert result["bmi"].isna().all()


def test_optional_source_overrides_fallback_when_present():
    source = pd.DataFrame(
        {"sampleName": ["S1"], "Gender": ["Female"], "BMI (kg/m2)": [21.5]},
        index=["c1"],
    )
    rules = {
        "sex": {"source": "Gender", "fallback_constant": "not_provided", "lower": True},
        "bmi": {"source": "BMI (kg/m2)", "fallback_constant": None, "dtype": "float"},
    }
    result = build_homogeneous_obs(source, "wang25", rules)
    assert result.loc["c1", "sex"] == "female"
    assert result.loc["c1", "bmi"] == 21.5


def test_nullable_strings_are_normalized_before_h5ad_write():
    class AnnDataLike:
        obs = pd.DataFrame({"cell_id": pd.Series(["c1", None], dtype="string")})
        var = pd.DataFrame({"feature": pd.Series(["g1", None], dtype="string")})

    adata = AnnDataLike()
    adata.obs.index = pd.Index(pd.array(["obs1", None], dtype="string"), name="cell_id")
    adata.var.index = pd.Index(pd.array(["var1", None], dtype="string"), name="feature_id")
    _normalize_nullable_strings_for_h5ad(adata)

    assert adata.obs["cell_id"].dtype == object
    assert adata.var["feature"].dtype == object
    assert adata.obs["cell_id"].tolist() == ["c1", None]
    assert adata.var["feature"].tolist() == ["g1", None]
    assert adata.obs.index.dtype == object
    assert adata.var.index.dtype == object
    assert adata.obs.index.tolist() == ["obs1", None]
    assert adata.var.index.tolist() == ["var1", None]


def test_duplicate_gene_symbols_are_aggregated():
    adata = ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        var=pd.DataFrame({"feature_name": ["A", "B", "A"]}),
    )

    repaired, duplicate_count = _repair_features(
        adata, {"source_column": "feature_name", "duplicate_policy": "sum"}
    )

    assert duplicate_count == 2
    assert repaired.var_names.tolist() == ["A", "B"]
    assert repaired.X.toarray().tolist() == [[4, 2], [10, 5]]


def test_missing_metadata_uses_not_provided_for_strings_and_nan_for_floats():
    adata = ad.AnnData(X=np.ones((2, 1)))
    adata.obs["sample"] = pd.Series(["sample_1", None], index=adata.obs_names, dtype="string")
    adata.obs["bmi"] = [20.0, np.nan]

    _normalize_missing_metadata(
        adata,
        {
            "missing_string": "not_provided",
            "required": {"sample": "string", "bmi": "float"},
        },
    )

    assert adata.obs["sample"].tolist() == ["sample_1", "not_provided"]
    assert adata.obs["bmi"].dtype == float
    assert np.isnan(adata.obs.loc[adata.obs_names[1], "bmi"])

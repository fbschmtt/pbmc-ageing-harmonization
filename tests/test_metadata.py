import pandas as pd

from pbmc_pipeline.harmonize import _normalize_nullable_strings_for_h5ad
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

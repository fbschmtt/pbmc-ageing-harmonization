import pandas as pd

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

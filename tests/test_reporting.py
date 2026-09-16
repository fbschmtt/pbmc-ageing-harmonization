import pandas as pd

from pbmc_pipeline.reporting import aifi_l2_concordance, gene_presence_indicators


def test_gene_presence_indicators_uses_merge_provenance_columns():
    var = pd.DataFrame({"available_in_b": [True, False], "available_in_a": [True, True], "other": [1, 2]})
    indicators = gene_presence_indicators(var)
    assert indicators.columns.tolist() == ["a", "b"]
    assert indicators.iloc[1].tolist() == [True, False]


def test_aifi_l2_concordance_normalizes_per_study_label():
    obs = pd.DataFrame({
        "aifi_l2_study_majority": ["T", "T", "B"],
        "aifi_l2_majority": ["T", "B", "B"],
    })
    matrix = aifi_l2_concordance(obs)
    assert matrix.loc["T", "T"] == 0.5
    assert matrix.loc["B", "B"] == 1.0

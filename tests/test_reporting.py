import pandas as pd

from pbmc_pipeline.reporting import (
    aifi_l2_concordance,
    gene_presence_indicators,
    pseudobulk_celltype_fractions,
    sample_cell_type_fractions,
    sample_cluster_fractions,
)


def test_gene_presence_indicators_uses_merge_provenance_columns():
    var = pd.DataFrame({"available_in_b": [True, False], "available_in_a": [True, True], "other": [1, 2]})
    indicators = gene_presence_indicators(var)
    assert indicators.columns.tolist() == ["a", "b"]
    assert indicators.iloc[1].tolist() == [True, False]


def test_aifi_l2_concordance_normalizes_per_study_label():
    obs = pd.DataFrame({
        "aifi_l2_majority": ["T", "T", "B"],
        "experimental_aifi_l2_majority": ["T", "B", "B"],
    })
    matrix = aifi_l2_concordance(obs)
    assert matrix.loc["T", "T"] == 0.5
    assert matrix.loc["B", "B"] == 1.0


def test_pseudobulk_celltype_fractions_weights_each_row_by_cell_count():
    obs = pd.DataFrame({
        "study": ["one", "one", "one"],
        "aifi_l2_majority": ["T", "T", "B"],
        "n_cells": [90, 10, 100],
    })

    fractions = pseudobulk_celltype_fractions(obs)

    assert fractions.loc["one", "T"] == 0.5
    assert fractions.loc["one", "B"] == 0.5


def test_sample_cell_type_fractions_uses_retained_pbmc_denominator():
    adata = type("Adata", (), {})()
    adata.obs = pd.DataFrame({
        "study": ["one", "one", "two"],
        "sample": ["a", "a", "b"],
        "age": [30.0, 30.0, 60.0],
        "n_cells_in_sample": [10, 10, 4],
    })

    fractions = sample_cell_type_fractions(
        adata,
        denominator="n_cells_in_sample",
        fraction_name="fraction_of_retained_pbmc",
    )

    assert fractions["fraction_of_retained_pbmc"].tolist() == [0.2, 0.25]


def test_sample_cluster_fractions_are_within_the_split_cell_type():
    adata = type("Adata", (), {})()
    adata.obs = pd.DataFrame({
        "study": ["one"] * 3,
        "sample": ["a"] * 3,
        "age": [30.0] * 3,
        "cluster": ["0", "0", "1"],
    })

    fractions = sample_cluster_fractions(adata)

    assert fractions.set_index("cluster")["fraction_within_cell_type"].to_dict() == {"0": 2 / 3, "1": 1 / 3}

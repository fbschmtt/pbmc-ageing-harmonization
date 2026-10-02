import pandas as pd

from pbmc_pipeline.covariates import (
    categorical_contrast_label,
    categorical_levels,
    normalize_cmv_status,
)
from pbmc_pipeline.reporting import (
    aifi_l2_concordance,
    gene_presence_indicators,
    pseudobulk_celltype_fractions,
    sample_cell_type_fractions,
    sample_cluster_fractions,
    signed_value_at_largest_absolute_magnitude,
)


def test_categorical_reference_levels_are_explicit_and_cmv_labels_are_harmonized():
    cmv = normalize_cmv_status(pd.Series(["Positive", "Negative", "no", "yes"]))

    assert cmv.tolist() == ["yes", "no", "no", "yes"]
    assert categorical_levels("sex", pd.Series(["male", "female"])) == ["female", "male"]
    assert categorical_levels("cmv", cmv) == ["no", "yes"]
    assert categorical_contrast_label("sex", "male", "female") == "SEX: male vs female"
    assert categorical_contrast_label("cmv", "yes", "no") == "CMV: yes vs no"


def test_signed_largest_magnitude_retains_direction_and_ignores_invalid_values():
    values = pd.Series([1.4, -3.2, 2.7, float("nan"), float("inf"), "invalid"])

    assert signed_value_at_largest_absolute_magnitude(values) == -3.2


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


def test_fraction_helpers_exclude_samples_below_the_minimum_denominator():
    adata = type("Adata", (), {})()
    adata.obs = pd.DataFrame({
        "study": ["one"] * 11,
        "sample": ["included"] * 10 + ["excluded"],
        "age": [30.0] * 10 + [60.0],
        "n_cells_in_sample": [10] * 10 + [1],
        "cluster": ["0"] * 10 + ["1"],
    })

    cell_type = sample_cell_type_fractions(
        adata,
        denominator="n_cells_in_sample",
        fraction_name="fraction",
        min_denominator_cells=10,
    )
    clusters = sample_cluster_fractions(adata, min_denominator_cells=10)

    assert cell_type["sample"].tolist() == ["included"]
    assert clusters["sample"].unique().tolist() == ["included"]

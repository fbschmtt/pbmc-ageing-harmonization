from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from pbmc_pipeline.age_trajectory import prepare_age_trajectory_metadata
from pbmc_pipeline.config import AgeTrajectorySettings, read_json
from pbmc_pipeline.differential_expression import (
    _design_with_covariates,
    _maximal_per_study_covariates,
    _prepare_metadata,
    _run_model,
    _write_age_model_diagnostics,
    combine_cell_type_differential_expression,
    list_pseudobulk_cell_types,
)
from pbmc_pipeline.synthetic_de import (
    _default_samples_per_study,
    create_synthetic_pseudobulk,
)

ROOT = Path(__file__).resolve().parents[1]
TRAJECTORY_SETTINGS = AgeTrajectorySettings.from_mapping(
    read_json(ROOT / "config/pipeline.json")["differential_expression"]["age_trajectory"]
)


def test_age_trajectory_requires_seven_supported_bins_and_excludes_age_90_plus() -> None:
    ages = (
        [24] * 10 + [34] * 10 + [44] * 10 + [54] * 10 + [64] * 10
        + [74] * 10 + [84] * 10 + [94] * 9
    )
    metadata = pd.DataFrame({
        "age": ages,
        "study": ["study_a"] * len(ages),
    })

    retained, support = prepare_age_trajectory_metadata(metadata, TRAJECTORY_SETTINGS)

    assert retained is not None
    assert retained["age_bin"].cat.categories.tolist() == [
        "20-30", "30-40", "40-50", "50-60", "60-70",
        "70-80", "80-90",
    ]
    assert support["dropped_bins"] == []
    assert support["samples_in_dropped_bins"] == 0
    assert support["samples_excluded_by_age_cutoff"] == 9
    assert support["manually_excluded_age_bins"] == ["90-100"]
    assert support["samples_in_manually_excluded_bins"] == 9
    assert support["age_bin_counts_before_filtering"]["90-100"] == 9
    assert support["n_samples_after_age_cutoff"] == 70
    assert support["n_samples_retained"] == 70


def test_age_trajectory_requires_reference_bin_and_minimum_bin_count() -> None:
    metadata = pd.DataFrame({"age": [34] * 10 + [44] * 10, "study": ["s"] * 20})

    retained, support = prepare_age_trajectory_metadata(metadata, TRAJECTORY_SETTINGS)

    assert retained is None
    assert support["status"] == "skipped"
    assert "reference age bin '20-30'" in support["reason"]


def test_test_mode_bypasses_cell_count_cutoff_but_keeps_adult_filter() -> None:
    obs = pd.DataFrame(
        {
            "study": ["study_a", "study_a", "study_a"],
            "study_site": ["site_1", "site_1", "site_1"],
            "sample": ["adult_low_cells", "adult_enough_cells", "under_20"],
            "age": [35, 50, 19],
            "sex": ["female", "male", "female"],
            "n_cells": [2, 12, 2],
            "total_counts": [200, 300, 150],
        },
        index=["pb_1", "pb_2", "pb_3"],
    )
    inclusion = {
        "minimum_age_years_inclusive": 20,
        "minimum_cells_per_pseudobulk_inclusive": 10,
    }

    standard, standard_counts, standard_excluded = _prepare_metadata(obs, inclusion)
    test, test_counts, excluded = _prepare_metadata(obs, inclusion, test_mode=True)

    assert standard["sample"].tolist() == ["adult_enough_cells"]
    assert standard["log10_total_counts"].tolist() == [np.log10(300)]
    assert standard_counts["n_cells_filter_bypassed"] == 0
    assert standard_excluded[-1]["reasons"] == ["age_below_minimum", "n_cells_below_minimum"]
    assert test["sample"].tolist() == ["adult_low_cells", "adult_enough_cells"]
    assert test_counts["n_cells_below_minimum"] == 2
    assert test_counts["n_cells_filter_bypassed"] == 1
    assert test_counts["eligible_before_duplicate_check"] == 2
    assert [item["pseudobulk_id"] for item in excluded] == ["pb_3"]
    assert excluded[0]["reasons"] == ["age_below_minimum"]


def test_synthetic_fixture_is_small_reproducible_and_estimable(tmp_path) -> None:
    first_path = create_synthetic_pseudobulk(
        tmp_path / "first.h5ad", trajectory_settings=TRAJECTORY_SETTINGS
    )
    second_path = create_synthetic_pseudobulk(
        tmp_path / "second.h5ad", trajectory_settings=TRAJECTORY_SETTINGS
    )
    first = ad.read_h5ad(first_path)
    second = ad.read_h5ad(second_path)

    assert first.shape == (168, 96)
    assert first.obs["study"].nunique() == 3
    assert first.obs["study_site"].nunique() == 4
    assert first.obs["sample"].nunique() == 28
    assert first.obs["aifi_l2_majority"].nunique() == 2
    assert first.obs["n_cells"].min() == 15
    assert (first.obs["age"] >= 20).all()
    assert first.obs.groupby(["study", "aifi_l2_majority"], observed=True).size().eq(28).all()
    assert (first.obs["total_counts"] > 0).all()
    age_bin_counts = (
        first.obs.assign(age_bin=(first.obs["age"] // 10 * 10).astype(int))
        .groupby(["aifi_l2_majority", "age_bin"], observed=True)
        .size()
    )
    assert age_bin_counts.eq(12).all()
    synthetic_metadata = first.uns["synthetic_test_data"]
    assert synthetic_metadata["age_trajectory_settings"] == TRAJECTORY_SETTINGS.to_mapping()
    assert len(synthetic_metadata["expected_age_trajectory_profiles"]) == 5
    assert _default_samples_per_study(
        replace(TRAJECTORY_SETTINGS, minimum_samples_per_bin=13), study_count=3
    ) == 35
    assert first.obs.groupby("study", observed=True)["bmi"].apply(
        lambda values: values.notna().any()
    ).to_dict() == {
        "synthetic_study_a": True,
        "synthetic_study_b": False,
        "synthetic_study_c": True,
    }
    assert first.obs.groupby("study", observed=True)["study_site"].nunique().to_dict() == {
        "synthetic_study_a": 2,
        "synthetic_study_b": 1,
        "synthetic_study_c": 1,
    }
    first_study_site_sex_counts = (
        first.obs.loc[first.obs["study"] == "synthetic_study_a"]
        .groupby(["aifi_l2_majority", "study_site", "sex"], observed=True)
        .size()
    )
    expected_samples_per_site_sex = (
        TRAJECTORY_SETTINGS.strict_age_cutoff_exclusive
        - TRAJECTORY_SETTINGS.reference_bin_start_age
    ) // TRAJECTORY_SETTINGS.bin_width_years
    assert first_study_site_sex_counts.eq(expected_samples_per_site_sex).all()
    assert first.obs.groupby("study", observed=True)["cmv"].apply(
        lambda values: values.notna().any()
    ).to_dict() == {
        "synthetic_study_a": True,
        "synthetic_study_b": False,
        "synthetic_study_c": False,
    }
    assert (first.X != second.X).nnz == 0
    assert first.uns["synthetic_test_data"]["not_biological_evidence"] is True


def test_per_study_design_uses_all_varying_covariates_once() -> None:
    metadata = pd.DataFrame({
        "age": [30, 40, 50, 60, 70],
        "sex": ["female", "male", "female", "male", "female"],
        "log10_total_counts": np.log10([100, 120, 140, 160, 180]),
        "bmi": [22.0, 24.0, np.nan, 28.0, 30.0],
        "cmv": ["no", "yes", pd.NA, "yes", "no"],
    })

    covariates = _maximal_per_study_covariates(metadata, ["bmi", "cmv"])

    assert covariates == ["age", "sex", "log10_total_counts", "bmi", "cmv"]
    assert _design_with_covariates(
        "~ study_site + age + sex + log10_total_counts", covariates
    ) == "~ study_site + age + sex + log10_total_counts + bmi + cmv"


def test_combined_covariate_fit_is_allowed_with_a_single_study(monkeypatch, tmp_path) -> None:
    import pbmc_pipeline.differential_expression as de

    index = [f"sample_{number}" for number in range(8)]
    obs = pd.DataFrame(
        {
            "study": ["aifi"] * 8,
            "study_site": ["aifi"] * 8,
            "sample": index,
            "age": [24, 31, 38, 45, 52, 60, 68, 76],
            "sex": ["female", "male"] * 4,
            "cmv": ["no", "yes", "no", pd.NA, "no", "yes", "no", "yes"],
            "n_cells": [20] * 8,
            "total_counts": [200 + 10 * number for number in range(8)],
            "log10_total_counts": np.log10([200 + 10 * number for number in range(8)]),
        },
        index=index,
    )
    var = pd.DataFrame({"available_in_aifi": [True, True]}, index=["G1", "G2"])
    counts = np.array([[100 + 5 * number, 100 + 5 * number] for number in range(8)])
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    captured = {}

    def fake_fit(counts, genes, metadata, **kwargs):
        captured["design"] = kwargs["design"]
        captured["metadata_index"] = metadata.index.tolist()
        return {
            "cmv": pd.DataFrame({
                "gene": ["G1", "G2"],
                "padj": [0.01, 0.9],
                "log2FoldChange": [1.2, 0.1],
                "contrast": ["yes vs no"] * 2,
            })
        }

    monkeypatch.setattr(de, "_fit_covariate_model", fake_fit)
    result = _run_model(
        adata=adata,
        obs=obs,
        metadata=obs.copy(),
        cell_type="CD14 monocyte",
        model_name="combined",
        study=None,
        design="~ study_site + age + sex + log10_total_counts + cmv",
        covariates=["cmv"],
        alpha=0.05,
        cpus=1,
        output_dir=tmp_path,
        test_mode=False,
        synthetic_test_data=True,
    )

    assert result["status"] == "complete"
    assert result["studies"] == ["aifi"]
    assert result["study_sample_counts"] == {"aifi": 7}
    assert result["design"] == "~ age + sex + log10_total_counts + cmv"
    assert captured["design"] == "~ age + sex + log10_total_counts + cmv"
    assert result["study_sites"] == ["aifi"]
    assert result["study_site_sample_counts"] == {"aifi": 7}
    assert captured["metadata_index"] == [index[position] for position in (0, 1, 2, 4, 5, 6, 7)]
    assert result["n_samples_before_complete_case_filter"] == 8
    assert result["n_samples_excluded_missing_design_covariates"] == 1
    assert result["missing_values_by_design_covariate"] == {
        "study_site": 0, "age": 0, "sex": 0, "log10_total_counts": 0, "cmv": 1,
    }


def test_age_model_diagnostics_record_the_shared_gene_universe_and_raw_counts(tmp_path) -> None:
    obs = pd.DataFrame(
        {
            "study": ["one", "two"],
            "study_site": ["site_one", "site_two"],
            "sample": ["one_sample", "two_sample"],
            "age": [30.0, 60.0],
            "sex": ["female", "male"],
            "n_cells": [20, 25],
            "total_counts": [10, 31],
            "log10_total_counts": np.log10([10, 31]),
        },
        index=["one::one_sample::T", "two::two_sample::T"],
    )
    adata = ad.AnnData(
        X=np.array([[2, 3, 5], [7, 11, 13]], dtype=np.int64),
        obs=obs,
        var=pd.DataFrame(
            {"available_in_one": [True, True, False], "available_in_two": [True, False, True]},
            index=["shared", "one_only", "two_only"],
        ),
    )

    output = _write_age_model_diagnostics(
        adata=adata,
        obs=obs,
        metadata=obs.copy(),
        cell_type="T",
        output_dir=tmp_path,
    )

    assert output == tmp_path / "t" / "age_model_diagnostics.csv"
    diagnostics = pd.read_csv(output)
    assert diagnostics["genes_in_age_intersection"].tolist() == [1, 1]
    assert diagnostics["counts_in_age_gene_intersection"].tolist() == [2, 7]
    assert diagnostics["total_counts"].tolist() == [10, 31]


def test_synthetic_fixture_requires_production_cell_cutoff(tmp_path) -> None:
    import pytest

    with pytest.raises(ValueError, match="configured DE cell cutoff"):
        create_synthetic_pseudobulk(
            tmp_path / "too_small.h5ad",
            trajectory_settings=TRAJECTORY_SETTINGS,
            cells_per_pseudobulk=9,
        )


def test_parallel_cell_type_results_rebuild_the_standard_de_manifest(tmp_path) -> None:
    input_path = create_synthetic_pseudobulk(
        tmp_path / "pseudobulk.h5ad", trajectory_settings=TRAJECTORY_SETTINGS
    )
    pipeline = read_json(ROOT / "config" / "pipeline.json")
    cell_types = list_pseudobulk_cell_types(input_path, pipeline)
    assert len(cell_types) == 2
    assert [slug for slug, _ in cell_types] == ["cd14-monocyte", "naive-cd4-t-cell"]

    task_results = []
    for slug, cell_type in cell_types:
        result_dir = tmp_path / "task_results" / slug
        result_dir.mkdir(parents=True)
        (result_dir / "cell_type_result.json").write_text(json.dumps({
            "cell_type": cell_type,
            "slug": slug,
            "n_pseudobulks": 24,
            "metadata_exclusions": {},
            "excluded_samples": [],
            "models": [],
        }))
        task_results.append(result_dir)

    report = combine_cell_type_differential_expression(
        input_path, tmp_path / "differential_expression", pipeline, task_results,
        synthetic_test_data=True,
    )

    assert report["status"] == "no_estimable_models"
    assert list(report["results_by_cell_type"]) == [slug for slug, _ in cell_types]
    assert (tmp_path / "differential_expression" / cell_types[0][0] / "cell_type_result.json").is_file()

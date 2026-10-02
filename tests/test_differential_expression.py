from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd

from pbmc_pipeline.config import read_json
from pbmc_pipeline.differential_expression import (
    _prepare_metadata,
    _run_model,
    combine_cell_type_differential_expression,
    list_pseudobulk_cell_types,
)
from pbmc_pipeline.synthetic_de import create_synthetic_pseudobulk

ROOT = Path(__file__).resolve().parents[1]


def test_test_mode_bypasses_cell_count_cutoff_but_keeps_adult_filter() -> None:
    obs = pd.DataFrame(
        {
            "study": ["study_a", "study_a", "study_a"],
            "sample": ["adult_low_cells", "adult_enough_cells", "under_20"],
            "age": [35, 50, 19],
            "sex": ["female", "male", "female"],
            "n_cells": [2, 12, 2],
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
    assert standard_counts["n_cells_filter_bypassed"] == 0
    assert standard_excluded[-1]["reasons"] == ["age_below_minimum", "n_cells_below_minimum"]
    assert test["sample"].tolist() == ["adult_low_cells", "adult_enough_cells"]
    assert test_counts["n_cells_below_minimum"] == 2
    assert test_counts["n_cells_filter_bypassed"] == 1
    assert test_counts["eligible_before_duplicate_check"] == 2
    assert [item["pseudobulk_id"] for item in excluded] == ["pb_3"]
    assert excluded[0]["reasons"] == ["age_below_minimum"]


def test_synthetic_fixture_is_small_reproducible_and_estimable(tmp_path) -> None:
    first_path = create_synthetic_pseudobulk(tmp_path / "first.h5ad")
    second_path = create_synthetic_pseudobulk(tmp_path / "second.h5ad")
    first = ad.read_h5ad(first_path)
    second = ad.read_h5ad(second_path)

    assert first.shape == (48, 91)
    assert first.obs["study"].nunique() == 3
    assert first.obs["sample"].nunique() == 8
    assert first.obs["aifi_l2_majority"].nunique() == 2
    assert first.obs["n_cells"].min() == 15
    assert (first.obs["age"] >= 20).all()
    assert first.obs.groupby(["study", "aifi_l2_majority"], observed=True).size().eq(8).all()
    assert first.obs.groupby("study", observed=True)["bmi"].apply(
        lambda values: values.notna().any()
    ).to_dict() == {
        "synthetic_study_a": True,
        "synthetic_study_b": False,
        "synthetic_study_c": True,
    }
    assert first.obs.groupby("study", observed=True)["cmv"].apply(
        lambda values: values.notna().any()
    ).to_dict() == {
        "synthetic_study_a": True,
        "synthetic_study_b": False,
        "synthetic_study_c": False,
    }
    assert (first.X != second.X).nnz == 0
    assert first.uns["synthetic_test_data"]["not_biological_evidence"] is True


def test_combined_covariate_fit_is_allowed_with_a_single_study(monkeypatch, tmp_path) -> None:
    import numpy as np

    import pbmc_pipeline.differential_expression as de

    index = [f"sample_{number}" for number in range(8)]
    obs = pd.DataFrame(
        {
            "study": ["aifi"] * 8,
            "sample": index,
            "age": [24, 31, 38, 45, 52, 60, 68, 76],
            "sex": ["female", "male"] * 4,
            "cmv": ["no", "yes"] * 4,
            "n_cells": [20] * 8,
        },
        index=index,
    )
    var = pd.DataFrame({"available_in_aifi": [True, True]}, index=["G1", "G2"])
    adata = ad.AnnData(X=np.full((8, 2), 100, dtype=np.int32), obs=obs, var=var)
    captured = {}

    def fake_fit(counts, genes, metadata, **kwargs):
        captured["design"] = kwargs["design"]
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
        design="~ study + age + sex + cmv",
        covariates=["cmv"],
        alpha=0.05,
        cpus=1,
        output_dir=tmp_path,
        test_mode=False,
        synthetic_test_data=True,
    )

    assert result["status"] == "complete"
    assert result["studies"] == ["aifi"]
    assert result["study_sample_counts"] == {"aifi": 8}
    assert result["design"] == "~ age + sex + cmv"
    assert captured["design"] == "~ age + sex + cmv"


def test_synthetic_fixture_requires_production_cell_cutoff(tmp_path) -> None:
    import pytest

    with pytest.raises(ValueError, match="configured DE cell cutoff"):
        create_synthetic_pseudobulk(tmp_path / "too_small.h5ad", cells_per_pseudobulk=9)


def test_parallel_cell_type_results_rebuild_the_standard_de_manifest(tmp_path) -> None:
    input_path = create_synthetic_pseudobulk(tmp_path / "pseudobulk.h5ad")
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

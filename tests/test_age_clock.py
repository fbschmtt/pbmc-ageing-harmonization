import json
import os
import runpy
import subprocess
import sys
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from pbmc_pipeline.age_clock import (
    ClockFitFailuresError,
    ClockMethod,
    WithinSampleRankTransformer,
    _clock_methods,
    _coefficient_coverage,
    _gene_coverage,
    _log_cpm,
    _metadata_for_type,
    _sample_weights,
    _training_gene_intersection_diagnostics,
    _training_gene_names,
    predict_age_clock,
    train_age_clocks,
    train_all_age_clocks,
)


class PseudobulkVar:
    def __init__(self, var: pd.DataFrame):
        self.var = var


def test_study_balanced_weights_give_each_study_equal_total_weight():
    methods = _clock_methods([0.1], [0.5])
    balanced = next(method for method in methods if method.equal_study_weight)
    weights = _sample_weights(pd.Series(["large"] * 3 + ["small"]), balanced)

    assert weights is not None
    assert weights.mean() == pytest.approx(1.0)
    assert weights[:3].sum() == pytest.approx(weights[3:].sum())


def test_sample_weighted_method_leaves_each_sample_equal():
    method = _clock_methods([0.1], [0.5])[0]

    assert _sample_weights(pd.Series(["large", "large", "small"]), method) is None


def test_rank_transform_uses_within_sample_percentiles_and_is_scale_invariant():
    transformer = WithinSampleRankTransformer().fit([[1, 2, 4, 8]])
    original = transformer.transform([[1, 2, 4, 8]])
    rescaled = transformer.transform([[10, 20, 40, 80]])

    assert original.tolist() == [[0.0, 1 / 3, 2 / 3, 1.0]]
    assert rescaled.tolist() == original.tolist()


def test_rank_clock_uses_sample_weighting_and_rank_pipeline():
    rank_method = next(
        method for method in _clock_methods([0.1], [0.5])
        if method.name == "rank_elastic_net"
    )

    assert _sample_weights(pd.Series(["large", "large", "small"]), rank_method) is None
    assert "rank" in rank_method.estimator_factory().named_steps


def test_gene_intersection_uses_only_the_requested_training_studies():
    pseudobulk = PseudobulkVar(pd.DataFrame({
        "available_in_train_a": [True, True, False],
        "available_in_train_b": [True, False, True],
        "available_in_heldout": [False, True, True],
    }, index=["shared", "a_only", "b_only"]))

    assert _training_gene_names(pseudobulk, ["train_a", "train_b"]) == ["shared"]


def test_gene_intersection_requires_training_study_availability_flags():
    pseudobulk = PseudobulkVar(pd.DataFrame({
        "available_in_train_a": [True],
    }, index=["gene"]))

    with pytest.raises(ValueError, match="availability is incomplete"):
        _training_gene_names(pseudobulk, ["train_a", "train_b"])


def test_training_gene_intersection_reports_genes_dropped_per_study():
    pseudobulk = PseudobulkVar(pd.DataFrame({
        "available_in_train_a": [True, True, False],
        "available_in_train_b": [True, False, True],
    }, index=["shared", "a_only", "b_only"]))

    diagnostic = _training_gene_intersection_diagnostics(
        pseudobulk, ["train_a", "train_b"], ["shared"]
    )

    assert diagnostic["n_genes_available_in_any_training_study"] == 3
    assert diagnostic["n_genes_in_final_model_universe"] == 1
    assert diagnostic["n_genes_dropped_from_training_union"] == 2
    assert diagnostic["by_study"] == [
        {
            "study": "train_a", "n_genes_available_before_intersection": 2,
            "n_genes_dropped_from_study_universe": 1,
        },
        {
            "study": "train_b", "n_genes_available_before_intersection": 2,
            "n_genes_dropped_from_study_universe": 1,
        },
    ]


@pytest.fixture
def clock_input(tmp_path):
    studies = np.repeat(["train_a", "train_b", "external_a", "external_b"], 4)
    rng = np.random.default_rng(12)
    matrix = rng.poisson(30, (len(studies), 6))
    obs = pd.DataFrame({
        "study": studies,
        "sample": [f"sample_{i}" for i in range(len(studies))],
        "age": np.tile([25.0, 40.0, 55.0, 70.0], 4),
        "n_cells": 15,
        "total_counts": matrix.sum(axis=1),
        "aifi_l2_majority": "CD14 monocyte",
    }, index=[f"bulk_{i}" for i in range(len(studies))])
    var = pd.DataFrame(
        {f"available_in_{study}": True for study in np.unique(studies)},
        index=[f"gene_{i}" for i in range(matrix.shape[1])],
    )
    adata = ad.AnnData(sparse.csr_matrix(matrix), obs=obs, var=var)
    config = {"differential_expression": {
        "split_by": "aifi_l2_majority",
        "sample_inclusion": {
            "minimum_age_years_inclusive": 20,
            "minimum_cells_per_pseudobulk_inclusive": 10,
        },
    }}
    input_path = tmp_path / "pseudobulk.h5ad"
    config_path = tmp_path / "config.json"
    adata.write_h5ad(input_path)
    config_path.write_text(json.dumps(config))
    return adata, config, input_path, config_path


def _train_fixture(clock_input, output_dir, **kwargs):
    _, _, input_path, config_path = clock_input
    return train_age_clocks(
        input_path, config_path, output_dir, ["external_a", "external_b"],
        cell_type_filter="cd14-monocyte", alphas=[0.1], l1_ratios=[0.5], **kwargs,
    )


def test_all_cell_type_runner_discovers_and_fits_available_types(clock_input, tmp_path):
    _, _, input_path, config_path = clock_input

    manifests = train_all_age_clocks(
        input_path, config_path, tmp_path / "clock", ["external_a", "external_b"],
        alphas=[0.1], l1_ratios=[0.5], cv=False,
    )

    assert list(manifests) == ["CD14 monocyte"]
    assert (
        tmp_path / "clock/cd14-monocyte/ageing_clock/elastic_net/model.joblib"
    ).is_file()


def test_all_cell_type_runner_skips_type_missing_from_studies(clock_input, tmp_path, caplog):
    adata, _, input_path, config_path = clock_input
    adata.obs["aifi_l2_majority"] = adata.obs["aifi_l2_majority"].astype(str)
    adata.obs.loc["bulk_8", "aifi_l2_majority"] = "Rare type"
    adata.write_h5ad(input_path)

    with caplog.at_level("ERROR"):
        manifests = train_all_age_clocks(
            input_path, config_path, tmp_path / "clock", ["external_a", "external_b"],
            alphas=[0.1], l1_ratios=[0.5], cv=False,
        )

    assert list(manifests) == ["CD14 monocyte"]
    assert "Skipping age clock for Rare type" in caplog.text
    assert "Training studies have no eligible samples" in caplog.text
    assert "External studies have no eligible samples" in caplog.text
    summary = json.loads((tmp_path / "clock/age_clock_all_types_summary.json").read_text())
    assert summary["successful_cell_types"] == ["CD14 monocyte"]
    assert summary["eligibility_skipped_cell_types"][0]["cell_type"] == "Rare type"
    assert summary["failed_cell_types"] == []


def test_all_cell_type_runner_continues_after_unexpected_fit_error(
    clock_input, tmp_path, monkeypatch, caplog,
):
    from pbmc_pipeline import age_clock

    adata, _, input_path, config_path = clock_input
    adata.obs["aifi_l2_majority"] = adata.obs["aifi_l2_majority"].astype(str)
    adata.obs.loc["bulk_8", "aifi_l2_majority"] = "Rare type"
    adata.write_h5ad(input_path)
    fit_age_clock = age_clock.fit_age_clocks

    def fail_rare_type(*args, cell_type_filter=None, **kwargs):
        if cell_type_filter == "Rare type":
            raise RuntimeError("deliberate fit failure")
        return fit_age_clock(*args, cell_type_filter=cell_type_filter, **kwargs)

    monkeypatch.setattr(age_clock, "fit_age_clocks", fail_rare_type)
    with caplog.at_level("ERROR"), pytest.raises(ClockFitFailuresError, match="1 cell type"):
        train_all_age_clocks(
            input_path, config_path, tmp_path / "clock", ["external_a", "external_b"],
            alphas=[0.1], l1_ratios=[0.5], cv=False,
        )

    assert (tmp_path / "clock/cd14-monocyte/ageing_clock/age_clock.json").is_file()
    assert "deliberate fit failure" in caplog.text
    summary = json.loads((tmp_path / "clock/age_clock_all_types_summary.json").read_text())
    assert summary["successful_cell_types"] == ["CD14 monocyte"]
    assert summary["failed_cell_types"][0]["cell_type"] == "Rare type"
    assert "RuntimeError: deliberate fit failure" in summary["failed_cell_types"][0]["traceback"]


def test_all_cell_type_runner_does_not_catch_publication_errors(clock_input, tmp_path, monkeypatch):
    from pbmc_pipeline import age_clock

    monkeypatch.setattr(
        age_clock,
        "publish_age_clock_fit",
        lambda *_: (_ for _ in ()).throw(OSError("deliberate publication failure")),
    )
    _, _, input_path, config_path = clock_input
    with pytest.raises(OSError, match="deliberate publication failure"):
        train_all_age_clocks(
            input_path, config_path, tmp_path / "clock", ["external_a", "external_b"],
            alphas=[0.1], l1_ratios=[0.5], cv=False,
        )


def test_module_cli_rank_model_loads_in_a_separate_process_context(clock_input, tmp_path):
    adata, _, input_path, config_path = clock_input
    output = tmp_path / "clock"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    subprocess.run([
        sys.executable, "-m", "pbmc_pipeline.age_clock", "--input", str(input_path),
        "--config", str(config_path), "--output-dir", str(output),
        "--external-studies", "external_a,external_b", "--cell-type", "cd14-monocyte",
        "--no-cv",
    ], env=env, check=True, capture_output=True, text=True)
    model = output / "cd14-monocyte/ageing_clock/rank_elastic_net/model.joblib"
    bundle = joblib.load(model)
    assert bundle["software_versions"]["scikit-learn"]
    assert bundle["configuration_sha256"]
    transformer = bundle["estimator"].named_steps["rank"]
    assert type(transformer).__module__ == "pbmc_pipeline.age_clock_features"
    predictions = predict_age_clock(model, adata, np.arange(adata.n_obs))
    assert predictions.shape == (adata.n_obs,)
    assert np.isfinite(predictions).all()


@pytest.mark.parametrize("reason", ["missing_cell_type", "below_cell_cutoff"])
def test_every_requested_external_study_needs_eligible_samples(clock_input, tmp_path, reason):
    adata, _, input_path, _ = clock_input
    external_b = adata.obs["study"] == "external_b"
    if reason == "missing_cell_type":
        adata.obs["aifi_l2_majority"] = adata.obs["aifi_l2_majority"].astype(str)
        adata.obs.loc[external_b, "aifi_l2_majority"] = "Other type"
    else:
        adata.obs.loc[external_b, "n_cells"] = 9
    adata.write_h5ad(input_path)
    with pytest.raises(ValueError, match="External studies have no eligible samples.*external_b") as error:
        _train_fixture(clock_input, tmp_path / "clock")
    assert "Per-study inclusion counts" in str(error.value)
    assert not list((tmp_path / "clock").rglob("model.joblib"))


def test_normalization_uses_model_genes_and_ignores_cached_totals(clock_input):
    adata, _, _, _ = clock_input
    rows = np.arange(adata.n_obs)
    genes = ["gene_1", "gene_0"]
    expected = np.log1p(
        adata[:, genes].X.toarray() * 1_000_000 / adata[:, genes].X.sum(axis=1)
    )
    np.testing.assert_allclose(_log_cpm(adata, rows, genes).toarray(), expected)
    adata.obs["total_counts"] = -123
    np.testing.assert_allclose(_log_cpm(adata, rows, genes).toarray(), expected)
    del adata.obs["total_counts"]
    np.testing.assert_allclose(_log_cpm(adata, rows, genes).toarray(), expected)
    with_extra_genes = adata.copy()
    expanded = with_extra_genes.X.toarray()
    expanded[:, 2:] *= 100
    with_extra_genes.X = sparse.csr_matrix(expanded)
    np.testing.assert_allclose(
        _log_cpm(with_extra_genes, rows, genes).toarray(), expected
    )


def test_inclusion_uses_actual_library_sizes_and_finite_ages(clock_input):
    adata, config, _, _ = clock_input
    counts = adata.X.tolil()
    counts[0, :] = 0
    adata.X = counts.tocsr()
    adata.obs.loc["bulk_1", "age"] = np.inf
    metadata, counts = _metadata_for_type(adata, config, "CD14 monocyte")
    assert "bulk_0" not in metadata.index
    assert "bulk_1" not in metadata.index
    assert counts["by_study"]["train_a"]["excluded_invalid_or_zero_library"] == 1
    assert counts["by_study"]["train_a"]["excluded_invalid_or_below_minimum_age"] == 1


def test_partial_gene_coverage_aligns_and_zero_fills_but_zero_overlap_fails(clock_input):
    adata, _, _, _ = clock_input
    partial = adata[:, ["gene_0", "gene_2"]].copy()
    rows = np.arange(partial.n_obs)
    coverage = _gene_coverage(partial, rows, ["gene_0", "gene_1"])
    assert all(record["fraction_available"] == 0.5 for record in coverage)
    normalized = _log_cpm(partial, rows, ["gene_0", "gene_1"]).toarray()
    assert np.all(normalized[:, 1] == 0)
    with pytest.raises(ValueError, match="No model genes are available"):
        _log_cpm(adata[:, ["gene_2"]], rows, ["gene_0", "gene_1"])
    adata.var["available_in_external_b"] = False
    with pytest.raises(ValueError, match="No model genes are available for study external_b"):
        _log_cpm(adata, rows, ["gene_0", "gene_1"])


def test_missing_coefficient_weight_uses_absolute_coefficient_mass(clock_input):
    adata, _, _, _ = clock_input
    adata.var.loc["gene_1", "available_in_external_a"] = False
    coverage = _coefficient_coverage(
        adata, ["gene_0", "gene_1", "gene_2"], np.array([2.0, -1.0, 7.0]),
        "external_a",
    )

    assert coverage["missing_coefficient_absolute_weight"] == pytest.approx(1.0)
    assert coverage["total_coefficient_absolute_weight"] == pytest.approx(10.0)
    assert coverage["missing_coefficient_weight_fraction"] == pytest.approx(0.1)
    assert coverage["missing_coefficient_weight_warning"] is True


def test_no_cv_rerun_removes_previous_cv_artifacts(clock_input, tmp_path):
    output = tmp_path / "clock"
    manifest = _train_fixture(clock_input, output)
    clock_root = output / "cd14-monocyte/ageing_clock"
    for method in manifest["methods"]:
        assert (clock_root / method["method"] / "cv_oof_predictions.csv").is_file()
    manifest = _train_fixture(clock_input, output, cv=False)
    assert manifest["cross_validation"] == "disabled"
    assert not list(clock_root.glob("*/cv_*.csv"))


def test_cv_candidate_parallelism_records_convergence_diagnostics(clock_input, tmp_path):
    _, _, input_path, config_path = clock_input
    manifest = train_age_clocks(
        input_path, config_path, tmp_path / "clock", ["external_a", "external_b"],
        cell_type_filter="cd14-monocyte", alphas=[0.1, 1.0], l1_ratios=[0.5], n_jobs=2,
    )

    assert manifest["cv_parallel_jobs"] == 2
    for method in manifest["methods"]:
        assert all("n_iterations" in fold for fold in method["cv_mae_by_held_out_study"])
        grid = pd.read_csv(Path(method["model"]).parent / "cv_grid_search.csv")
        assert {"n_iterations_train_a", "n_iterations_train_b"}.issubset(grid.columns)


def test_external_expression_age_and_availability_do_not_change_fits(clock_input, tmp_path):
    adata, _, input_path, _ = clock_input
    first = _train_fixture(clock_input, tmp_path / "first")
    external = adata.obs["study"].astype(str).str.startswith("external")
    adata.X[external.to_numpy(), :] *= 7
    adata.obs.loc[external, "age"] = 95.0 - adata.obs.loc[external, "age"]
    adata.var.loc["gene_0", "available_in_external_a"] = False
    adata.write_h5ad(input_path)
    second = _train_fixture(clock_input, tmp_path / "second")
    for first_method, second_method in zip(first["methods"], second["methods"]):
        first_bundle = joblib.load(first_method["model"])
        second_bundle = joblib.load(second_method["model"])
        assert first_bundle["genes"] == second_bundle["genes"]
        assert first_method["best_parameters"] == second_method["best_parameters"]
        for step, attribute in [("scale", "mean_"), ("scale", "scale_"), ("regression", "coef_")]:
            np.testing.assert_array_equal(
                getattr(first_bundle["estimator"].named_steps[step], attribute),
                getattr(second_bundle["estimator"].named_steps[step], attribute),
            )


@pytest.mark.parametrize("cv", [True, False])
def test_generic_method_parameters_and_report_templates(clock_input, tmp_path, monkeypatch, cv):
    method = ClockMethod(
        name="ridge", display_name="Ridge", description="Ridge regression.",
        estimator_factory=lambda: Pipeline([
            ("scale", StandardScaler(with_mean=False)), ("regression", Ridge()),
        ]),
        parameter_grid={"regression__alpha": [1.0], "regression__max_iter": [50]},
        default_parameters={"regression__alpha": 1.0},
    )
    monkeypatch.setattr("pbmc_pipeline.age_clock._clock_methods", lambda *_: [method])
    manifest = _train_fixture(clock_input, tmp_path / "clock", cv=cv)
    if cv:
        assert isinstance(manifest["methods"][0]["best_parameters"]["regression__max_iter"], int)
    import matplotlib.pyplot as plt

    monkeypatch.setattr("IPython.display.display", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        "pbmc_pipeline.report_theme.display_accessible_figure",
        lambda figure, **kwargs: plt.close(figure),
    )
    report_root = tmp_path / "clock/cd14-monocyte/ageing_clock"
    if not cv:
        for filename in ("cv_oof_predictions.csv", "cv_grid_search.csv"):
            (report_root / "ridge" / filename).write_text("stale file from an old run")
    project_root = Path(__file__).resolve().parents[1]
    monkeypatch.setenv("AGE_CLOCK_REPORT_DIR", str(report_root / "ridge"))
    runpy.run_path(str(project_root / "reports/age_clock_report.py"))
    monkeypatch.setenv("AGE_CLOCK_REPORT_DIR", str(report_root))
    runpy.run_path(str(project_root / "reports/age_clock_comparison_report.py"))

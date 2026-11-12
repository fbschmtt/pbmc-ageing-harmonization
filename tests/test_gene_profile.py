from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from formulaic_contrasts import FormulaicContrasts

from pbmc_pipeline.gene_profile_report import (
    format_model_results,
    generate_gene_profile_report,
)
from scripts.profile_gene import (
    _annotate_study_slope_coefficients,
    _collect_cross_cell_type_trajectories,
    _exploratory_infant_adult_contrast,
    _fit_age_bin_profile,
    _fit_profile,
    _fit_study_age_bin_trajectories,
    _fit_study_slope_profile,
    _prepare_model_for_cache,
)

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = ROOT / "config" / "pipeline.json"


def test_study_slope_coefficient_annotations_match_reported_contrasts() -> None:
    coefficients = pd.DataFrame({
        "term": ["Intercept", "age", "age:study[T.study_b]"],
        "coefficient": [1.0, 0.2, -0.05],
    })
    slopes = pd.DataFrame({"study": ["study_a", "study_b"], "slope": [0.2, 0.15]})

    annotated = _annotate_study_slope_coefficients(
        coefficients, slopes, ["study_a", "study_b"],
    )

    assert annotated.loc[annotated["term"] == "age", "parameter_role"].item() == (
        "Per-year slope for the reference study"
    )
    interaction = annotated.loc[annotated["term"] == "age:study[T.study_b]"].iloc[0]
    assert interaction["parameter_role"] == "Difference from the reference-study per-year slope"
    assert "not a standalone slope" in interaction["interpretation"]


def test_cross_cell_type_trajectories_read_completed_age_bin_fits(tmp_path) -> None:
    de_dir = tmp_path / "per_type"
    output_dir = de_dir / "cd14-monocyte" / "differential_expression"
    (output_dir / "combined").mkdir(parents=True)
    (output_dir / "run_metadata.json").write_text(json.dumps({
        "cell_type": "CD14 monocyte",
        "models": [{
            "purpose": "age_trajectory", "status": "complete",
            "results_by_covariate": {"age_bin": "combined/age_bin.csv"},
        }],
    }))
    pd.DataFrame({
        "cell_type": ["stale source value", "stale source value", "stale source value"],
        "gene": ["GENE1", "GENE1", "GENE2"],
        "age_bin": ["20-30", "30-40", "20-30"],
        "log2FoldChange": [0.0, 0.2, 0.0],
        "omnibus_padj": [0.01, 0.01, 0.8],
    }).to_csv(output_dir / "combined" / "age_bin.csv", index=False)

    trajectories, support = _collect_cross_cell_type_trajectories(de_dir, "GENE1")

    assert trajectories["age_bin"].tolist() == ["20-30", "30-40"]
    assert trajectories["cell_type"].tolist() == ["CD14 monocyte", "CD14 monocyte"]
    assert support["status"].tolist() == ["complete"]


def test_study_age_bin_trajectories_align_missing_reference_and_keep_gaps() -> None:
    metadata = pd.DataFrame({
        "study": ["reference_a"] * 10 + ["reference_b"] * 10 + ["missing_reference"] * 10,
        "age_bin": (
            ["20-30"] * 5 + ["30-40"] * 5
            + ["20-30"] * 5 + ["30-40"] * 5
            + ["30-40"] * 5 + ["50-60"] * 5
        ),
    })
    expression = np.array(
        [2.0] * 5 + [3.0] * 5 + [4.0] * 5 + [6.0] * 5
        + [10.0] * 5 + [12.0] * 5
    )

    trajectories, support = _fit_study_age_bin_trajectories(
        metadata, expression, minimum_samples_per_bin=5,
        reference_bin_start_age=20, bin_width=10,
    )

    missing_reference = trajectories.loc[
        trajectories["study"] == "missing_reference"
    ].set_index("age_bin")
    assert missing_reference.loc["30-40", "aligned_expression"] == 3.0
    assert missing_reference.loc["50-60", "aligned_expression"] == 5.0
    assert (missing_reference["alignment_bin"] == "30-40").all()
    reference_baselines = trajectories.loc[
        trajectories["age_bin"] == "20-30", "aligned_expression",
    ]
    assert reference_baselines.tolist() == [3.0, 3.0]
    missing_gap = support.loc[
        (support["study"] == "missing_reference") & (support["age_bin"] == "40-50")
    ].iloc[0]
    assert missing_gap["status"] == "missing"
    assert missing_gap["n_samples"] == 0


def test_study_age_bin_trajectories_exclude_bins_below_cutoff_and_handle_no_reference() -> None:
    metadata = pd.DataFrame({
        "study": ["study_a"] * 4 + ["study_b"] * 5,
        "age_bin": ["20-30"] * 4 + ["30-40"] * 5,
    })

    trajectories, support = _fit_study_age_bin_trajectories(
        metadata, np.arange(9, dtype=float), minimum_samples_per_bin=5,
        reference_bin_start_age=20, bin_width=10,
    )

    assert trajectories.empty
    assert set(support.loc[support["status"] == "below minimum", "study"]) == {"study_a"}
    assert set(support.loc[support["status"] == "unaligned", "study"]) == {"study_b"}


def test_profile_model_cache_drops_series_metadata_before_h5ad_write(tmp_path) -> None:
    obs = pd.DataFrame({"status": pd.Series(["yes", "no"], dtype="string")})
    model = ad.AnnData(
        X=np.array([[1], [2]], dtype=np.int64),
        obs=obs,
        var=pd.DataFrame(index=["GENE1"]),
    )
    model.layers["fitted"] = np.array([[1.1], [1.9]])
    model.obsm["design_matrix"] = np.array([[1.0], [0.0]])
    model.uns["trend_coeffs"] = pd.Series([1.0, 2.0])
    model.uns["vst_trend_coeffs"] = pd.Series([3.0, 4.0])

    _prepare_model_for_cache(model)
    cache_path = tmp_path / "fitted_model.h5ad"
    model.write_h5ad(cache_path)
    restored = ad.read_h5ad(cache_path)

    assert restored.uns == {}
    assert restored.layers["fitted"].tolist() == [[1.1], [1.9]]
    assert restored.obsm["design_matrix"].tolist() == [[1.0], [0.0]]
    assert restored.obs["status"].tolist() == ["yes", "no"]


@pytest.mark.parametrize('dtype', ['string', object])
@pytest.mark.parametrize('values', [['no', None], [None, None]])
def test_profile_cache_preserves_missing_text_metadata(tmp_path, dtype, values) -> None:
    obs = pd.DataFrame({
        'smoking_status': pd.Series(values, dtype=dtype),
        'age': [30.0, 50.0],
        'size_factors': [0.8, 1.2],
    })
    model = ad.AnnData(X=np.array([[2], [3]]), obs=obs)
    missing = model.obs['smoking_status'].isna().tolist()

    _prepare_model_for_cache(model)
    model.write_h5ad(tmp_path / 'model.h5ad')
    restored = ad.read_h5ad(tmp_path / 'model.h5ad')

    assert restored.obs['smoking_status'].isna().tolist() == missing
    assert restored.obs['smoking_status'].dropna().tolist() == obs['smoking_status'].dropna().tolist()
    assert restored.obs['age'].tolist() == [30.0, 50.0]
    assert restored.obs['size_factors'].tolist() == [0.8, 1.2]
    np.testing.assert_array_equal(restored.X, [[2], [3]])


def age_bin_profile_inputs(ages, sexes, *, counts=None):
    pipeline = json.loads(PIPELINE_PATH.read_text())
    spec = pipeline['differential_expression']
    spec['study_comparison_minimum_samples_per_bin'] = 2
    obs = pd.DataFrame({
        'study': 'study_a', 'study_site': 'site_a',
        'sample': [f'sample_{i}' for i in range(len(ages))],
        'age': ages, 'sex': sexes, 'total_counts': 100.0, 'n_cells': 20,
    }, index=[f'pb_{i}' for i in range(len(ages))])
    values = counts if counts is not None else np.array([[i + 2, i + 3] for i in range(len(ages))])
    adata = ad.AnnData(
        X=values, obs=obs,
        var=pd.DataFrame({'available_in_study_a': True}, index=['GENE1', 'GENE2']),
    )
    return SimpleNamespace(
        pipeline=pipeline, spec=spec, cell_type='CD14 monocyte', type_adata=adata,
        obs=obs, design_terms=['age', 'sex'], required_terms=['age', 'sex'],
        all_samples_geometric_mean_total_counts=100.0,
        valid_all_type_totals=np.ones(len(ages), dtype=bool),
    )


def test_optional_age_bin_fit_skips_rank_deficiency_and_keeps_study_means(tmp_path) -> None:
    prepared = age_bin_profile_inputs([20, 21, 30, 31], ['female', 'female', 'male', 'male'])

    result = _fit_age_bin_profile(prepared, 'GENE1')

    assert result.summary['status'] == 'skipped'
    assert 'not estimable for 4 samples' in result.summary['reason']
    assert result.diagnostic is None and result.comparison is None
    assert result.coefficients.empty and result.predictions.empty and result.samples.empty
    assert result.study_trajectories['age_bin'].tolist() == ['20-30', '30-40']
    assert result.study_support['status'].tolist() == ['aligned']
    for name in ('samples', 'coefficients', 'predictions', 'correlations', 'correlated_residuals'):
        output = tmp_path / f'{name}.csv'
        getattr(result, name).to_csv(output, index=False)
        assert pd.read_csv(output).empty


def test_optional_age_bin_fit_skips_saturated_design() -> None:
    prepared = age_bin_profile_inputs([20, 30, 40], ['female'] * 3)

    result = _fit_age_bin_profile(prepared, 'GENE1')

    assert result.summary['status'] == 'skipped'
    assert 'not estimable for 3 samples' in result.summary['reason']


def test_optional_age_bin_fit_retains_estimable_results() -> None:
    prepared = age_bin_profile_inputs([20, 21, 22, 30, 31, 32], ['female', 'male', 'female'] * 2)

    result = _fit_age_bin_profile(prepared, 'GENE1')

    assert result.summary['status'] == 'complete'
    assert 'reason' not in result.summary
    assert len(result.samples) == 6 and len(result.predictions) == 2
    assert np.isfinite(result.samples['residual']).all()
    assert result.diagnostic['n_samples'] == 6
    np.testing.assert_allclose(result.samples['naive_log_expression'], np.log1p(np.arange(2, 8)))


def test_optional_age_bin_fit_keeps_invalid_counts_fatal() -> None:
    prepared = age_bin_profile_inputs([20, 21], ['female', 'male'], counts=np.full((2, 2), 1.5))

    with pytest.raises(ValueError, match='integer'):
        _fit_age_bin_profile(prepared, 'GENE1')


def test_optional_study_slopes_skip_saturated_design() -> None:
    prepared = SimpleNamespace(
        model_metadata=pd.DataFrame({'study': ['a', 'a', 'b', 'b'], 'age': [20, 30, 20, 30]}),
        varying_terms=['age'],
    )

    result = _fit_study_slope_profile(prepared, np.array([1.0, 2.0, 2.0, 4.0]))

    assert result.summary['status'] == 'skipped'
    assert 'not estimable for 4 samples' in result.summary['reason']
    assert result.coefficients.empty and result.slopes.empty and result.predictions.empty
    assert np.isnan(result.residuals).all()


def test_profile_gene_rejects_study_specific_gene_before_fitting(tmp_path) -> None:
    rows = []
    for study in ("synthetic_study_a", "synthetic_study_b"):
        for index in range(4):
            rows.append({
                "study": study,
                "study_site": f"{study}_site",
                "sample": f"{study}_sample_{index}",
                "age": 40 + index * 10,
                "sex": "female" if index % 2 == 0 else "male",
                "n_cells": 15,
                "total_counts": 100 + index,
                "aifi_l2_majority": "CD14 monocyte",
            })
    obs = pd.DataFrame(rows, index=[f"pb_{i}" for i in range(len(rows))])
    var = pd.DataFrame(
        {
            "available_in_synthetic_study_a": [True],
            "available_in_synthetic_study_b": [False],
        },
        index=["SYNTH_STUDY_A_ONLY"],
    )
    input_path = tmp_path / "pseudobulk.h5ad"
    ad.AnnData(X=np.full((len(obs), 1), 20, dtype=np.int64), obs=obs, var=var).write_h5ad(
        input_path,
    )
    with pytest.raises(ValueError, match="not available in every study"):
        _fit_profile(
            input_path=input_path,
            config_path=PIPELINE_PATH,
            cell_type="CD14 monocyte",
            gene="SYNTH_STUDY_A_ONLY",
            output_dir=tmp_path / "profiles",
            cpus=1,
        )


def test_profile_report_renderer_requires_fit_summary(tmp_path) -> None:
    with pytest.raises(FileNotFoundError, match="fit_summary.json"):
        generate_gene_profile_report(
            profile_dir=tmp_path / "missing",
            template_path=ROOT / "reports" / "gene_profile_report.py",
            project_root=ROOT,
        )


def test_profile_report_renderer_captures_routine_tool_output(tmp_path, monkeypatch) -> None:
    profile_dir = tmp_path / "profile"
    profile_dir.mkdir()
    (profile_dir / "fit_summary.json").write_text(json.dumps({
        "gene": "GENE1", "cell_type": "CD14 monocyte",
    }))
    template = tmp_path / "gene_profile_report.py"
    template.write_text("# report template")
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))
        if "--to" in command and command[command.index("--to") + 1] == "html":
            output_name = next(arg.removeprefix("--output=") for arg in command
                               if arg.startswith("--output="))
            output_dir = Path(next(arg.removeprefix("--output-dir=") for arg in command
                                   if arg.startswith("--output-dir=")))
            (output_dir / output_name).write_text("<html><body>Report</body></html>")
        return subprocess.CompletedProcess(
            command, 0, stdout="routine progress\n", stderr="routine warning\n",
        )

    monkeypatch.setattr("pbmc_pipeline.notebook.subprocess.run", record)
    monkeypatch.setattr("pbmc_pipeline.notebook.set_notebook_title", lambda *_: None)
    monkeypatch.setattr(
        "pbmc_pipeline.notebook.populate_html_toc", lambda *_, **__: None,
    )

    report = generate_gene_profile_report(
        profile_dir=profile_dir,
        template_path=template,
        project_root=tmp_path,
    )

    assert report == profile_dir / "report.html"
    assert len(calls) == 3
    assert all(call[1]["capture_output"] for call in calls)
    render_log = (profile_dir / "report_render.log").read_text()
    assert "routine progress" in render_log
    assert "routine warning" in render_log


def test_profile_report_renderer_prints_and_logs_tool_output_on_failure(
    tmp_path, monkeypatch, capsys,
) -> None:
    profile_dir = tmp_path / "profile"
    profile_dir.mkdir()
    (profile_dir / "fit_summary.json").write_text(json.dumps({"gene": "G", "cell_type": "T"}))
    template = tmp_path / "gene_profile_report.py"
    template.write_text("# report template")

    def fail(command, **kwargs):
        return subprocess.CompletedProcess(
            command, 1, stdout="partial output\n", stderr="render failed\n",
        )

    monkeypatch.setattr("pbmc_pipeline.notebook.subprocess.run", fail)

    with pytest.raises(subprocess.CalledProcessError, match="returned non-zero exit status 1"):
        generate_gene_profile_report(
            profile_dir=profile_dir,
            template_path=template,
            project_root=tmp_path,
        )

    captured = capsys.readouterr()
    assert "partial output" in captured.out
    assert "render failed" in captured.err
    render_log = (profile_dir / "report_render.log").read_text()
    assert "partial output" in render_log
    assert "render failed" in render_log


def test_gene_profile_results_keep_small_p_values_in_scientific_notation() -> None:
    result = format_model_results(pd.DataFrame({
        "term": ["age"],
        "coefficient": [0.125678],
        "p_value": [1.2e-12],
        "pvalue": [2.3e-10],
        "padj": [3.4e-8],
    }))

    assert result.loc[0, "coefficient"] == 0.1257
    assert result.loc[0, "p_value"] == "1.20e-12"
    assert result.loc[0, "pvalue"] == "2.30e-10"
    assert result.loc[0, "padj"] == "3.40e-08"


def test_exploratory_infant_adult_contrast_compares_supported_bins() -> None:
    metadata = pd.DataFrame({
        "age_bin": pd.Categorical(
            ["0-10"] * 4 + ["20-30"] * 4,
            categories=["0-10", "20-30"],
            ordered=True,
        ),
    })
    design = "~ age_bin"
    design_matrix = FormulaicContrasts(metadata, design).design_matrix.to_numpy()
    response = np.array([1.0, 1.2, 0.8, 1.0, 4.0, 4.2, 3.8, 4.0])
    fit = sm.OLS(response, design_matrix).fit()

    result = _exploratory_infant_adult_contrast(
        metadata, design, ["0-10", "20-30"], fit, bin_width=10,
    )

    assert result is not None
    assert result["comparison"] == "20-30 versus 0-10 (approximate ages 25 versus 5)"
    assert result["age_bin_effect_log2_count_plus_1"] == pytest.approx(3.0 / np.log(2))
    assert result["age_bin_p_value"] < 0.001


def test_exploratory_infant_adult_contrast_skips_missing_bins() -> None:
    metadata = pd.DataFrame({
        "age_bin": pd.Categorical(["20-30", "30-40"], categories=["20-30", "30-40"]),
    })
    design = "~ age_bin"
    design_matrix = FormulaicContrasts(metadata, design).design_matrix.to_numpy()
    fit = sm.OLS([1.0, 2.0], design_matrix).fit()

    assert _exploratory_infant_adult_contrast(
        metadata, design, ["20-30", "30-40"], fit, bin_width=10,
    ) is None


def test_prepared_profile_model_mismatch_never_refits(tmp_path, monkeypatch):
    from scripts.profile_gene import _load_or_fit_deseq2_model

    prepared = tmp_path / 'prepared'
    prepared.mkdir()
    (prepared / 'manifest.json').write_text('{"cache_key":"from-a-different-input"}')
    (prepared / 'fitted_model.h5ad').write_bytes(b'not opened on mismatch')
    (prepared / 'age_results.csv').write_text('gene,padj\nGENE1,0.01\n')
    config = tmp_path / 'pipeline.json'
    config.write_text('{}')

    def unexpected_fit(*args, **kwargs):
        pytest.fail('A gene task must not refit a mismatched immutable model')

    monkeypatch.setattr('pydeseq2.dds.DeseqDataSet', unexpected_fit)
    with pytest.raises(ValueError, match='prepared gene-profile model does not match'):
        _load_or_fit_deseq2_model(
            input_path=tmp_path / 'input.h5ad', config_path=config,
            cache_parent=tmp_path / 'cache', prepared_model_dir=prepared,
            cell_type='CD14 monocyte', genes=pd.Index(['GENE1']),
            counts_df=pd.DataFrame({'GENE1': [5, 7]}, index=['sample1', 'sample2']),
            metadata=pd.DataFrame(index=['sample1', 'sample2']),
            design_matrix=pd.DataFrame({'Intercept': [1.0, 1.0], 'age': [20.0, 30.0]}),
            alpha=0.05, cpus=1,
        )
    assert not (tmp_path / 'cache').exists()

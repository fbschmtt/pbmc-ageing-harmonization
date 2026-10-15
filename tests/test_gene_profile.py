from __future__ import annotations

import json
import subprocess
from pathlib import Path

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
from scripts.profile_gene import _exploratory_infant_adult_contrast, _fit_profile

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = ROOT / "config" / "pipeline.json"
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

    monkeypatch.setattr("pbmc_pipeline.gene_profile_report.subprocess.run", record)
    monkeypatch.setattr("pbmc_pipeline.gene_profile_report.set_notebook_title", lambda *_: None)
    monkeypatch.setattr(
        "pbmc_pipeline.gene_profile_report.populate_html_toc", lambda *_, **__: None,
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
    (profile_dir / "fit_summary.json").write_text("{}")
    template = tmp_path / "gene_profile_report.py"
    template.write_text("# report template")

    def fail(command, **kwargs):
        return subprocess.CompletedProcess(
            command, 1, stdout="partial output\n", stderr="render failed\n",
        )

    monkeypatch.setattr("pbmc_pipeline.gene_profile_report.subprocess.run", fail)

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

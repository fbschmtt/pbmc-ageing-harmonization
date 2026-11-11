import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline.config import read_json
from pbmc_pipeline.de_diagnostics import (
    copy_diagnostics,
    prepare_de_diagnostics,
    read_table,
    validate_diagnostic_sources,
)
from pbmc_pipeline.differential_expression_report import generate_differential_expression_report

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def fraction_inputs(tmp_path):
    directory = tmp_path / "models"
    directory.mkdir()
    # A planted ten-percentage-point-per-decade slope in an estimable adult model.
    ages = np.arange(20, 70, 2)
    obs = pd.DataFrame({
        "study": "synthetic", "study_site": "site", "sample": [f"s{i}" for i in range(len(ages))],
        "age": ages, "n_cells": ages - 10, "n_cells_in_sample": 100,
        "n_cells_in_sample_l1_parent": 100, "aifi_l1_parent_for_l2": "Monocyte",
    })
    obs.to_csv(directory / "abundance_input.tsv", sep="\t", index=False)
    (directory / "run_metadata.json").write_text(json.dumps({"cell_type": "CD14 monocyte", "models": []}))
    (directory / "cell_type_result.json").write_text(json.dumps({"cell_type": "CD14 monocyte"}))
    config = tmp_path / "pipeline.json"
    config.write_text((ROOT / "config/pipeline.json").read_text())
    return directory, config


def test_preparation_publishes_expected_fraction_effect_and_reuses_unchanged_analysis(
    fraction_inputs, tmp_path, monkeypatch,
):
    models, config = fraction_inputs
    analysis = tmp_path / "analysis"
    metadata = prepare_de_diagnostics(models, analysis, config)
    coefficients = read_table(analysis / "fraction_model_coefficients.tsv", sep="\t")
    assert coefficients.loc[coefficients.covariate == "age", "estimate"].item() == pytest.approx(0.1)
    assert metadata["residuals"] is None
    assert metadata["sources"]["abundance_input.tsv"]
    copy_diagnostics(analysis, models, read_json(config))

    def unexpected_fit(*args, **kwargs):
        pytest.fail("unchanged diagnostics must be reused without fitting")

    monkeypatch.setattr("pbmc_pipeline.de_diagnostics.fit_parent_fraction_models", unexpected_fit)
    assert prepare_de_diagnostics(models, tmp_path / "reused", config) == metadata
    validate_diagnostic_sources(models, analysis)
    with (models / "abundance_input.tsv").open("a") as table:
        table.write("\n")
    with pytest.raises(ValueError, match="inputs differ"):
        validate_diagnostic_sources(models, analysis)


def test_render_only_reads_prepared_tables_and_detects_corruption(fraction_inputs, tmp_path, monkeypatch):
    models, config = fraction_inputs
    analysis = tmp_path / "analysis"
    prepare_de_diagnostics(models, analysis, config)
    calls = []
    monkeypatch.setattr("pbmc_pipeline.differential_expression_report.render_notebook",
                        lambda **kwargs: calls.append(kwargs) or kwargs["output_dir"] / "report.html")

    def unexpected_fit(*args, **kwargs):
        pytest.fail("rendering must not run diagnostic calculations")

    monkeypatch.setattr("pbmc_pipeline.de_diagnostics.prepare_de_diagnostics", unexpected_fit)
    output = tmp_path / "report"
    generate_differential_expression_report(
        result_dir=models, analysis_dir=analysis, output_dir=output, config_path=config,
        template_path=ROOT / "reports/differential_expression_report.py", project_root=ROOT,
    )
    assert calls[0]["environment"]["DE_REPORT_ANALYSIS_DIR"] == str(output)
    assert (output / "fraction_model_coefficients.tsv").is_file()
    (analysis / "fraction_model_coefficients.tsv").write_text("changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        copy_diagnostics(analysis, output, read_json(config))


def test_prepared_diagnostics_reject_changed_configuration(fraction_inputs, tmp_path):
    models, config = fraction_inputs
    analysis = tmp_path / "analysis"
    prepare_de_diagnostics(models, analysis, config)
    pipeline = read_json(config)
    pipeline["differential_expression"]["alpha"] = 0.01
    with pytest.raises(ValueError, match="configuration changed"):
        copy_diagnostics(analysis, tmp_path / "report", pipeline)


def test_changed_enrichment_inputs_invalidate_prepared_diagnostics(fraction_inputs, tmp_path, monkeypatch):
    models, config = fraction_inputs
    resources = tmp_path / "go"
    resources.mkdir()
    library = resources / "go_biological_process_2023.gmt"
    library.write_text("first library")
    monkeypatch.setenv("PBMC_GO_ENRICHMENT_DATA_DIR", str(resources))
    prepared = prepare_de_diagnostics(models, tmp_path / "analysis", config)
    copy_diagnostics(tmp_path / "analysis", models, read_json(config))
    library.write_text("changed library")
    refreshed = prepare_de_diagnostics(models, tmp_path / "refreshed", config)
    assert refreshed["enrichment_resources"] != prepared["enrichment_resources"]


def test_preparation_repairs_a_corrupt_published_diagnostic_table(fraction_inputs, tmp_path):
    models, config = fraction_inputs
    prepared = prepare_de_diagnostics(models, tmp_path / "analysis", config)
    copy_diagnostics(tmp_path / "analysis", models, read_json(config))
    (models / "fraction_model_coefficients.tsv").write_text("corrupt")
    repaired = prepare_de_diagnostics(models, tmp_path / "repaired", config)
    assert repaired["artifacts"] == prepared["artifacts"]
    coefficients = read_table(tmp_path / "repaired/fraction_model_coefficients.tsv", sep="\t")
    assert coefficients.loc[coefficients.covariate == "age", "estimate"].item() == pytest.approx(0.1)

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest
import seaborn as sns

from pbmc_pipeline.cell_type_report import generate_cell_type_report


def test_cell_type_report_runner_executes_notebook_then_exports_html(tmp_path, monkeypatch):
    input_path = tmp_path / "type.h5ad"
    analysis_path = tmp_path / "analysis.h5ad"
    template = tmp_path / "cell_type_report.py"
    config = tmp_path / "pipeline.json"
    for path in (input_path, analysis_path, template, config):
        path.write_text("{}")
    config.write_text((Path(__file__).resolve().parents[1] / "config/pipeline.json").read_text())
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))
        if "--execute" in command:
            output_dir = Path(next(arg.removeprefix("--output-dir=") for arg in command if arg.startswith("--output-dir=")))
            (output_dir / "executed.ipynb").write_text('{"cells": [], "metadata": {}, "nbformat": 4, "nbformat_minor": 5}')
        if "--to" in command and command[command.index("--to") + 1] == "html":
            output_dir = Path(next(arg.removeprefix("--output-dir=") for arg in command if arg.startswith("--output-dir=")))
            output_name = next(arg.removeprefix("--output=") for arg in command if arg.startswith("--output="))
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / output_name).write_text(
                '<nav class="report-toc"><ul><!-- REPORT_TOC_PLACEHOLDER -->'
                '</ul></nav><h2 id="main-section">Main section</h2>'
                '<h3 id="subsection">Subsection</h3>'
                '<h4 id="sub-subsection">Sub-subsection</h4>'
            )

        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("pbmc_pipeline.notebook.subprocess.run", record)
    html = generate_cell_type_report(
        input_path=input_path,
        analysis_path=analysis_path,
        output_dir=tmp_path / "report",
        template_path=template,
        config_path=config,
    )

    assert html == tmp_path / "report" / "report.html"
    assert len(calls) == 3
    assert "jupytext" in calls[0][0]
    assert "--execute" in calls[1][0]
    assert calls[1][1]["env"]["CELL_TYPE_H5AD"] == str(input_path.resolve())
    assert calls[1][1]["env"]["CELL_TYPE_ANALYSIS_H5AD"] == str(analysis_path.resolve())
    assert "--HTMLExporter.exclude_input=True" in calls[2][0]
    assert "--HTMLExporter.exclude_output_prompt=True" in calls[2][0]
    rendered_html = html.read_text()
    assert "REPORT_TOC_PLACEHOLDER" not in rendered_html
    assert '<a href="#main-section">Main section</a>' in rendered_html
    assert '<a href="#subsection">Subsection</a>' in rendered_html
    assert (
        '<a href="#subsection">Subsection</a><ul>'
        '<li><a href="#sub-subsection">Sub-subsection</a></li></ul>'
    ) in rendered_html


def test_cluster_fraction_plot_draws_pooled_and_only_estimable_study_trends(monkeypatch):
    from pbmc_pipeline.cell_type_figures import plot_cluster_fractions

    calls = []
    monkeypatch.setattr(sns, "regplot", lambda **kwargs: calls.append(kwargs))
    figures = []
    data = pd.DataFrame({
        "cluster": ["0"] * 5, "study": ["a", "a", "b", "b", "c"],
        "age": [20, 40, 30, 50, 60], "fraction_within_cell_type": [0, 1, 0.5, 0, 0.2],
    })
    try:
        result = plot_cluster_fractions(data, palette={"a": "red", "b": "blue", "c": "green"}, display_figure=lambda figure, *_: figures.append(figure))
        assert [call["color"] for call in calls] == ["red", "blue", "#1f1f1f"] * 2
        assert [len(call["data"]) for call in calls] == [2, 2, 5] * 2
        fits = result[0]["fits"]
        assert fits[0]["study"] is None
        assert fits[0]["meanFraction"] == pytest.approx(0.34)
        assert fits[1]["slopePpPerDecade"] == pytest.approx(50)
        assert fits[3]["slopePpPerDecade"] is None
    finally:
        for figure in figures:
            plt.close(figure)

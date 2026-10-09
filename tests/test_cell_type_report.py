from pathlib import Path

from pbmc_pipeline.cell_type_report import generate_cell_type_report


def test_cell_type_report_runner_executes_notebook_then_exports_html(tmp_path, monkeypatch):
    input_path = tmp_path / "type.h5ad"
    analysis_path = tmp_path / "analysis.h5ad"
    template = tmp_path / "cell_type_report.py"
    config = tmp_path / "pipeline.json"
    for path in (input_path, analysis_path, template, config):
        path.write_text("{}")
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))
        if "--to" in command and command[command.index("--to") + 1] == "html":
            output_dir = Path(command[command.index("--output-dir") + 1])
            output_name = command[command.index("--output") + 1]
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / output_name).write_text(
                '<nav class="report-toc"><ul><!-- REPORT_TOC_PLACEHOLDER -->'
                '</ul></nav><h2 id="main-section">Main section</h2>'
                '<h3 id="subsection">Subsection</h3>'
                '<h4 id="sub-subsection">Sub-subsection</h4>'
            )

    monkeypatch.setattr("pbmc_pipeline.cell_type_report.subprocess.run", record)
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

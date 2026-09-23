from pbmc_pipeline.cell_type_report import generate_cell_type_report


def test_cell_type_report_runner_executes_notebook_then_exports_html(tmp_path, monkeypatch):
    input_path = tmp_path / "type.h5ad"
    template = tmp_path / "cell_type_report.py"
    config = tmp_path / "pipeline.json"
    for path in (input_path, template, config):
        path.write_text("{}")
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))

    monkeypatch.setattr("pbmc_pipeline.cell_type_report.subprocess.run", record)
    html = generate_cell_type_report(
        input_path=input_path,
        output_dir=tmp_path / "report",
        template_path=template,
        config_path=config,
    )

    assert html == tmp_path / "report" / "report.html"
    assert len(calls) == 3
    assert "jupytext" in calls[0][0]
    assert "--execute" in calls[1][0]
    assert calls[1][1]["env"]["CELL_TYPE_H5AD"] == str(input_path.resolve())

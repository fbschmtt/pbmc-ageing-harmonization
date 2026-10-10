from pbmc_pipeline.qc import generate_qc_report


def test_qc_runner_materializes_executes_and_exports_report(tmp_path, monkeypatch):
    root = tmp_path
    template = root / "reports" / "qc_report.py"
    template.parent.mkdir()
    template.write_text("{}")
    input_path = root / "study.h5ad"
    input_path.touch()
    run_report = root / "study.json"
    run_report.write_text("{}")
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))
        if "--execute" in command:
            output_dir = root / "output" / "harmonization_qc" / "study"
            (output_dir / "executed.ipynb").write_text(
                '{"cells": [], "metadata": {}, "nbformat": 4, "nbformat_minor": 5}'
            )
        if "html" in command:
            output_dir = root / "output" / "harmonization_qc" / "study"
            (output_dir / "report.html").write_text(
                '<nav><ul><!-- REPORT_TOC_PLACEHOLDER --></ul></nav>'
                '<h2 id="run-provenance">Run provenance and dimensions</h2>'
            )

    monkeypatch.setattr("pbmc_pipeline.qc.subprocess.run", record)
    html = generate_qc_report(
        root=root,
        input_path=input_path,
        run_report_path=run_report,
        output_dir=root / "output" / "harmonization_qc" / "study",
        study="study",
    )

    assert html == root / "output" / "harmonization_qc" / "study" / "report.html"
    assert len(calls) == 3
    assert "jupytext" in calls[0][0]
    assert "--execute" in calls[1][0]
    assert "html" in calls[2][0]
    assert calls[0][1]["env"]["QC_STUDY"] == "study"
    assert "--HTMLExporter.exclude_input=True" in calls[2][0]
    assert "--HTMLExporter.exclude_input_prompt=True" in calls[2][0]
    assert "--HTMLExporter.exclude_output_prompt=True" in calls[2][0]
    assert '<a href="#run-provenance">Run provenance and dimensions</a>' in html.read_text()

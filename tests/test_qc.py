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

    monkeypatch.setattr("pbmc_pipeline.qc.subprocess.run", record)
    html = generate_qc_report(
        root=root,
        input_path=input_path,
        run_report_path=run_report,
        output_dir=root / "output" / "qc" / "study",
        study="study",
    )

    assert html == root / "output" / "qc" / "study" / "report.html"
    assert len(calls) == 3
    assert "jupytext" in calls[0][0]
    assert "--execute" in calls[1][0]
    assert "html" in calls[2][0]
    assert calls[0][1]["env"]["QC_STUDY"] == "study"

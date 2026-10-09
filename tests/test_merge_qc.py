import json

from pbmc_pipeline.merge_qc import generate_merge_qc_report


def test_merge_qc_runner_uses_a_writable_temporary_matplotlib_cache(tmp_path, monkeypatch):
    template = tmp_path / "merge_qc_report.py"
    template.write_text("{}")
    pseudobulk = tmp_path / "merged.h5ad"
    pseudobulk.touch()
    run_report = tmp_path / "merged.json"
    run_report.write_text(json.dumps({"kind": "pseudobulk_merge"}))
    studies_config = tmp_path / "studies.json"
    studies_config.write_text("{}")
    calls = []

    def record(command, **kwargs):
        calls.append((command, kwargs))
        if command[command.index("--to") + 1] == "html":
            report = tmp_path / "output" / "qc" / "merged" / "report.html"
            report.write_text(
                '<nav><!-- REPORT_TOC_PLACEHOLDER --></nav>'
                '<h2 id="merge-summary">Merge Summary</h2>'
            )

    monkeypatch.setattr("pbmc_pipeline.merge_qc.subprocess.run", record)
    html = generate_merge_qc_report(
        root=tmp_path,
        inputs=[pseudobulk],
        reports=[run_report],
        output_dir=tmp_path / "output" / "qc" / "merged",
        template_path=template,
        studies_config_path=studies_config,
    )

    assert html == tmp_path / "output" / "qc" / "merged" / "report.html"
    assert len(calls) == 3
    matplotlib_cache = calls[0][1]["env"]["MPLCONFIGDIR"]
    assert matplotlib_cache.startswith("/tmp/pbmc-merge-qc-notebook-")
    assert matplotlib_cache.endswith("/matplotlib")
    assert '<a href="#merge-summary">Merge Summary</a>' in html.read_text()

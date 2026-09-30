import json
import sys

from pbmc_pipeline.manifest import main


def test_manifest_records_requested_selected_and_skipped_studies(tmp_path, monkeypatch):
    config = tmp_path / "pipeline.json"
    artifact = tmp_path / "merged.h5ad"
    report = tmp_path / "merged.json"
    output = tmp_path / "run_manifest.json"
    config.write_text("{}")
    artifact.write_bytes(b"artifact")
    report.write_text("{}")
    monkeypatch.setattr(sys, "argv", [
        "pbmc-manifest", "--config", str(config), "--artifact", str(artifact), "--report", str(report),
        "--requested-studies", "aifi,terekhova23,wang25", "--studies", "aifi",
        "--skipped-studies", "terekhova23,wang25", "--pipeline-revision", "test",
        "--python-image", "python:test", "--r-image", "r:test", "--output", str(output),
    ])

    main()

    manifest = json.loads(output.read_text())
    assert manifest["requested_studies"] == ["aifi", "terekhova23", "wang25"]
    assert manifest["selected_studies"] == ["aifi"]
    assert manifest["skipped_studies"] == ["terekhova23", "wang25"]

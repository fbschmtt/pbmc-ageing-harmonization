import json

import pytest

from pbmc_pipeline.cell_type_manifest import write_cell_type_manifest


def test_cell_type_manifest_validates_and_records_every_report(tmp_path):
    input_path = tmp_path / "merged.h5ad"
    config_path = tmp_path / "pipeline.json"
    model_path = tmp_path / "l2.pkl"
    for path in (input_path, model_path):
        path.write_bytes(b"content")
    config_path.write_text(json.dumps({"cell_type_analysis": {"analysis_version": 1}}))
    report_dir = tmp_path / "naive-cd4-t-cell"
    report_dir.mkdir()
    (report_dir / "analysis.h5ad").write_bytes(b"analysis")
    (report_dir / "report.html").write_text("<html></html>")
    (report_dir / "executed.ipynb").write_text("{}")
    (report_dir / "report.json").write_text(json.dumps({
        "status": "complete", "cell_type": "Naive CD4 T cell", "aifi_l1_parent": "T cell",
        "n_cells": 50, "n_clusters": 2,
    }))

    document = write_cell_type_manifest(
        input_path=input_path, config_path=config_path, model_paths=[model_path],
        report_dirs=[report_dir], output_path=tmp_path / "manifest.json",
    )

    assert document["analysis_specification"]["analysis_version"] == 1
    record = document["cell_types"]
    assert len(record) == 1
    assert record[0].items() >= {
        "slug": "naive-cd4-t-cell", "cell_type": "Naive CD4 T cell",
        "aifi_l1_parent": "T cell", "status": "complete", "n_cells": 50,
        "n_clusters": 2,
    }.items()
    assert len(record[0]["analysis_h5ad_sha256"]) == 64
    assert len(record[0]["report_html_sha256"]) == 64


def test_cell_type_manifest_rejects_incomplete_report_directory(tmp_path):
    input_path = tmp_path / "merged.h5ad"
    config_path = tmp_path / "pipeline.json"
    model_path = tmp_path / "l2.pkl"
    for path in (input_path, model_path):
        path.write_bytes(b"content")
    config_path.write_text(json.dumps({"cell_type_analysis": {"analysis_version": 1}}))
    report_dir = tmp_path / "incomplete"
    report_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="required report outputs"):
        write_cell_type_manifest(
            input_path=input_path, config_path=config_path, model_paths=[model_path],
            report_dirs=[report_dir], output_path=tmp_path / "manifest.json",
        )

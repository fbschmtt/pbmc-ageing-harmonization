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
    (report_dir / "age_trajectory_clusters.csv").write_text(
        "gene,trajectory_cluster\nGENE1,1\n"
    )
    (report_dir / "age_trajectory_cluster_means.csv").write_text(
        "trajectory_cluster,n_trajectories\n1,1\n"
    )
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
    assert set(record[0]["age_trajectory_outputs"]) == {
        "age_trajectory_clusters.csv", "age_trajectory_cluster_means.csv",
    }
    assert all(
        len(item["sha256"]) == 64
        for item in record[0]["age_trajectory_outputs"].values()
    )


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


@pytest.fixture
def manifest_run(tmp_path):
    input_path = tmp_path / "merged.h5ad"
    input_path.write_bytes(b"input")
    config = tmp_path / "pipeline.json"
    config.write_text(json.dumps({"cell_type_analysis": {"analysis_version": 1}}))
    split = tmp_path / "cell_types.json"
    split.write_text(json.dumps({"cell_types": [
        {"slug": "cd14-monocyte", "cell_type": "CD14 monocyte"},
        {"slug": "naive-cd4-t-cell", "cell_type": "Naive CD4 T cell"},
    ]}))
    return {"input_path": input_path, "config_path": config, "model_paths": [], "report_dirs": [],
            "output_path": tmp_path / "manifest.json", "split_manifest_path": split}


@pytest.mark.parametrize("selected", [["typo"], ["cd14-monocyte", "typo"], [""]])
def test_manifest_rejects_unknown_selections_without_publishing(manifest_run, selected):
    with pytest.raises(ValueError, match="Requested cell types are absent.*Available"):
        write_cell_type_manifest(**manifest_run, selected_slugs=selected)
    assert not manifest_run["output_path"].exists()


def test_manifest_keeps_real_analysis_and_report_failures_visible(manifest_run):
    document = write_cell_type_manifest(**manifest_run, analysed_slugs=["cd14-monocyte"])
    assert document["status"] == "partial"
    assert document["failed_cell_types"] == [
        {"slug": "cd14-monocyte", "cell_type": "CD14 monocyte", "failed_stage": "report"},
        {"slug": "naive-cd4-t-cell", "cell_type": "Naive CD4 T cell", "failed_stage": "analysis"},
    ]
    selected = write_cell_type_manifest(**manifest_run, selected_slugs=["cd14-monocyte"])
    assert selected["status"] == "partial"
    assert [record["slug"] for record in selected["failed_cell_types"]] == ["cd14-monocyte"]


def test_manifest_cannot_publish_an_empty_run_as_complete(manifest_run):
    manifest_run["split_manifest_path"].write_text('{"cell_types": []}')
    with pytest.raises(ValueError, match="non-empty"):
        write_cell_type_manifest(**manifest_run)

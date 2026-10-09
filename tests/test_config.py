import copy
import json
from pathlib import Path

import pytest

from pbmc_pipeline.config import (
    AgeTrajectorySettings,
    ConfigurationError,
    load_configuration,
    read_json,
    validate_configuration,
    validate_document,
)

ROOT = Path(__file__).resolve().parents[1]


def test_configuration_is_complete():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    assert set(studies["studies"]) == {
        "aida25", "aifi", "fachrul26", "nehar_belaid26", "onek1k", "perez22",
        "terekhova23", "wang25"
    }
    assert pipeline["schema_version"] == schema["schema_version"] == 1


def test_age_trajectory_report_and_clustering_choices_are_explicit_and_validated():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    age_trajectory = pipeline["differential_expression"]["age_trajectory"]
    settings = AgeTrajectorySettings.from_mapping(age_trajectory)
    assert settings.to_mapping() == age_trajectory
    assert age_trajectory["minimum_bins"] == 7
    assert age_trajectory["cluster_fdr_threshold"] == 0.001
    assert age_trajectory["de_fdr_threshold"] == 0.05
    assert age_trajectory["max_clusters"] == 7
    assert age_trajectory["linkage_method"] == "ward"
    assert age_trajectory["distance_metric"] == "euclidean"
    assert age_trajectory["minimum_umap_trajectories"] == 4
    assert age_trajectory["residual_umap_neighbors"] == 50
    assert age_trajectory["report_top_n_genes"] == 30
    assert age_trajectory["cross_report_top_n_genes"] == 100
    assert pipeline["differential_expression"]["expression_atlas"]["minimum_study_cell_type_total_counts"] == 1_000_000

    input_sources = read_json(ROOT / pipeline["input_sources"])
    invalid = copy.deepcopy(pipeline)
    invalid["differential_expression"]["age_trajectory"]["umap_min_dist"] = 1.1
    with pytest.raises(ConfigurationError, match="umap_min_dist"):
        validate_configuration(invalid, studies, schema, input_sources)

    incomplete = age_trajectory.copy()
    incomplete.pop("umap_min_dist")
    with pytest.raises(ConfigurationError, match="missing=.*umap_min_dist"):
        AgeTrajectorySettings.from_mapping(incomplete)


def test_cell_type_parent_mapping_matches_the_shipped_aifi_l2_model(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLTYPIST_FOLDER", str(tmp_path / "celltypist"))
    from celltypist import models

    pipeline, _, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    model = models.Model.load(str(ROOT / pipeline["models"]["aifi_l2"]))

    assert set(pipeline["cell_type_analysis"]["l2_parent_l1"]) == set(model.cell_types)


def test_all_json_files_are_valid():
    for path in (ROOT / "config").glob("*.json"):
        with path.open(encoding="utf-8") as handle:
            assert isinstance(json.load(handle), dict)


def test_test_inputs_are_unique():
    pipeline, document, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    names = [study["test_input"] for study in document["studies"].values()]
    assert len(names) == len(set(names))
    assert pipeline["test_input_root"] != "input_data"
    assert document["studies"]["wang25"]["conversion"]["test_source"].startswith(
        f'{pipeline["test_input_root"]}/'
    )


def test_new_cellxgene_studies_use_gene_symbol_columns():
    _, document, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    for study_id in ("fachrul26", "perez22"):
        assert document["studies"][study_id]["features"] == {
            "source_column": "feature_name", "duplicate_policy": "sum"
        }


def test_perez_healthy_only_selection_is_explicitly_permitted_and_recorded():
    _, document, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    perez = document["studies"]["perez22"]

    assert perez["preparation"]["allow_cell_subset"] is True
    assert "healthy cells only" in perez["provenance"]["selection"]


def test_configuration_rejects_invalid_feature_policy():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    input_sources = read_json(ROOT / pipeline["input_sources"])
    broken = copy.deepcopy(studies)
    broken["studies"]["fachrul26"]["features"]["duplicate_policy"] = "keep"

    with pytest.raises(ConfigurationError, match="duplicate_policy"):
        validate_configuration(pipeline, broken, schema, input_sources)


def test_configuration_rejects_unmanifested_expression_input():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    input_sources = read_json(ROOT / pipeline["input_sources"])
    broken = copy.deepcopy(input_sources)
    broken["artifacts"] = [
        artifact for artifact in broken["artifacts"] if "fachrul26" not in artifact["studies"]
    ]

    with pytest.raises(ConfigurationError, match="fachrul26"):
        validate_configuration(pipeline, studies, schema, broken)


def test_study_registry_schema_rejects_an_undeclared_field():
    _, studies, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    broken = copy.deepcopy(studies)
    broken["studies"]["aifi"]["unexpected"] = True

    with pytest.raises(ConfigurationError, match="Additional properties"):
        validate_document(
            broken,
            read_json(ROOT / "config/studies.schema.json"),
            document_name="config/studies.json",
        )

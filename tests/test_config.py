import copy
import json
from pathlib import Path

import pytest

from pbmc_pipeline.config import (
    AgeTrajectorySettings,
    ConfigurationError,
    ResidualWGCNASettings,
    load_configuration,
    read_json,
    validate_configuration,
    validate_document,
)

ROOT = Path(__file__).resolve().parents[1]


def test_configuration_is_complete():
    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    assert set(studies["studies"]) == {
        "aida25", "aifi", "cima", "fachrul26", "nehar_belaid26", "onek1k", "perez22",
        "terekhova23", "terekhova26", "wang25"
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
    assert pipeline["differential_expression"][
        "study_comparison_minimum_samples_per_bin"
    ] == 5
    assert pipeline["differential_expression"]["expression_atlas"]["minimum_study_cell_type_total_counts"] == 1_000_000
    wgcna = pipeline["cell_type_analysis"]["residual_wgcna"]
    assert ResidualWGCNASettings.from_mapping(wgcna).to_mapping() == wgcna
    assert wgcna["maximum_genes"] == 2000
    assert wgcna["minimum_module_size"] == 20

    input_sources = read_json(ROOT / pipeline["input_sources"])
    invalid = copy.deepcopy(pipeline)
    invalid["differential_expression"]["age_trajectory"]["umap_min_dist"] = 1.1
    with pytest.raises(ConfigurationError, match="umap_min_dist"):
        validate_configuration(invalid, studies, schema, input_sources)

    invalid_study_comparison = copy.deepcopy(pipeline)
    invalid_study_comparison["differential_expression"][
        "study_comparison_minimum_samples_per_bin"
    ] = 1
    with pytest.raises(ConfigurationError, match="study_comparison_minimum_samples_per_bin"):
        validate_configuration(invalid_study_comparison, studies, schema, input_sources)

    incomplete = age_trajectory.copy()
    incomplete.pop("umap_min_dist")
    with pytest.raises(ConfigurationError, match="missing=.*umap_min_dist"):
        AgeTrajectorySettings.from_mapping(incomplete)


@pytest.mark.external_models
def test_cell_type_parent_mapping_matches_the_shipped_aifi_l2_model(tmp_path, monkeypatch):
    monkeypatch.setenv("CELLTYPIST_FOLDER", str(tmp_path / "celltypist"))
    from celltypist import models

    pipeline, _, _ = load_configuration(ROOT, Path("config/pipeline.json"))
    model_path = ROOT / pipeline["models"]["aifi_l2"]
    assert model_path.is_file(), "Acquire the AIFI models described in aifi_models/README before make test-models"
    model = models.Model.load(str(model_path))

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


def test_staged_analysis_config_needs_no_registries_and_preserves_digest(tmp_path):
    from pbmc_pipeline.config import config_digest, load_analysis_configuration

    pipeline = read_json(ROOT / "config/pipeline.json")
    pipeline["studies_config"] = "not-staged/studies.json"
    pipeline["obs_schema"] = "not-staged/schema.json"
    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(pipeline))
    loaded = load_analysis_configuration(staged)
    assert loaded == pipeline
    assert config_digest(loaded) == config_digest(pipeline)
    assert list(tmp_path.iterdir()) == [staged]


@pytest.mark.parametrize(("section", "key", "value", "message"), [
    ("cell_type_analysis", "min_cells", 1, "min_cells"),
    ("cell_type_analysis", "clustering", [], "clustering must be an object"),
    ("cell_type_analysis", "residual_wgcna", {}, "missing="),
    ("differential_expression", "alpha", 2, "alpha"),
    ("differential_expression", "abundance", None, "abundance must be an object"),
    ("differential_expression", "sample_inclusion", [], "sample_inclusion must be an object"),
    ("differential_expression", "age_trajectory", {}, "missing="),
])
def test_core_and_downstream_configuration_share_analysis_rules(tmp_path, section, key, value, message):
    from pbmc_pipeline.config import load_analysis_configuration

    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    pipeline[section][key] = value
    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(pipeline))
    with pytest.raises(ConfigurationError, match=message):
        validate_configuration(pipeline, studies, schema)
    with pytest.raises(ConfigurationError, match=message):
        load_analysis_configuration(staged)


@pytest.mark.parametrize("document", [[], {}, {"schema_version": 2}])
def test_analysis_loader_rejects_missing_or_unsupported_schema(tmp_path, document):
    from pbmc_pipeline.config import load_analysis_configuration

    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(document))
    with pytest.raises(ConfigurationError, match="schema_version"):
        load_analysis_configuration(staged)


def test_analysis_loader_reports_malformed_json(tmp_path):
    from pbmc_pipeline.config import load_analysis_configuration

    staged = tmp_path / "pipeline.json"
    staged.write_text('{"schema_version":')
    with pytest.raises(ConfigurationError, match="pipeline.json: invalid JSON"):
        load_analysis_configuration(staged)


@pytest.mark.parametrize("value", [0, True, "ten"])
def test_clock_and_de_share_sample_inclusion_validation(tmp_path, value):
    from pbmc_pipeline.config import load_analysis_configuration

    pipeline = read_json(ROOT / "config/pipeline.json")
    pipeline["differential_expression"]["sample_inclusion"]["minimum_cells_per_pseudobulk_inclusive"] = value
    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(pipeline))
    for sections in (("age_clock",), ("differential_expression",)):
        with pytest.raises(ConfigurationError, match="minimum_cells_per_pseudobulk"):
            load_analysis_configuration(staged, sections=sections)


def test_clock_validation_only_requires_consumed_settings(tmp_path):
    from pbmc_pipeline.config import load_analysis_configuration

    pipeline = read_json(ROOT / "config/pipeline.json")
    spec = pipeline["differential_expression"]
    staged_document = {"schema_version": 1, "differential_expression": {
        "split_by": spec["split_by"], "sample_inclusion": spec["sample_inclusion"],
    }}
    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(staged_document))
    assert load_analysis_configuration(staged, sections=("age_clock",)) == staged_document


@pytest.mark.parametrize("command", ["cell_type", "differential_expression", "age_clock", "de_diagnostics"])
def test_invalid_analysis_configuration_fails_before_reading_data(tmp_path, monkeypatch, command):
    from pbmc_pipeline import age_clock, cell_type, de_diagnostics, differential_expression

    pipeline = read_json(ROOT / "config/pipeline.json")
    pipeline["cell_type_analysis"]["min_cells"] = 1
    pipeline["differential_expression"]["sample_inclusion"]["minimum_cells_per_pseudobulk_inclusive"] = 0
    config = tmp_path / "pipeline.json"
    config.write_text(json.dumps(pipeline))
    missing = tmp_path / "missing.h5ad"
    output = tmp_path / "output"
    with pytest.raises(ConfigurationError):
        if command == "age_clock":
            age_clock.fit_age_clocks(missing, config, ["held-out"], cell_type_filter="cd14-monocyte")
        elif command == "de_diagnostics":
            de_diagnostics.prepare_de_diagnostics(missing, output, config)
        else:
            arguments = [command]
            if command == "cell_type":
                arguments += ["analyse"]
            arguments += ["--input", str(missing), "--config", str(config), "--output-dir", str(output)]
            monkeypatch.setattr("sys.argv", arguments)
            {"cell_type": cell_type.main, "differential_expression": differential_expression.main}[command]()
    assert not output.exists()


@pytest.mark.parametrize(("key", "value"), [
    ("target_sum", 0), ("target_sum", float("inf")), ("hvg", []),
    ("exclude_vdj_regex", "["), ("random_seed", True),
    ("neighbors", {"n_neighbors": 50, "n_pcs": 0}), ("output_compression", "typo"),
])
def test_cell_type_processing_rules_apply_to_core_and_staged_configs(tmp_path, key, value):
    from pbmc_pipeline.config import load_analysis_configuration

    pipeline, studies, schema = load_configuration(ROOT, Path("config/pipeline.json"))
    pipeline["processing"][key] = value
    staged = tmp_path / "pipeline.json"
    staged.write_text(json.dumps(pipeline))
    with pytest.raises(ConfigurationError, match=key):
        validate_configuration(pipeline, studies, schema)
    with pytest.raises(ConfigurationError, match=key):
        load_analysis_configuration(staged, sections=("cell_type_analysis", "processing"))

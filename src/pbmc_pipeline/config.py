from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


class ConfigurationError(ValueError):
    """Raised when pipeline configuration is incomplete or inconsistent."""


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_configuration(root: Path, pipeline_path: Path) -> tuple[dict, dict, dict]:
    pipeline = read_json(root / pipeline_path)
    studies = read_json(root / pipeline["studies_config"])
    schema = read_json(root / pipeline["obs_schema"])
    input_sources_path = pipeline.get("input_sources")
    if not isinstance(input_sources_path, str) or not input_sources_path:
        raise ConfigurationError("pipeline.input_sources must be a non-empty path")
    input_sources = read_json(root / input_sources_path)
    validate_configuration(pipeline, studies, schema, input_sources)
    return pipeline, studies, schema


def validate_configuration(
    pipeline: dict,
    studies: dict,
    schema: dict,
    input_sources: dict | None = None,
) -> None:
    if pipeline.get("schema_version") != 1 or studies.get("schema_version") != 1:
        raise ConfigurationError("Only configuration schema_version 1 is supported")
    if not studies.get("studies"):
        raise ConfigurationError("No studies are configured")
    test_inputs: list[str] = []
    merge = pipeline.get("merge", {})
    groupby = merge.get("pseudobulk", {}).get("groupby", [])
    if not groupby or not {"sample", "aifi_l2_majority"}.issubset(groupby):
        raise ConfigurationError(
            "merge.pseudobulk.groupby must include sample and aifi_l2_majority"
        )
    single_cell = merge.get("single_cell", {})
    if single_cell.get("gene_join") not in {"inner"}:
        raise ConfigurationError("merge.single_cell.gene_join currently supports only 'inner'")
    if single_cell.get("annotation_level") != "l2":
        raise ConfigurationError("merge.single_cell.annotation_level must be 'l2'")
    integration = single_cell.get("integration", {})
    if integration.get("method") != "harmony":
        raise ConfigurationError("merge.single_cell.integration.method must be 'harmony'")
    if not isinstance(integration.get("batch_key"), str) or not integration["batch_key"]:
        raise ConfigurationError("merge.single_cell.integration.batch_key must be a non-empty string")
    if integration.get("adjusted_basis") != "X_pca_harmony":
        raise ConfigurationError(
            "merge.single_cell.integration.adjusted_basis must be 'X_pca_harmony'"
        )

    for study_id, study in studies["studies"].items():
        missing_keys = {"input", "test_input", "counts_source", "preparation", "annotation"} - study.keys()
        if missing_keys:
            raise ConfigurationError(f"{study_id}: missing keys {sorted(missing_keys)}")
        if not isinstance(study["input"], str) or not study["input"]:
            raise ConfigurationError(f"{study_id}: input must be a non-empty path")
        if not isinstance(study["test_input"], str) or not study["test_input"]:
            raise ConfigurationError(f"{study_id}: test_input must be a non-empty filename")
        test_inputs.append(study["test_input"])
        counts_source = study["counts_source"]
        if not isinstance(counts_source, str) or not counts_source:
            raise ConfigurationError(f"{study_id}: counts_source must be a non-empty string")
        if counts_source not in {"X", "raw"} and not counts_source.startswith("layer:"):
            raise ConfigurationError(f"{study_id}: unsupported counts_source {counts_source!r}")
        if counts_source == "layer:":
            raise ConfigurationError(f"{study_id}: layer counts_source must name a layer")

        features = study.get("features", {})
        if not isinstance(features, dict):
            raise ConfigurationError(f"{study_id}: features must be an object")
        source_column = features.get("source_column")
        if source_column is not None and (not isinstance(source_column, str) or not source_column):
            raise ConfigurationError(f"{study_id}: features.source_column must be a non-empty string")
        duplicate_policy = features.get("duplicate_policy", "error")
        if duplicate_policy not in {"error", "sum"}:
            raise ConfigurationError(
                f"{study_id}: features.duplicate_policy must be 'error' or 'sum'"
            )
        split = features.get("split")
        if split is not None and (
            not isinstance(split, dict)
            or not isinstance(split.get("separator"), str)
            or not isinstance(split.get("index"), int)
        ):
            raise ConfigurationError(
                f"{study_id}: features.split requires string separator and integer index"
            )

        if not isinstance(study["preparation"].get("adapter"), str):
            raise ConfigurationError(f"{study_id}: preparation.adapter must be a string")
        adapter = study["preparation"]["adapter"]
        if importlib.util.find_spec(f"pbmc_pipeline.studies.{adapter}") is None:
            raise ConfigurationError(f"{study_id}: preparation adapter module does not exist: {adapter}")
        if not isinstance(study["preparation"].get("dependencies", []), list):
            raise ConfigurationError(f"{study_id}: preparation.dependencies must be a list")
        if not all(isinstance(path, str) and path for path in study["preparation"].get("dependencies", [])):
            raise ConfigurationError(f"{study_id}: preparation.dependencies must contain paths")
        if not isinstance(study["preparation"].get("allow_cell_subset", False), bool):
            raise ConfigurationError(f"{study_id}: preparation.allow_cell_subset must be a boolean")

        annotation = study["annotation"]
        if not isinstance(annotation, dict) or annotation.get("method") not in {"existing", "celltypist"}:
            raise ConfigurationError(f"{study_id}: annotation.method must be existing or celltypist")
        if annotation["method"] == "celltypist":
            levels = annotation.get("levels")
            if not isinstance(levels, list) or not levels or not set(levels).issubset({"l1", "l2", "l3"}):
                raise ConfigurationError(
                    f"{study_id}: celltypist annotation.levels must contain l1, l2, or l3"
                )
        conversion = study.get("conversion")
        if conversion is not None:
            if not isinstance(conversion, dict):
                raise ConfigurationError(f"{study_id}: conversion must be an object")
            required_conversion = {"source", "test_source", "output", "format", "assay"}
            missing_conversion = required_conversion - conversion.keys()
            if missing_conversion:
                raise ConfigurationError(
                    f"{study_id}: conversion missing keys {sorted(missing_conversion)}"
                )
            if not all(isinstance(conversion[key], str) and conversion[key] for key in required_conversion):
                raise ConfigurationError(f"{study_id}: conversion values must be non-empty strings")

    if len(test_inputs) != len(set(test_inputs)):
        raise ConfigurationError("Study test_input filenames must be unique")

    if input_sources is not None:
        _validate_input_sources(studies, input_sources)


def _validate_input_sources(studies: dict, input_sources: dict) -> None:
    if input_sources.get("schema_version") != 1:
        raise ConfigurationError("Only input_sources schema_version 1 is supported")
    artifacts = input_sources.get("artifacts")
    if not isinstance(artifacts, list):
        raise ConfigurationError("input_sources.artifacts must be a list")

    study_ids = set(studies["studies"])
    artifact_ids: set[str] = set()
    paths_by_study: dict[str, set[str]] = {study_id: set() for study_id in study_ids}
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ConfigurationError("input_sources artifacts must be objects")
        artifact_id = artifact.get("id")
        artifact_studies = artifact.get("studies")
        artifact_path = artifact.get("path")
        if not isinstance(artifact_id, str) or not artifact_id:
            raise ConfigurationError("input_sources artifact ids must be non-empty strings")
        if artifact_id in artifact_ids:
            raise ConfigurationError(f"Duplicate input_sources artifact id: {artifact_id}")
        artifact_ids.add(artifact_id)
        if not isinstance(artifact_studies, list) or not artifact_studies:
            raise ConfigurationError(f"{artifact_id}: studies must be a non-empty list")
        unknown = set(artifact_studies) - study_ids
        if unknown:
            raise ConfigurationError(f"{artifact_id}: unknown studies {sorted(unknown)}")
        if not isinstance(artifact_path, str) or not artifact_path:
            raise ConfigurationError(f"{artifact_id}: path must be a non-empty string")
        for study_id in artifact_studies:
            paths_by_study[study_id].add(artifact_path)

    for study_id, study in studies["studies"].items():
        if not paths_by_study[study_id]:
            raise ConfigurationError(f"{study_id}: no input_sources artifact is declared")
        expression_path = study.get("conversion", {}).get("source", study["input"])
        if expression_path not in paths_by_study[study_id]:
            raise ConfigurationError(
                f"{study_id}: expression input is absent from input_sources: {expression_path}"
            )
        missing_dependencies = set(study["preparation"].get("dependencies", [])) - paths_by_study[study_id]
        if missing_dependencies:
            raise ConfigurationError(
                f"{study_id}: preparation dependencies absent from input_sources: "
                f"{sorted(missing_dependencies)}"
            )


def config_digest(*documents: dict) -> str:
    payload = json.dumps(documents, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()

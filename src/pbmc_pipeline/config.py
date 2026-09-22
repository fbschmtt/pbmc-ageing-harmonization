from __future__ import annotations

import hashlib
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
    validate_configuration(pipeline, studies, schema)
    return pipeline, studies, schema


def validate_configuration(pipeline: dict, studies: dict, schema: dict) -> None:
    if pipeline.get("schema_version") != 1 or studies.get("schema_version") != 1:
        raise ConfigurationError("Only configuration schema_version 1 is supported")
    if not studies.get("studies"):
        raise ConfigurationError("No studies are configured")
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
        if not isinstance(study["preparation"].get("adapter"), str):
            raise ConfigurationError(f"{study_id}: preparation.adapter must be a string")
        if not isinstance(study["preparation"].get("dependencies", []), list):
            raise ConfigurationError(f"{study_id}: preparation.dependencies must be a list")
        if not isinstance(study["preparation"].get("allow_cell_subset", False), bool):
            raise ConfigurationError(f"{study_id}: preparation.allow_cell_subset must be a boolean")


def config_digest(*documents: dict) -> str:
    payload = json.dumps(documents, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()

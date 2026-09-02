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

    required = set(schema["required"])
    prediction_columns = {f"aifi_{level}_majority" for level in ("l1", "l2", "l3")}
    for study_id, study in studies["studies"].items():
        missing_keys = {"input", "test_input", "counts_source", "metadata", "annotation"} - study.keys()
        if missing_keys:
            raise ConfigurationError(f"{study_id}: missing keys {sorted(missing_keys)}")
        supplied = set(study["metadata"])
        if study["annotation"]["method"] == "celltypist":
            supplied |= prediction_columns
        missing_columns = required - supplied - {"study"}
        if missing_columns:
            raise ConfigurationError(
                f"{study_id}: homogeneous metadata is missing {sorted(missing_columns)}"
            )


def config_digest(*documents: dict) -> str:
    payload = json.dumps(documents, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


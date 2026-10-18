"""Materialize study-specific metadata into one canonical row per cell."""
from __future__ import annotations

import importlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .source import read_study_obs


class PreparationError(ValueError):
    """Raised when a prepared cell-metadata artifact is not identity-safe."""


PREDICTION_COLUMNS = {f"aifi_{level}_majority" for level in ("l1", "l2", "l3")}


def prepare_study(
    root: Path,
    study_id: str,
    study: dict[str, Any],
    schema: dict[str, Any],
    expression_path: Path,
    output_path: Path,
    report_path: Path,
    source_obs_output_path: Path | None = None,
) -> dict[str, Any]:
    """Run one named adapter and write canonical plus source-metadata cell tables."""
    source_obs = read_study_obs(root, expression_path, study)
    source_ids = pd.Index(source_obs.index.astype(str), name="cell_id")
    source_obs.index = source_ids

    adapter_name = study["preparation"]["adapter"]
    module = importlib.import_module(f"pbmc_pipeline.studies.{adapter_name}")
    prepared = module.prepare_cells(source_obs, root)
    prepared["study_doi"] = study["publication_doi"]
    allow_cell_subset = study["preparation"].get("allow_cell_subset", False)
    prepared = validate_prepared_cells(
        prepared, source_ids, schema, allow_cell_subset=allow_cell_subset
    )
    prepared["study_doi"] = pd.Categorical(
        prepared["study_doi"], categories=[study["publication_doi"]]
    )
    sample_batches = prepared[["sample", "batch_single_cell"]].drop_duplicates()
    fallback_prefix = f"{study_id}__sample__"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(output_path, index_label="cell_id", compression="infer")
    if source_obs_output_path is not None:
        source_obs_output_path.parent.mkdir(parents=True, exist_ok=True)
        source_obs.loc[prepared.index].to_csv(
            source_obs_output_path, index_label="cell_id", compression="infer"
        )
    report = {
        "study": study_id,
        "adapter": adapter_name,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "expression_input": str(expression_path),
        "output": str(output_path),
        "source_n_cells": len(source_ids),
        "n_cells": len(prepared),
        "excluded_n_cells": len(source_ids) - len(prepared),
        "columns": list(prepared.columns),
        "technical_batch": {
            "n_batches": int(sample_batches["batch_single_cell"].nunique()),
            "n_sample_specific_fallbacks": int(
                sample_batches["batch_single_cell"].astype("string")
                .str.startswith((fallback_prefix, f"{study_id}__sample_fallback__"))
                .sum()
            ),
        },
        "source_obs_output": None if source_obs_output_path is None else str(source_obs_output_path),
        "source_obs_columns": list(source_obs.columns),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def read_prepared_cells(
    path: Path,
    expression_ids,
    schema: dict[str, Any],
    *,
    allow_cell_subset: bool = False,
) -> pd.DataFrame:
    """Read and recheck a materialized metadata table before attaching it to counts."""
    prepared = pd.read_csv(path, index_col="cell_id")
    prepared.index = pd.Index(prepared.index.astype(str), name="cell_id")
    return validate_prepared_cells(
        prepared,
        pd.Index(expression_ids.astype(str), name="cell_id"),
        schema,
        allow_cell_subset=allow_cell_subset,
    )


def validate_prepared_cells(
    prepared,
    expression_ids,
    schema: dict[str, Any],
    *,
    allow_cell_subset: bool = False,
) -> pd.DataFrame:
    """Require exact cell coverage and normalize the canonical metadata contract."""
    if not isinstance(prepared, pd.DataFrame):
        raise PreparationError("Study adapter must return a pandas DataFrame")
    result = prepared.copy()
    result.index = pd.Index(result.index.astype(str), name="cell_id")
    if not result.index.is_unique:
        raise PreparationError("Prepared cell metadata contains duplicate cell_id values")
    if not expression_ids.is_unique:
        raise PreparationError("Expression input contains duplicate cell identifiers")
    extra_ids = result.index.difference(expression_ids)
    missing_ids = expression_ids.difference(result.index)
    if len(extra_ids) or (len(missing_ids) and not allow_cell_subset):
        raise PreparationError(
            "Prepared cell metadata must cover exactly the expression cells "
            f"(missing={len(missing_ids)}, extra={len(extra_ids)})"
        )
    selected_ids = expression_ids[expression_ids.isin(result.index)]
    result = result.loc[selected_ids].copy()
    required = set(schema["required"]) - PREDICTION_COLUMNS
    missing_columns = required - set(result.columns)
    if missing_columns:
        raise PreparationError(f"Prepared metadata is missing {sorted(missing_columns)}")
    for column in ("sample", "subject"):
        values = result[column]
        if values.isna().any() or values.astype("string").str.strip().eq("").any():
            raise PreparationError(f"Prepared metadata has missing biological identifier {column}")
    result = normalize_metadata_frame(result, schema)
    result = normalize_sample_nested_batch(result)
    batch_counts = result.groupby("sample", observed=True)["batch_single_cell"].nunique()
    if batch_counts.gt(1).any():
        raise PreparationError(
            "Technical batch assignment must be constant within each sample"
        )
    return result


def normalize_sample_nested_batch(frame: pd.DataFrame) -> pd.DataFrame:
    """Make batch_single_cell a sample-nested technical grouping.

    Keep source batch values when each sample has one unambiguous value and
    the source identifies at least one batch shared by multiple samples.
    Samples with missing or conflicting batch values receive a sample-specific
    fallback. If the source has no shared grouping (including a study-wide
    constant or one unique value per sample), every sample becomes its own
    batch.
    """
    result = frame.copy()
    sample = result["sample"].astype("string").str.strip()
    batch = result["batch_single_cell"].astype("string").str.strip()
    unknown = {"", "nan", "none", "null", "not_provided", "unknown", "unclear"}
    source_batch = batch.mask(batch.str.lower().isin(unknown))

    pairs = pd.DataFrame({"sample": sample, "batch": source_batch}).dropna().drop_duplicates()
    values_by_sample = pairs.groupby("sample", observed=True)["batch"].agg(list)
    assigned: dict[str, str] = {}
    fallback_samples: set[str] = set()
    existing_fallbacks: set[str] = set()
    for sample_id in sample.drop_duplicates().astype(str):
        values = values_by_sample.get(sample_id, [])
        if len(values) == 1:
            assigned[sample_id] = str(values[0])
            if str(values[0]) in {
                f"{result['study'].iloc[0]}__sample__{sample_id}",
                f"{result['study'].iloc[0]}__sample_fallback__{sample_id}",
            }:
                existing_fallbacks.add(sample_id)
        else:
            fallback_samples.add(sample_id)

    # A source field that never groups distinct samples contributes no shared
    # technical correlation information; represent those samples directly.
    assigned_counts = pd.Series(list(assigned.values()), dtype="string").value_counts()
    has_shared_source_batch = assigned_counts.gt(1).any()
    if not has_shared_source_batch:
        fallback_samples.update(assigned)

    study = str(result["study"].iloc[0])
    source_values = set(source_batch.dropna().astype(str))
    for sample_id in fallback_samples:
        if sample_id in existing_fallbacks:
            continue
        fallback = f"{study}__sample__{sample_id}"
        if fallback in source_values:
            fallback = f"{study}__sample_fallback__{sample_id}"
        assigned[sample_id] = fallback
    result["batch_single_cell"] = sample.map(assigned).astype("string")
    if result["batch_single_cell"].isna().any():
        raise PreparationError("Could not assign a technical batch to every sample")
    return result


def normalize_metadata_frame(frame: pd.DataFrame, schema: dict[str, Any]) -> pd.DataFrame:
    """Use one missing-value convention across prepared and harmonized metadata."""
    result = frame.copy()
    for column, expected_type in schema["required"].items():
        if column not in result:
            continue
        if expected_type == "string":
            result[column] = result[column].astype("string").fillna(schema["missing_string"])
        elif expected_type == "float":
            result[column] = pd.to_numeric(result[column], errors="raise").astype(float)
        else:
            raise PreparationError(f"Unsupported metadata type for {column}: {expected_type}")
    return result

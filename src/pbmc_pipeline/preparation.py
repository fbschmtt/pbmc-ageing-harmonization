"""Materialize study-specific metadata into one canonical row per cell."""
from __future__ import annotations

import importlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


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
) -> dict[str, Any]:
    """Run one named adapter and write its canonical one-row-per-cell table."""
    import anndata as ad

    source = ad.read_h5ad(expression_path, backed="r")
    try:
        source_obs = source.obs.copy()
        source_ids = pd.Index(source.obs_names.astype(str), name="cell_id")
    finally:
        source.file.close()

    adapter_name = study["preparation"]["adapter"]
    module = importlib.import_module(f"pbmc_pipeline.studies.{adapter_name}")
    prepared = module.prepare_cells(source_obs, root)
    allow_cell_subset = study["preparation"].get("allow_cell_subset", False)
    prepared = validate_prepared_cells(
        prepared, source_ids, schema, allow_cell_subset=allow_cell_subset
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(output_path, index_label="cell_id", compression="infer")
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
    return normalize_metadata_frame(result, schema)


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

from __future__ import annotations

from typing import Any


class ValidationError(ValueError):
    """Raised when an input or output violates the pipeline contract."""


COUNT_INTEGER_TOLERANCE = 1e-2


def validate_counts(adata) -> None:
    import numpy as np
    from scipy import sparse

    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
    if values.size == 0:
        raise ValidationError("Count matrix contains no non-zero values")
    if not np.isfinite(values).all() or values.min() < 0:
        raise ValidationError("Count matrix must be finite and nonnegative")
    max_integer_deviation = float(np.abs(values - np.rint(values)).max())
    if max_integer_deviation > COUNT_INTEGER_TOLERANCE:
        raise ValidationError(
            "Matrix is not count-like: maximum deviation from the nearest integer is "
            f"{max_integer_deviation:.3g}, exceeding {COUNT_INTEGER_TOLERANCE:g}"
        )


def validate_output(adata, schema: dict[str, Any]) -> dict[str, Any]:
    import pandas as pd

    if not adata.obs_names.is_unique:
        raise ValidationError("Cell identifiers are duplicated")
    if not adata.var_names.is_unique:
        raise ValidationError("Feature identifiers are duplicated")
    missing = set(schema["required"]) - set(adata.obs.columns)
    if missing:
        raise ValidationError(f"Required output columns are absent: {sorted(missing)}")
    for column in schema["non_nullable"]:
        if adata.obs[column].isna().any():
            raise ValidationError(f"Required identifier {column} contains missing values")
    for column, allowed in schema.get("controlled_values", {}).items():
        observed = set(adata.obs[column].dropna().astype(str).unique())
        unexpected = observed - set(allowed)
        if unexpected:
            raise ValidationError(f"{column} has unexpected values: {sorted(unexpected)}")
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "missing_by_column": {
            column: int(pd.isna(adata.obs[column]).sum()) for column in schema["required"]
        },
    }

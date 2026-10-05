"""Shared normalization and presentation rules for categorical covariates."""
from __future__ import annotations

import pandas as pd


def normalize_cmv_status(values: pd.Series) -> pd.Series:
    """Normalize observed negative/positive CMV labels to no/yes."""
    normalized = values.astype("string").str.strip().str.lower()
    return normalized.replace({"negative": "no", "positive": "yes"})


def categorical_levels(covariate: str, values: pd.Series) -> list[str]:
    """Return treatment levels with the intended reference level first."""
    levels = sorted(values.astype("string").dropna().unique().tolist())
    if covariate == "age_bin":
        levels.sort(key=lambda value: int(value.split("-", maxsplit=1)[0]))
        return levels
    preferred = {"sex": "female", "cmv": "no"}.get(covariate)
    if preferred not in levels and covariate == "cmv" and "negative" in levels:
        preferred = "negative"
    if preferred in levels:
        return [preferred, *(level for level in levels if level != preferred)]
    return levels


def categorical_contrast_label(covariate: str, level: str, reference: str) -> str:
    """Format a categorical coefficient with its explicit reference group."""
    return f"{covariate.upper()}: {level} vs {reference}"

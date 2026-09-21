"""Pure data-preparation helpers used by executable QC reports."""
from __future__ import annotations

import pandas as pd


def gene_presence_indicators(var: pd.DataFrame) -> pd.DataFrame:
    """Return study-by-gene presence booleans from merge provenance columns."""
    columns = sorted(column for column in var if column.startswith("available_in_"))
    if not columns:
        return pd.DataFrame(index=var.index)
    indicators = var[columns].fillna(False).astype(bool).copy()
    indicators.columns = [column.removeprefix("available_in_") for column in columns]
    return indicators


def aifi_l2_concordance(obs: pd.DataFrame) -> pd.DataFrame:
    """Cross-tab original per-study and newly merged AIFI-L2 labels."""
    required = {"aifi_l2_study_majority", "aifi_l2_majority"}
    if not required.issubset(obs):
        return pd.DataFrame()
    return pd.crosstab(
        obs["aifi_l2_study_majority"].astype(str),
        obs["aifi_l2_majority"].astype(str),
        normalize="index",
    )


def pseudobulk_celltype_fractions(obs: pd.DataFrame) -> pd.DataFrame:
    """Return within-study cell fractions, weighting pseudobulk rows by ``n_cells``."""
    required = {"study", "aifi_l2_majority", "n_cells"}
    missing = required - set(obs)
    if missing:
        raise ValueError(f"Pseudobulk observations are missing columns: {sorted(missing)}")
    composition = pd.pivot_table(
        obs.assign(
            study=obs["study"].astype(str),
            aifi_l2_majority=obs["aifi_l2_majority"].astype(str),
        ),
        index="study",
        columns="aifi_l2_majority",
        values="n_cells",
        aggfunc="sum",
        fill_value=0,
    )
    return composition.div(composition.sum(axis=1), axis=0)

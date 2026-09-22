from __future__ import annotations

import pandas as pd

from ..metadata import MetadataError, build_homogeneous_obs


def safe_left_join(cells: pd.DataFrame, table: pd.DataFrame, on: list[str]) -> pd.DataFrame:
    """Join sample-level data without changing source-cell order or cardinality."""
    if table.duplicated(on).any():
        count = int(table.duplicated(on, keep=False).sum())
        raise MetadataError(f"Auxiliary metadata has {count} duplicate rows for {on}")
    index_name = cells.index.name or "cell_id"
    joined = cells.reset_index(names=index_name).merge(
        table, on=on, how="left", validate="many_to_one", suffixes=("", "_joined")
    ).set_index(index_name)
    if not joined.index.equals(cells.index):
        raise MetadataError("Auxiliary metadata join changed cell order or identity")
    return joined


def canonical(source: pd.DataFrame, study_id: str, rules: dict) -> pd.DataFrame:
    return build_homogeneous_obs(source, study_id, rules)

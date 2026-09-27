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
    """Cross-tab per-study AIFI-L2 labels against experimental merged labels."""
    required = {"aifi_l2_majority", "experimental_aifi_l2_majority"}
    if not required.issubset(obs):
        return pd.DataFrame()
    return pd.crosstab(
        obs["aifi_l2_majority"].astype(str),
        obs["experimental_aifi_l2_majority"].astype(str),
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


def sample_cell_type_fractions(
    adata, *, denominator: str, fraction_name: str, min_denominator_cells: int = 1
) -> pd.DataFrame:
    """Summarize eligible samples for one split type against a cell-count denominator."""
    if not isinstance(min_denominator_cells, int) or min_denominator_cells < 1:
        raise ValueError("min_denominator_cells must be a positive integer")
    required = {"study", "sample", "age", denominator}
    missing = required - set(adata.obs)
    if missing:
        raise ValueError(f"Sample cell-type fractions require columns: {sorted(missing)}")
    grouped = adata.obs.groupby(["study", "sample"], observed=True)
    result = grouped.agg(
        age=("age", "first"),
        n_cells_type=("study", "size"),
        denominator_cells=(denominator, "first"),
        n_unique_ages=("age", "nunique"),
        n_unique_denominators=(denominator, "nunique"),
    ).reset_index()
    if (result["n_unique_ages"] > 1).any() or (result["n_unique_denominators"] > 1).any():
        raise ValueError("Age and sample cell count must be constant within each study × sample")
    result[fraction_name] = result["n_cells_type"] / result["denominator_cells"]
    if (result[fraction_name] > 1).any():
        raise ValueError("Cell-type count exceeds its configured denominator")
    return result.loc[
        result["denominator_cells"] >= min_denominator_cells,
    ].drop(columns=["n_unique_ages", "n_unique_denominators"])


def sample_cluster_fractions(adata, *, min_denominator_cells: int = 1) -> pd.DataFrame:
    """Return eligible samples' local-cluster fractions within their split cell type."""
    if not isinstance(min_denominator_cells, int) or min_denominator_cells < 1:
        raise ValueError("min_denominator_cells must be a positive integer")
    required = {"study", "sample", "age", "cluster"}
    missing = required - set(adata.obs)
    if missing:
        raise ValueError(f"Sample cluster fractions require columns: {sorted(missing)}")
    sample_keys = ["study", "sample", "age"]
    counts = (
        adata.obs.groupby([*sample_keys, "cluster"], observed=True)
        .size()
        .rename("n_cells_cluster")
        .reset_index()
    )
    totals = (
        adata.obs.groupby(sample_keys, observed=True)
        .size()
        .rename("n_cells_type")
        .reset_index()
    )
    result = counts.merge(totals, on=sample_keys, validate="many_to_one")
    result["fraction_within_cell_type"] = result["n_cells_cluster"] / result["n_cells_type"]
    return result.loc[result["n_cells_type"] >= min_denominator_cells]

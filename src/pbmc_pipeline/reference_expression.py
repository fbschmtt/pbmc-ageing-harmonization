"""Descriptive study-balanced reference-bin pseudobulk expression."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

_UNKNOWN = {"", "not_provided", "heterogeneous during bulk", "nan", "none", "unknown"}


def reference_expression(adata, spec: dict, *, cell_type: str, test_mode: bool = False):
    """Return per-study and equal-study reference-bin summaries in CPM."""
    obs = adata.obs
    inclusion = spec["sample_inclusion"]
    settings = spec["age_trajectory"]
    start = settings["reference_bin_start_age"]
    end = start + settings["bin_width_years"]
    study = obs["study"].astype("string").str.strip()
    sample = obs["sample"].astype("string").str.strip()
    age = pd.to_numeric(obs["age"], errors="coerce")
    cells = pd.to_numeric(obs["n_cells"], errors="coerce")
    totals = pd.to_numeric(obs["total_counts"], errors="coerce")
    valid = (
        study.notna() & ~study.str.lower().isin(_UNKNOWN)
        & sample.notna() & ~sample.str.lower().isin(_UNKNOWN)
        & np.isfinite(age) & (age >= start) & (age < end)
        & np.isfinite(cells) & (cells > 0)
        & np.isfinite(totals) & (totals > 0)
    )
    if not test_mode:
        valid &= cells >= inclusion["minimum_cells_per_pseudobulk_inclusive"]
    exclusions = inclusion.get("exclude_subject_prefixes_by_study", {})
    for name, prefixes in exclusions.items():
        if (study == name).fillna(False).any():
            if "subject" not in obs:
                raise ValueError("subject metadata is required to apply configured exclusions")
            subject = obs["subject"].astype("string").str.strip().str.upper()
            valid &= ~((study == name) & subject.str.startswith(tuple(p.upper() for p in prefixes), na=False)).fillna(False)
    selected = obs.loc[valid, ["study", "sample"]].copy()
    if selected.duplicated(["study", "sample"]).any():
        raise ValueError("Reference expression requires one pseudobulk per study/sample and cell type")
    genes = adata.var_names.astype(str).to_numpy()
    rows = []
    for name in sorted(study[valid].unique()):
        column = f"available_in_{name}"
        if column not in adata.var:
            raise ValueError(f"Reference expression requires gene availability column {column}")
        available = adata.var[column].fillna(False).astype(bool).to_numpy()
        mask = (valid & (study == name)).to_numpy(dtype=bool)
        matrix = adata.X[mask, :]
        matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
        if not np.isfinite(matrix).all() or (matrix < 0).any():
            raise ValueError("Reference expression counts must be finite and non-negative")
        library = matrix.sum(axis=1)
        if (library <= 0).any():
            raise ValueError("Reference expression requires positive measured library totals")
        mean_log = np.log1p(matrix[:, available] / library[:, None] * 1e6).mean(axis=0)
        rows.append(pd.DataFrame({
            "cell_type": cell_type, "gene": genes[available], "study": name,
            "reference_bin": f"{start}-{end}", "n_samples": int(mask.sum()),
            "reference_cpm": np.expm1(mean_log),
        }))
    columns = ["cell_type", "gene", "study", "reference_bin", "n_samples", "reference_cpm"]
    per_study = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=columns)
    combined_columns = ["cell_type", "gene", "reference_bin", "n_studies", "n_samples", "reference_cpm"]
    if per_study.empty:
        return per_study, pd.DataFrame(columns=combined_columns)
    summary = per_study.assign(_log=np.log1p(per_study["reference_cpm"])).groupby(
        ["cell_type", "gene", "reference_bin"], as_index=False, sort=True,
    ).agg(n_studies=("study", "size"), n_samples=("n_samples", "sum"), _log=("_log", "mean"))
    summary["reference_cpm"] = np.expm1(summary.pop("_log"))
    return per_study, summary[combined_columns]

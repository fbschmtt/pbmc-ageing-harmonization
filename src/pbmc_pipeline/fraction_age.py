"""Study-specific, reference-aligned age-bin models of sample cell fractions."""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .config import AgeTrajectorySettings
from .covariates import categorical_levels, normalize_cmv_status, normalize_smoking_status

_UNKNOWN = {"", "not_provided", "heterogeneous during bulk", "unknown", "nan", "none"}
_COLUMNS = [
    "study", "denominator", "age_bin", "reference_bin", "n_samples_bin", "n_samples_model",
    "estimate", "ci_low", "ci_high", "pvalue", "reference_mean_fraction",
    "reference_adjusted_fraction", "covariates", "design",
]


def fit_fraction_age_bins(
    abundance: pd.DataFrame, *, settings: AgeTrajectorySettings,
    min_denominator_cells: int, minimum_samples_per_bin: int,
) -> tuple[pd.DataFrame, list[dict]]:
    """Fit separate adult OLS fraction models and align each reference bin to zero.

    A row is one study/sample, including zero counts. Each model uses all
    covariates with 95% coverage and variation, and complete cases thereafter.
    Standard OLS intervals are pointwise, not adjusted across contrasts.
    """
    required = {"study", "sample", "age", "n_cells", "n_cells_in_sample",
                "n_cells_in_sample_l1_parent", "study_site"}
    missing = required - set(abundance)
    if missing:
        return pd.DataFrame(columns=_COLUMNS), [{
            "status": "skipped", "reason": f"missing abundance columns: {sorted(missing)}; "
            "rerun DE to regenerate abundance input",
        }]
    if abundance.duplicated(["study", "sample"]).any():
        raise ValueError("Abundance input must contain one row per study and sample")
    df = abundance.copy()
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df["n_cells"] = pd.to_numeric(df["n_cells"], errors="coerce")
    for column in ("n_cells_in_sample", "n_cells_in_sample_l1_parent"):
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if "cmv" in df:
        df["cmv"] = normalize_cmv_status(df["cmv"])
    if "smoking_status" in df:
        df["smoking_status"] = normalize_smoking_status(df["smoking_status"])
    ref_start = settings.reference_bin_start_age
    width = settings.bin_width_years
    reference = f"{ref_start}-{ref_start + width}"
    df = df.loc[df["age"].ge(ref_start) & df["age"].lt(settings.strict_age_cutoff_exclusive)].copy()
    df["age_start"] = (np.floor(df["age"] / width) * width).astype(int)
    candidates = ("study_site", "sex", "bmi", "cmv", "smoking_status")
    results: list[dict] = []
    statuses: list[dict] = []
    for denominator, denom_column in (
        ("retained_pbmc", "n_cells_in_sample"),
        ("l1_parent", "n_cells_in_sample_l1_parent"),
    ):
        for study, group in df.groupby("study", sort=True, observed=True):
            group = group.loc[
                group[denom_column].ge(min_denominator_cells)
                & group["n_cells"].ge(0) & group["n_cells"].le(group[denom_column])
            ].copy()
            status = {"study": str(study), "denominator": denominator, "reference_bin": reference}
            if group.empty:
                statuses.append({**status, "status": "skipped", "reason": "no eligible denominator"})
                continue
            selected = []
            valid_masks = {}
            for column in candidates:
                if column not in group:
                    continue
                if column == "bmi":
                    values = pd.to_numeric(group[column], errors="coerce")
                    valid = values.notna() & np.isfinite(values)
                else:
                    values = group[column].astype("string").str.strip().str.lower()
                    valid = values.notna() & ~values.isin(_UNKNOWN)
                if valid.mean() >= 0.95 and values.loc[valid].nunique() >= 2:
                    selected.append(column)
                    valid_masks[column] = valid
                    group[column] = values
            for column in selected:
                group = group.loc[valid_masks[column]].copy()
            counts = group["age_start"].value_counts().to_dict()
            starts = sorted(int(start) for start, count in counts.items()
                            if count >= minimum_samples_per_bin)
            support = {f"{int(start)}-{int(start) + width}": int(count)
                       for start, count in sorted(counts.items())}
            status.update({"bin_counts": support, "covariates": selected,
                           "minimum_samples_per_bin": minimum_samples_per_bin})
            if ref_start not in starts or len(starts) < 2:
                statuses.append({**status, "status": "skipped", "reason":
                                 "reference bin and at least one other bin must meet minimum support"})
                continue
            group = group.loc[group["age_start"].isin(starts)].copy()
            pieces = [pd.Series(1.0, index=group.index, name="intercept")]
            for start in starts[1:]:
                pieces.append((group["age_start"] == start).astype(float).rename(f"age_bin_{start}"))
            for column in selected:
                if column == "bmi":
                    pieces.append(pd.to_numeric(group[column]).rename(column))
                else:
                    values = group[column].astype("string")
                    levels = categorical_levels(column, values)
                    categorical = pd.Series(pd.Categorical(values, categories=levels), index=group.index)
                    dummies = pd.get_dummies(categorical, prefix=column, drop_first=True, dtype=float)
                    pieces.extend(dummies[col].rename(col) for col in dummies)
            design = pd.concat(pieces, axis=1).astype(float)
            if len(group) <= design.shape[1] or np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
                statuses.append({**status, "status": "skipped", "reason": "rank-deficient design or insufficient residual degrees of freedom"})
                continue
            fit = sm.OLS((group["n_cells"] / group[denom_column]).to_numpy(), design).fit()
            intervals = fit.conf_int()
            design_label = " + ".join(design.columns)
            reference_rows = group.loc[group["age_start"] == ref_start]
            reference_mean = float((reference_rows["n_cells"] / reference_rows[denom_column]).mean())
            # Average the predictions at the reference bin over the complete
            # model cohort. This retains mean BMI and observed category/site
            # proportions, rather than setting treatment indicators to zero.
            reference_design = design.copy()
            reference_design.loc[:, [f"age_bin_{start}" for start in starts[1:]]] = 0.0
            reference_adjusted = float(fit.predict(reference_design).mean())
            for start in starts:
                term = f"age_bin_{start}"
                is_ref = start == ref_start
                results.append({
                    "study": str(study), "denominator": denominator,
                    "age_bin": f"{start}-{start + width}", "reference_bin": reference,
                    "n_samples_bin": counts[start], "n_samples_model": len(group),
                    "estimate": 0.0 if is_ref else float(fit.params[term]),
                    "ci_low": 0.0 if is_ref else float(intervals.loc[term, 0]),
                    "ci_high": 0.0 if is_ref else float(intervals.loc[term, 1]),
                    "pvalue": np.nan if is_ref else float(fit.pvalues[term]),
                    "reference_mean_fraction": reference_mean if is_ref else np.nan,
                    "reference_adjusted_fraction": reference_adjusted if is_ref else np.nan,
                    "covariates": ";".join(selected), "design": design_label,
                })
            statuses.append({**status, "status": "complete", "n_samples_model": len(group)})
    return pd.DataFrame(results, columns=_COLUMNS), statuses

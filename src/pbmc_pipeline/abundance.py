"""Adult sample-fraction models and conditional binomial sampling diagnostics."""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .covariates import categorical_levels, normalize_cmv_status, normalize_smoking_status


def _bootstrap_binomial_residual_sd(design, denominators, fitted_probabilities, *, rng, n_bootstrap):
    """Estimate post-model residual SD expected from binomial cell sampling alone.

    The fitted probabilities and design are fixed inputs.  Thus this is a
    conditional measurement-noise reference: only finite cell sampling varies
    between replicates, while the observed covariate structure is retained.
    """
    n_samples, n_terms = design.shape
    if n_samples <= n_terms:
        return np.nan
    denominators = np.asarray(denominators, dtype=np.int64)
    fitted_probabilities = np.asarray(fitted_probabilities, dtype=float)
    simulated = rng.binomial(
        denominators, fitted_probabilities, size=(n_bootstrap, n_samples),
    ) / denominators
    orthonormal_design, _ = np.linalg.qr(design, mode="reduced")
    # Refit the same fixed design because the observed quantity is also a
    # post-covariate-model residual SD.  This removes chance alignment between
    # simulated binomial noise and the covariates, just as in the observed fit.
    residuals = simulated - (simulated @ orthonormal_design) @ orthonormal_design.T
    return float(np.sqrt(np.mean(np.sum(residuals**2, axis=1) / (n_samples - n_terms))))


def fit_parent_fraction_models(
    fractions, adata, *, fraction_column, bootstrap_replicates, bootstrap_seed,
):
    """Fit adult-only, study-specific OLS models and residual diagnostics."""
    candidate_covariates = ["age", "sex", "bmi", "cmv", "smoking_status"]
    present = [column for column in candidate_covariates if column in adata.obs]
    metadata_columns = [column for column in present if column != "age"]
    metadata = adata.obs.groupby(["study", "sample"], observed=True)[metadata_columns].first().reset_index()
    # Repeated cell metadata must describe one sample consistently.
    for column in present:
        n_unique = adata.obs.groupby(["study", "sample"], observed=True)[column].nunique(dropna=False)
        if (n_unique > 1).any():
            raise ValueError(f"{column} must be constant within each study × sample")
    data = fractions.merge(metadata, on=["study", "sample"], how="left", validate="one_to_one")
    data["age"] = pd.to_numeric(data["age"], errors="coerce")
    if "cmv" in present:
        data["cmv"] = normalize_cmv_status(data["cmv"])
    if "smoking_status" in present:
        data["smoking_status"] = normalize_smoking_status(data["smoking_status"])
    data = data.loc[data["age"] >= 20].copy()
    results = []
    residual_diagnostics = []
    unavailable = {
        "", "not_provided", "heterogeneous during bulk", "unknown", "nan", "none",
    }
    bootstrap_rng = np.random.default_rng(bootstrap_seed)

    def valid_values(values, column):
        if column == "smoking_status":
            return normalize_smoking_status(values).notna()
        if column in {"age", "bmi"}:
            numeric = pd.to_numeric(values, errors="coerce")
            return numeric.notna() & np.isfinite(numeric)
        normalized = values.astype("string").str.strip().str.lower()
        return normalized.notna() & ~normalized.isin(unavailable)

    for study, subset in data.groupby("study", sort=True, observed=True):
        eligible = subset.copy()
        selected = []
        for column in candidate_covariates:
            if column not in present:
                continue
            valid = valid_values(eligible[column], column)
            coverage = float(valid.mean()) if len(eligible) else np.nan
            if coverage >= 0.95:
                if column in {"sex", "cmv", "smoking_status"}:
                    levels = categorical_levels(column, eligible.loc[valid, column])
                    variable = len(levels) >= 2
                else:
                    variable = pd.to_numeric(eligible.loc[valid, column], errors="coerce").nunique() >= 2
                if variable:
                    selected.append(column)

        if not selected:
            continue
        valid_rows = pd.Series(True, index=eligible.index)
        for column in selected:
            valid_rows &= valid_values(eligible[column], column)
        model_data = eligible.loc[valid_rows].copy()
        if model_data.empty:
            continue
        pieces = [pd.Series(1.0, index=model_data.index, name="intercept")]
        term_covariates = {"intercept": "intercept"}
        for column in selected:
            if column in {"age", "bmi"}:
                pieces.append(pd.to_numeric(model_data[column], errors="coerce").rename(column))
                term_covariates[column] = column
            else:
                values = model_data[column].astype("string").str.strip().str.lower()
                levels = categorical_levels(column, values)
                values = pd.Series(pd.Categorical(values, categories=levels), index=model_data.index)
                dummies = pd.get_dummies(values, prefix=column, drop_first=True, dtype=float)
                for term in dummies:
                    pieces.append(dummies[term].rename(term))
                    term_covariates[term] = column
        design = pd.concat(pieces, axis=1).astype(float)
        if len(model_data) <= design.shape[1] or np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
            continue
        design_array = design.to_numpy()
        fit = sm.OLS(model_data[fraction_column].to_numpy(dtype=float), design_array).fit()
        confidence = fit.conf_int()
        mean_fraction = float(model_data[fraction_column].mean())
        unbounded_fitted_probabilities = np.asarray(fit.predict(design_array), dtype=float)
        n_clipped_probabilities = int(
            ((unbounded_fitted_probabilities < 0) | (unbounded_fitted_probabilities > 1)).sum()
        )
        # OLS is retained as the reported fraction model.  Its fitted values
        # provide donor-specific conditional means for the binomial null; clip
        # only because binomial probabilities must lie in [0, 1].  The count is
        # published so any material boundary correction is visible.
        fitted_probabilities = np.clip(unbounded_fitted_probabilities, 0, 1)
        residual_diagnostics.append({
            "study": str(study),
            "residual_sd": float(np.sqrt(fit.mse_resid)),
            "binomial_sampling_sd": _bootstrap_binomial_residual_sd(
                design_array,
                model_data["denominator_cells"].to_numpy(),
                fitted_probabilities,
                rng=bootstrap_rng,
                n_bootstrap=bootstrap_replicates,
            ),
            "n_samples": len(model_data),
            "mean_fraction": mean_fraction,
            "n_clipped_binomial_probabilities": n_clipped_probabilities,
            "fraction_clipped_binomial_probabilities": n_clipped_probabilities / len(model_data),
            "covariates": ";".join(selected),
            "model_terms": ";".join(design.columns.drop("intercept")),
        })
        for index, term in enumerate(design.columns):
            if term == "intercept":
                continue
            covariate = term_covariates[term]
            if covariate in {"sex", "cmv", "smoking_status"}:
                levels = categorical_levels(
                    covariate,
                    model_data[covariate].astype("string").str.strip().str.lower(),
                )
                contrast = f"{term.removeprefix(f'{covariate}_')} vs {levels[0]}"
            else:
                contrast = "per decade" if covariate == "age" else "per BMI unit"
            effect_scale = 10 if covariate == "age" else 1
            results.append({
                "study": str(study), "covariate": covariate, "term": term,
                "contrast": contrast, "estimate": fit.params[index] * effect_scale,
                "ci_low": confidence[index, 0] * effect_scale,
                "ci_high": confidence[index, 1] * effect_scale,
                "n_samples": len(model_data),
            })
    coefficient_columns = [
        "study", "covariate", "term", "contrast", "estimate", "ci_low", "ci_high",
        "n_samples",
    ]
    diagnostic_columns = [
        "study", "residual_sd", "binomial_sampling_sd", "n_samples", "mean_fraction",
        "n_clipped_binomial_probabilities", "fraction_clipped_binomial_probabilities",
        "covariates", "model_terms",
    ]
    # Keep the report-facing schema stable even when no study has enough adult
    # samples or covariate variation for an estimable model.
    return (
        pd.DataFrame(results, columns=coefficient_columns),
        pd.DataFrame(residual_diagnostics, columns=diagnostic_columns),
    )


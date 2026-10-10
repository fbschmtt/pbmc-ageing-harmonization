"""Refit the merged DESeq2 model for one cell type and profile one gene."""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from pbmc_pipeline.config import read_json
from pbmc_pipeline.covariates import categorical_levels
from pbmc_pipeline.differential_expression import (
    _UNKNOWN,
    _availability_mask,
    _has_covariate_variation,
    _load_pseudobulk,
    _matrix_to_counts,
    _prepare_metadata,
)
from pbmc_pipeline.gene_profile_report import generate_gene_profile_report
from pbmc_pipeline.utils import slugify


def _ols_result_table(fit, term_names: list[str]) -> pd.DataFrame:
    confidence_intervals = np.asarray(fit.conf_int(), dtype=float)
    return pd.DataFrame({
        "term": term_names,
        "coefficient": np.asarray(fit.params, dtype=float),
        "standard_error": np.asarray(fit.bse, dtype=float),
        "t_value": np.asarray(fit.tvalues, dtype=float),
        "p_value": np.asarray(fit.pvalues, dtype=float),
        "ci_lower": confidence_intervals[:, 0],
        "ci_upper": confidence_intervals[:, 1],
    })


def _in_sample_r_squared(observed, fitted) -> float | None:
    observed = np.asarray(observed, dtype=float)
    fitted = np.asarray(fitted, dtype=float)
    valid = np.isfinite(observed) & np.isfinite(fitted)
    observed = observed[valid]
    fitted = fitted[valid]
    total_sum_squares = float(np.sum(np.square(observed - observed.mean())))
    if len(observed) < 2 or total_sum_squares <= 0:
        return None
    residual_sum_squares = float(np.sum(np.square(observed - fitted)))
    return 1.0 - residual_sum_squares / total_sum_squares


def _residual_sd(residuals) -> float | None:
    residuals = np.asarray(residuals, dtype=float)
    residuals = residuals[np.isfinite(residuals)]
    return float(np.std(residuals, ddof=1)) if len(residuals) > 1 else None


def _fit_profile(
    *, input_path: Path, config_path: Path, cell_type: str, gene: str,
    output_dir: Path, cpus: int,
) -> Path:
    import statsmodels.api as sm
    from formulaic_contrasts import FormulaicContrasts
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    pipeline = read_json(config_path)
    fit_diagnostics_rows = []
    adata, spec, labels = _load_pseudobulk(input_path, pipeline)
    available_types = sorted(labels.astype(str).unique())
    if cell_type not in available_types:
        matches = [label for label in available_types if slugify(label) == slugify(cell_type)]
        if len(matches) == 1:
            cell_type = matches[0]
    if cell_type not in available_types:
        raise ValueError(f"Cell type {cell_type!r} is not present in {input_path}")
    if gene not in adata.var_names:
        raise ValueError(f"Gene {gene!r} is not present in {input_path}")

    type_positions = np.flatnonzero(labels.astype(str).to_numpy() == cell_type)
    output_dir = output_dir / slugify(cell_type) / slugify(gene)
    type_adata = adata[type_positions, :]
    obs = type_adata.obs.copy()
    all_type_total_counts = pd.to_numeric(
        obs["total_counts"], errors="coerce",
    ).to_numpy(dtype=float)
    valid_all_type_totals = (
        np.isfinite(all_type_total_counts) & (all_type_total_counts > 0)
    )
    if not valid_all_type_totals.any():
        raise ValueError("No positive finite total_counts are available for this cell type")
    all_samples_geometric_mean_total_counts = float(np.exp(
        np.mean(np.log(all_type_total_counts[valid_all_type_totals]))
    ))
    metadata, exclusions, _ = _prepare_metadata(obs, spec["sample_inclusion"])
    if metadata.empty:
        raise ValueError(f"No samples for {cell_type!r} pass configured inclusion filters")

    configured_design = spec["merged_design"]
    design_terms = [
        term.strip()
        for term in configured_design.removeprefix("~").split("+")
        if term.strip()
    ]
    core_covariates = [
        covariate for covariate in ("age", "sex", "log10_total_counts")
        if _has_covariate_variation(metadata, covariate)
    ] or ["age", "sex"]
    required_terms = [term for term in design_terms if term != "study_site"]
    valid = pd.Series(True, index=metadata.index)
    for term in design_terms:
        values = metadata[term]
        if term in {"age", "log10_total_counts", "bmi"}:
            values = pd.to_numeric(values, errors="coerce")
            valid &= values.notna() & np.isfinite(values)
        else:
            values = values.astype("string").str.strip().str.lower()
            valid &= values.notna() & ~values.isin(_UNKNOWN)
    metadata = metadata.loc[valid].copy()
    if metadata.empty:
        raise ValueError("No samples have complete metadata for the configured merged design")
    for covariate in core_covariates:
        if not _has_covariate_variation(metadata, covariate):
            raise ValueError(f"The merged design is not estimable: {covariate!r} does not vary")

    while True:
        studies = sorted(metadata["study"].astype(str).unique())
        gene_mask = _availability_mask(
            type_adata.var, studies, context=f"{cell_type}/merged"
        )
        fit_genes = type_adata.var_names[gene_mask]
        if gene not in fit_genes:
            raise ValueError(
                f"Gene {gene!r} is not available in every study used by the merged model "
                f"({', '.join(studies)})"
            )
        row_positions = obs.index.get_indexer(metadata.index)
        col_positions = type_adata.var_names.get_indexer(fit_genes)
        counts = _matrix_to_counts(
            type_adata.X[row_positions, :][:, col_positions],
            context=f"{cell_type}/merged profile",
        )
        nonzero_genes = counts.sum(axis=0) > 0
        counts = counts[:, nonzero_genes]
        fit_genes = fit_genes[nonzero_genes]
        if gene not in fit_genes:
            raise ValueError(f"Gene {gene!r} has zero counts across included samples")
        nonzero_libraries = counts.sum(axis=1) > 0
        if nonzero_libraries.all():
            break
        metadata = metadata.iloc[np.flatnonzero(nonzero_libraries)].copy()
        if metadata.empty:
            raise ValueError("All samples have zero counts over the merged model gene universe")

    varying_terms = [
        term for term in required_terms
        if term in core_covariates or _has_covariate_variation(metadata, term)
    ]
    use_site = "study_site" in design_terms and metadata["study_site"].nunique() > 1
    design = "~ " + " + ".join((["study_site"] if use_site else []) + varying_terms)

    model_metadata = metadata.copy()
    for column in ("study_site", "sex", "cmv"):
        if column in design and column in model_metadata:
            model_metadata[column] = pd.Categorical(
                model_metadata[column],
                categories=categorical_levels(column, model_metadata[column]),
            )
    design_matrix = FormulaicContrasts(model_metadata, design).design_matrix
    numeric_design = design_matrix.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(numeric_design))
    if rank < numeric_design.shape[1] or numeric_design.shape[0] <= rank:
        raise ValueError(
            f"Merged design {design!r} is not estimable for {len(model_metadata)} samples"
        )
    gene_position = fit_genes.get_loc(gene)

    # Normalize to the geometric mean of the included samples' full-library
    # totals, then fit log(normalized count + 1) with the same design matrix.
    total_counts = pd.to_numeric(model_metadata["total_counts"], errors="coerce").to_numpy(
        dtype=float,
    )
    if not np.isfinite(total_counts).all() or (total_counts <= 0).any():
        raise ValueError("Naive normalization requires positive finite sample total_counts")
    geometric_mean_total_counts = float(np.exp(np.mean(np.log(total_counts))))
    naive_size_factors = total_counts / geometric_mean_total_counts
    naive_normalized_count = counts[:, gene_position] / naive_size_factors
    naive_log_expression = np.log(naive_normalized_count + 1.0)
    naive_fit = sm.OLS(naive_log_expression, numeric_design).fit()
    naive_residuals = np.asarray(naive_fit.resid, dtype=float)
    fit_diagnostics_rows.append({
        "model": "Naive linear age model",
        "n_samples": len(naive_log_expression),
        "residual_sd": _residual_sd(naive_residuals),
        "residual_sd_scale": "log(normalized count + 1)",
        "r_squared": float(naive_fit.rsquared),
        "r_squared_label": "In-sample R²",
        "fixed_effects_r_squared": np.nan,
    })
    naive_model_result = _ols_result_table(
        naive_fit, design_matrix.columns.astype(str).tolist(),
    )

    # Fit an age-bin model on all otherwise eligible samples, including those
    # below the configured continuous-age model's minimum age.
    age_bin_inclusion = {
        **spec["sample_inclusion"], "minimum_age_years_inclusive": -np.inf,
    }
    age_bin_metadata, age_bin_exclusions, _ = _prepare_metadata(obs, age_bin_inclusion)
    if age_bin_metadata.empty:
        raise ValueError("No samples have valid metadata for the age-bin model")
    bin_width = int(pipeline["differential_expression"]["age_trajectory"]["bin_width_years"])
    age_bin_starts = (
        np.floor(age_bin_metadata["age"] / bin_width).astype(int) * bin_width
    )
    age_bin_metadata["age_bin"] = age_bin_starts.map(
        lambda start: f"{start}-{start + bin_width}"
    )

    while True:
        age_bin_studies = sorted(age_bin_metadata["study"].astype(str).unique())
        age_bin_gene_mask = _availability_mask(
            type_adata.var, age_bin_studies, context=f"{cell_type}/age-bin profile",
        )
        age_bin_genes = type_adata.var_names[age_bin_gene_mask]
        if gene not in age_bin_genes:
            raise ValueError(
                f"Gene {gene!r} is not available in every study used by the age-bin model "
                f"({', '.join(age_bin_studies)})"
            )
        age_bin_rows = obs.index.get_indexer(age_bin_metadata.index)
        age_bin_columns = type_adata.var_names.get_indexer(age_bin_genes)
        age_bin_counts = _matrix_to_counts(
            type_adata.X[age_bin_rows, :][:, age_bin_columns],
            context=f"{cell_type}/age-bin profile",
        )
        nonzero_age_bin_genes = age_bin_counts.sum(axis=0) > 0
        age_bin_counts = age_bin_counts[:, nonzero_age_bin_genes]
        age_bin_genes = age_bin_genes[nonzero_age_bin_genes]
        if gene not in age_bin_genes:
            raise ValueError(f"Gene {gene!r} has zero counts in the age-bin model samples")
        nonzero_age_bin_libraries = age_bin_counts.sum(axis=1) > 0
        if nonzero_age_bin_libraries.all():
            break
        age_bin_metadata = age_bin_metadata.iloc[
            np.flatnonzero(nonzero_age_bin_libraries)
        ].copy()
        if age_bin_metadata.empty:
            raise ValueError("All samples have zero counts over the age-bin model gene universe")

    observed_age_bins = sorted(
        age_bin_metadata["age_bin"].astype(str).unique(),
        key=lambda value: int(value.split("-", maxsplit=1)[0]),
    )
    age_bin_model_metadata = age_bin_metadata.copy()
    age_bin_model_metadata["age_bin"] = pd.Categorical(
        age_bin_model_metadata["age_bin"], categories=observed_age_bins, ordered=True,
    )
    age_bin_terms = ["age_bin"]
    if "study_site" in design_terms and _has_covariate_variation(
        age_bin_model_metadata, "study_site",
    ):
        age_bin_terms.insert(0, "study_site")
    age_bin_terms.extend(
        term for term in required_terms
        if term != "age" and _has_covariate_variation(age_bin_model_metadata, term)
    )
    age_bin_design = "~ " + " + ".join(age_bin_terms)
    for column in ("study_site", "sex", "cmv"):
        if column in age_bin_design and column in age_bin_model_metadata:
            age_bin_model_metadata[column] = pd.Categorical(
                age_bin_model_metadata[column],
                categories=categorical_levels(column, age_bin_model_metadata[column]),
            )
    age_bin_design_matrix = FormulaicContrasts(
        age_bin_model_metadata, age_bin_design,
    ).design_matrix
    age_bin_numeric_design = age_bin_design_matrix.to_numpy(dtype=float)
    age_bin_rank = int(np.linalg.matrix_rank(age_bin_numeric_design))
    if age_bin_rank < age_bin_numeric_design.shape[1] or len(age_bin_model_metadata) <= age_bin_rank:
        raise ValueError(
            f"Age-bin design {age_bin_design!r} is not estimable for "
            f"{len(age_bin_model_metadata)} samples"
        )
    age_bin_total_counts = pd.to_numeric(
        age_bin_model_metadata["total_counts"], errors="coerce",
    ).to_numpy(dtype=float)
    if not np.isfinite(age_bin_total_counts).all() or (age_bin_total_counts <= 0).any():
        raise ValueError("Age-bin normalization requires positive finite sample total_counts")
    # Estimate the size-factor reference across every cell-type sample first;
    # only then restrict to samples that pass the age-bin model filters.
    age_bin_geometric_mean = all_samples_geometric_mean_total_counts
    age_bin_size_factors = age_bin_total_counts / age_bin_geometric_mean
    age_bin_target_counts = age_bin_counts[:, age_bin_genes.get_loc(gene)]
    age_bin_normalized_counts = age_bin_target_counts / age_bin_size_factors
    age_bin_log_expression = np.log(age_bin_normalized_counts + 1.0)
    age_bin_fit = sm.OLS(age_bin_log_expression, age_bin_numeric_design).fit()
    age_bin_residuals = np.asarray(age_bin_fit.resid, dtype=float)
    fit_diagnostics_rows.append({
        "model": "Naive age-bin model",
        "n_samples": len(age_bin_log_expression),
        "residual_sd": _residual_sd(age_bin_residuals),
        "residual_sd_scale": "log(normalized count + 1)",
        "r_squared": float(age_bin_fit.rsquared),
        "r_squared_label": "In-sample R²",
        "fixed_effects_r_squared": np.nan,
    })
    age_bin_model_result = _ols_result_table(
        age_bin_fit, age_bin_design_matrix.columns.astype(str).tolist(),
    )
    age_bin_samples = age_bin_model_metadata.copy()
    age_bin_samples.insert(0, "pseudobulk_id", age_bin_samples.index.astype(str))
    age_bin_samples["raw_count"] = age_bin_target_counts
    age_bin_samples["naive_normalized_count"] = age_bin_normalized_counts
    age_bin_samples["naive_log_expression"] = age_bin_log_expression
    age_bin_samples["fitted_expression"] = np.asarray(age_bin_fit.fittedvalues, dtype=float)
    age_bin_samples["residual"] = np.asarray(age_bin_fit.resid, dtype=float)
    age_bin_predictions = []
    for age_bin in observed_age_bins:
        standardized = age_bin_model_metadata.copy()
        standardized["age_bin"] = pd.Categorical(
            [age_bin] * len(standardized), categories=observed_age_bins, ordered=True,
        )
        standardized_design = FormulaicContrasts(
            standardized, age_bin_design,
        ).design_matrix.to_numpy(dtype=float)
        average_design = standardized_design.mean(axis=0)
        prediction = age_bin_fit.t_test(average_design[None, :])
        interval = np.asarray(prediction.conf_int(), dtype=float)[0]
        age_bin_predictions.append({
            "age_bin": age_bin,
            "age_bin_start": int(age_bin.split("-", maxsplit=1)[0]),
            "age_bin_end": int(age_bin.split("-", maxsplit=1)[1]),
            "n_samples": int((age_bin_samples["age_bin"].astype(str) == age_bin).sum()),
            "fitted_expression": float(np.asarray(prediction.effect).reshape(-1)[0]),
            "ci_lower": float(interval[0]),
            "ci_upper": float(interval[1]),
        })
    age_bin_predictions = pd.DataFrame(age_bin_predictions)

    # Allow study-specific age slopes with study intercepts. Study site is
    # omitted here because it is nested within study and would be redundant.
    slope_model_metadata = model_metadata.copy()
    study_levels = categorical_levels("study", slope_model_metadata["study"])
    slope_model_metadata["study"] = pd.Categorical(
        slope_model_metadata["study"], categories=study_levels,
    )
    slope_terms = ["age * study"]
    slope_terms.extend(term for term in varying_terms if term != "age")
    slope_design = "~ " + " + ".join(slope_terms)
    slope_design_matrix = FormulaicContrasts(
        slope_model_metadata, slope_design,
    ).design_matrix
    slope_numeric_design = slope_design_matrix.to_numpy(dtype=float)
    slope_rank = int(np.linalg.matrix_rank(slope_numeric_design))
    if slope_rank < slope_numeric_design.shape[1] or len(slope_model_metadata) <= slope_rank:
        raise ValueError(
            f"Study-slope design {slope_design!r} is not estimable for "
            f"{len(slope_model_metadata)} samples"
        )
    study_slope_fit = sm.OLS(naive_log_expression, slope_numeric_design).fit()
    study_slope_residuals = np.asarray(study_slope_fit.resid, dtype=float)
    fit_diagnostics_rows.append({
        "model": "Naive age model with study-specific slopes",
        "n_samples": len(naive_log_expression),
        "residual_sd": _residual_sd(study_slope_residuals),
        "residual_sd_scale": "log(normalized count + 1)",
        "r_squared": float(study_slope_fit.rsquared),
        "r_squared_label": "In-sample R²",
        "fixed_effects_r_squared": np.nan,
    })
    study_slope_model_result = _ols_result_table(
        study_slope_fit, slope_design_matrix.columns.astype(str).tolist(),
    )
    study_slope_rows = []
    study_slope_lines = []
    for study in study_levels:
        contrast_metadata = pd.concat(
            [slope_model_metadata.iloc[[0]].copy(), slope_model_metadata.iloc[[0]].copy()],
            ignore_index=True,
        )
        contrast_metadata["age"] = [0.0, 1.0]
        contrast_metadata["study"] = pd.Categorical(
            [study, study], categories=study_levels,
        )
        contrast_design = FormulaicContrasts(
            contrast_metadata, slope_design,
        ).design_matrix.to_numpy(dtype=float)
        contrast = contrast_design[1] - contrast_design[0]
        slope_test = study_slope_fit.t_test(contrast[None, :])
        slope_interval = np.asarray(slope_test.conf_int(), dtype=float)[0]
        study_slope_rows.append({
            "study": study,
            "slope": float(np.asarray(slope_test.effect).reshape(-1)[0]),
            "standard_error": float(np.asarray(slope_test.sd).reshape(-1)[0]),
            "t_value": float(np.asarray(slope_test.tvalue).reshape(-1)[0]),
            "p_value": float(np.asarray(slope_test.pvalue).reshape(-1)[0]),
            "ci_lower": float(slope_interval[0]),
            "ci_upper": float(slope_interval[1]),
        })
        group = slope_model_metadata.loc[slope_model_metadata["study"].astype(str) == study]
        age_grid = np.linspace(float(group["age"].min()), float(group["age"].max()), 100)
        line_metadata = pd.concat(
            [slope_model_metadata.iloc[[0]].copy()] * len(age_grid), ignore_index=True,
        )
        line_metadata["age"] = age_grid
        line_metadata["study"] = pd.Categorical(
            [study] * len(age_grid), categories=study_levels,
        )
        if "log10_total_counts" in line_metadata:
            line_metadata["log10_total_counts"] = float(
                slope_model_metadata["log10_total_counts"].mean()
            )
        for column in ("sex", "cmv"):
            if column in line_metadata and column in slope_design:
                reference = categorical_levels(column, slope_model_metadata[column])[0]
                line_metadata[column] = pd.Categorical(
                    [reference] * len(age_grid),
                    categories=categorical_levels(column, slope_model_metadata[column]),
                )
        line_design = FormulaicContrasts(
            line_metadata, slope_design,
        ).design_matrix.to_numpy(dtype=float)
        study_slope_lines.append(pd.DataFrame({
            "study": study,
            "age": age_grid,
            "fitted_expression": line_design @ np.asarray(study_slope_fit.params),
        }))
    study_slope_result = pd.DataFrame(study_slope_rows)
    study_slope_predictions = pd.concat(study_slope_lines, ignore_index=True)

    # Keep study baselines as fixed effects while partially pooling age slopes
    # with a random slope deviation for each study. Centering age makes the
    # fixed study intercepts refer to a typical sample age and improves fit
    # conditioning; it does not change the per-year slope interpretation.
    mixed_model_metadata = slope_model_metadata.copy()
    mean_age = float(mixed_model_metadata["age"].mean())
    mixed_model_metadata["age_centered"] = mixed_model_metadata["age"] - mean_age
    for column in ("sex", "cmv"):
        if column in mixed_model_metadata and column in varying_terms:
            mixed_model_metadata[column] = pd.Categorical(
                mixed_model_metadata[column],
                categories=categorical_levels(column, mixed_model_metadata[column]),
            )
    mixed_covariates = [
        f"C({term})" if term in {"sex", "cmv"} else term
        for term in varying_terms if term != "age"
    ]
    mixed_formula_terms = ["age_centered", "C(study)", *mixed_covariates]
    mixed_formula = "naive_log_expression ~ " + " + ".join(mixed_formula_terms)
    mixed_model_metadata["naive_log_expression"] = naive_log_expression
    import statsmodels.formula.api as smf

    mixed_fit = None
    mixed_fit_errors = []
    mixed_fit_warnings = []
    for method in ("lbfgs", "powell"):
        try:
            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter("always")
                candidate = smf.mixedlm(
                    mixed_formula,
                    data=mixed_model_metadata,
                    groups=mixed_model_metadata["study"],
                    re_formula="0 + age_centered",
                ).fit(reml=True, method=method, disp=False)
            mixed_fit_warnings.extend(
                f"{warning.category.__name__}: {warning.message}"
                for warning in caught_warnings
            )
            mixed_fit = candidate
            if candidate.converged:
                break
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            mixed_fit_errors.append(f"{method}: {type(exc).__name__}: {exc}")
    if mixed_fit is None:
        raise RuntimeError(
            "The study random-slope mixed model could not be fit: "
            + "; ".join(mixed_fit_errors)
        )
    mixed_fixed = mixed_fit.fe_params
    mixed_confidence = mixed_fit.conf_int().loc[mixed_fixed.index]
    mixed_standard_errors = mixed_fit.bse_fe.loc[mixed_fixed.index].copy()
    mixed_statistics = mixed_fit.tvalues.loc[mixed_fixed.index].copy()
    mixed_p_values = mixed_fit.pvalues.loc[mixed_fixed.index].copy()
    mixed_uncertainty_note = "MixedLM Wald intervals"
    random_slope_variance = float(mixed_fit.cov_re.iloc[0, 0])
    if random_slope_variance <= 1e-10 or not np.isfinite(mixed_standard_errors).all():
        # At a zero-variance boundary MixedLM's Hessian can be singular even
        # though the fixed-effect estimates are identifiable. Use the
        # equivalent fixed-effects OLS covariance for those mean parameters.
        boundary_fixed_fit = smf.ols(
            mixed_formula, data=mixed_model_metadata,
        ).fit()
        mixed_standard_errors = boundary_fixed_fit.bse.loc[mixed_fixed.index]
        mixed_statistics = boundary_fixed_fit.tvalues.loc[mixed_fixed.index]
        mixed_p_values = boundary_fixed_fit.pvalues.loc[mixed_fixed.index]
        mixed_confidence = boundary_fixed_fit.conf_int().loc[mixed_fixed.index]
        mixed_uncertainty_note = (
            "Fixed-effect OLS uncertainty because the random-slope variance is at its boundary"
        )
    mixed_model_result = pd.DataFrame({
        "term": mixed_fixed.index.astype(str),
        "coefficient": mixed_fixed.to_numpy(dtype=float),
        "standard_error": mixed_standard_errors.to_numpy(dtype=float),
        "test_statistic": mixed_statistics.to_numpy(dtype=float),
        "p_value": mixed_p_values.to_numpy(dtype=float),
        "ci_lower": mixed_confidence.iloc[:, 0].to_numpy(dtype=float),
        "ci_upper": mixed_confidence.iloc[:, 1].to_numpy(dtype=float),
    })
    fixed_age_slope = float(mixed_fixed["age_centered"])
    mixed_slope_rows = []
    mixed_slope_lines = []
    random_slope_deviations = {}
    if random_slope_variance > 1e-10:
        try:
            random_slope_deviations = {
                study: float(mixed_fit.random_effects[study].iloc[0])
                for study in study_levels
            }
        except (ValueError, np.linalg.LinAlgError) as exc:
            mixed_fit_warnings.append(
                "Random effects unavailable from singular covariance; "
                f"using zero conditional deviations ({type(exc).__name__})."
            )
    for study in study_levels:
        random_slope_deviation = random_slope_deviations.get(study, 0.0)
        mixed_slope_rows.append({
            "study": study,
            "fixed_age_slope": fixed_age_slope,
            "random_slope_deviation": random_slope_deviation,
            "conditional_age_slope": fixed_age_slope + random_slope_deviation,
        })
        group = mixed_model_metadata.loc[mixed_model_metadata["study"].astype(str) == study]
        age_grid = np.linspace(float(group["age"].min()), float(group["age"].max()), 100)
        line_metadata = pd.DataFrame({
            "age_centered": age_grid - mean_age,
            "study": pd.Categorical([study] * len(age_grid), categories=study_levels),
        })
        for term in varying_terms:
            if term in {"age", "study_site"}:
                continue
            if term in {"sex", "cmv"}:
                levels = categorical_levels(term, mixed_model_metadata[term])
                line_metadata[term] = pd.Categorical([levels[0]] * len(age_grid), categories=levels)
            elif term == "log10_total_counts":
                line_metadata[term] = float(mixed_model_metadata[term].mean())
            else:
                line_metadata[term] = mixed_model_metadata[term].iloc[0]
        fixed_prediction = np.asarray(mixed_fit.predict(line_metadata), dtype=float)
        mixed_slope_lines.append(pd.DataFrame({
            "study": study,
            "age": age_grid,
            "fitted_expression": fixed_prediction + random_slope_deviation * (age_grid - mean_age),
        }))
    mixed_study_slopes = pd.DataFrame(mixed_slope_rows)
    mixed_study_predictions = pd.concat(mixed_slope_lines, ignore_index=True)
    mixed_fixed_sample_fitted = np.asarray(
        mixed_fit.predict(mixed_model_metadata), dtype=float,
    )
    mixed_sample_slope_deviation = mixed_model_metadata["study"].astype(str).map(
        random_slope_deviations,
    ).fillna(0.0).to_numpy(dtype=float)
    mixed_model_residuals = (
        naive_log_expression
        - mixed_fixed_sample_fitted
        - mixed_sample_slope_deviation * mixed_model_metadata["age_centered"].to_numpy(dtype=float)
    )
    mixed_fixed_effects_r_squared = _in_sample_r_squared(
        naive_log_expression, mixed_fixed_sample_fitted,
    )
    mixed_conditional_fitted = naive_log_expression - mixed_model_residuals
    fit_diagnostics_rows.append({
        "model": "Naive mixed age model",
        "n_samples": len(naive_log_expression),
        "residual_sd": _residual_sd(mixed_model_residuals),
        "residual_sd_scale": "log(normalized count + 1), conditional on study slopes",
        "r_squared": _in_sample_r_squared(
            naive_log_expression, mixed_conditional_fitted,
        ),
        "r_squared_label": "Conditional in-sample R²",
        "fixed_effects_r_squared": mixed_fixed_effects_r_squared,
    })
    mixed_model_diagnostics = {
        "design": mixed_formula,
        "n_samples": len(mixed_model_metadata),
        "n_studies": len(study_levels),
        "studies": study_levels,
        "mean_age_years": mean_age,
        "converged": bool(mixed_fit.converged),
        "optimizer": str(mixed_fit.method if hasattr(mixed_fit, "method") else "unknown"),
        "random_slope_sd_per_year": float(np.sqrt(max(random_slope_variance, 0.0))),
        "residual_sd": float(np.sqrt(float(mixed_fit.scale))),
        "fixed_effect_uncertainty": mixed_uncertainty_note,
        "warnings": [*mixed_fit_errors, *mixed_fit_warnings],
        "interpretation": (
            "The fixed age coefficient is the across-study mean slope; each "
            "conditional study slope adds its partially pooled random slope deviation."
        ),
    }

    counts_df = pd.DataFrame(counts, index=metadata.index, columns=fit_genes)
    dds = DeseqDataSet(
        counts=counts_df, metadata=model_metadata, design=design_matrix,
        n_cpus=cpus, quiet=True,
    )
    dds.deseq2()
    fitted_design = dds.obsm["design_matrix"]
    if "age" not in fitted_design.columns:
        raise ValueError("The configured merged model did not produce an age coefficient")
    age_contrast = np.zeros(fitted_design.shape[1], dtype=float)
    age_contrast[list(fitted_design.columns).index("age")] = 1.0
    age_stats = DeseqStats(
        dds, contrast=age_contrast, alpha=float(spec["alpha"]),
        n_cpus=cpus, quiet=True,
    )
    age_stats.summary()
    gene_result = age_stats.results_df.loc[[gene]].reset_index(names="gene")

    # Capture final-model diagnostics before adding the VST layer for downstream use.
    variance = dds.var.copy()
    fitted_means = np.asarray(dds.obsm["_mu_LFC"], dtype=float)
    gene_dispersion = float(variance.loc[gene, "dispersions"])
    gene_mean = fitted_means[:, gene_position]
    gene_variance = gene_mean + gene_dispersion * np.square(gene_mean)
    if not np.isfinite(gene_variance).all() or (gene_variance <= 0).any():
        raise ValueError(f"The fitted model did not produce valid Pearson residuals for {gene!r}")
    gene_pearson_residuals = (
        counts[:, gene_position] - gene_mean
    ) / np.sqrt(gene_variance)
    fit_diagnostics_rows.append({
        "model": "PyDESeq2",
        "n_samples": len(gene_pearson_residuals),
        "residual_sd": _residual_sd(gene_pearson_residuals),
        "residual_sd_scale": "Pearson residuals",
        "r_squared": _in_sample_r_squared(counts[:, gene_position], gene_mean),
        "r_squared_label": "Raw-count in-sample R²",
        "fixed_effects_r_squared": np.nan,
    })
    all_dispersions = pd.to_numeric(
        variance["dispersions"], errors="coerce",
    ).to_numpy(dtype=float)
    all_gene_variances = fitted_means + all_dispersions[None, :] * np.square(fitted_means)
    valid_residual_genes = (
        np.isfinite(all_dispersions) & (all_dispersions >= 0)
        & np.isfinite(fitted_means).all(axis=0)
        & np.isfinite(all_gene_variances).all(axis=0)
        & (all_gene_variances > 0).all(axis=0)
    )
    residual_gene_names = fit_genes[valid_residual_genes]
    residual_matrix = (
        counts[:, valid_residual_genes] - fitted_means[:, valid_residual_genes]
    ) / np.sqrt(all_gene_variances[:, valid_residual_genes])
    target_residual_position = residual_gene_names.get_loc(gene)
    target_residuals = residual_matrix[:, target_residual_position]
    centered_residuals = residual_matrix - residual_matrix.mean(axis=0)
    centered_target = target_residuals - target_residuals.mean()
    denominators = np.sqrt(
        np.sum(np.square(centered_target))
        * np.sum(np.square(centered_residuals), axis=0)
    )
    correlations = np.full(len(residual_gene_names), np.nan, dtype=float)
    np.divide(
        centered_target @ centered_residuals,
        denominators,
        out=correlations,
        where=denominators > 0,
    )
    correlations[target_residual_position] = np.nan
    correlation_table = pd.DataFrame({
        "gene": residual_gene_names.astype(str),
        "pearson_r": correlations,
    }).dropna(subset=["pearson_r"])
    correlation_table["absolute_pearson_r"] = correlation_table["pearson_r"].abs()
    correlation_table = correlation_table.sort_values(
        "absolute_pearson_r", ascending=False,
    ).reset_index(drop=True)
    top_correlated = correlation_table.head(4)
    top_residual_tables = []
    for correlated in top_correlated.itertuples(index=False):
        correlated_position = residual_gene_names.get_loc(correlated.gene)
        sample_residuals = model_metadata.copy()
        sample_residuals.insert(0, "pseudobulk_id", sample_residuals.index.astype(str))
        sample_residuals["target_gene_pearson_residual"] = target_residuals
        sample_residuals["correlated_gene"] = correlated.gene
        sample_residuals["correlated_gene_pearson_residual"] = (
            residual_matrix[:, correlated_position]
        )
        sample_residuals["pearson_r"] = correlated.pearson_r
        top_residual_tables.append(sample_residuals)
    top_correlated_residuals = (
        pd.concat(top_residual_tables, ignore_index=True)
        if top_residual_tables else pd.DataFrame()
    )
    dds.vst()
    normalized = np.asarray(dds.layers["normed_counts"])
    vst_counts = np.asarray(dds.layers["vst_counts"])
    gene_samples = model_metadata.copy()
    gene_samples.insert(0, "pseudobulk_id", gene_samples.index.astype(str))
    gene_samples["raw_count"] = counts[:, gene_position]
    gene_samples["fitted_mean"] = gene_mean
    gene_samples["pearson_residual"] = gene_pearson_residuals
    gene_samples["normalized_count"] = normalized[:, gene_position]
    gene_samples["vst_count"] = vst_counts[:, gene_position]
    gene_samples["naive_size_factor"] = naive_size_factors
    gene_samples["naive_normalized_count"] = naive_normalized_count
    gene_samples["naive_log_expression"] = naive_log_expression
    gene_samples["naive_fitted"] = np.asarray(naive_fit.fittedvalues, dtype=float)
    gene_samples["naive_residual"] = np.asarray(naive_fit.resid, dtype=float)
    gene_samples["study_slope_residual"] = study_slope_residuals
    gene_samples["mixed_model_residual"] = mixed_model_residuals

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_samples.to_csv(output_dir / "gene_sample_counts.csv", index=False)
    gene_result.to_csv(output_dir / "gene_age_result.csv", index=False)
    naive_model_result.to_csv(output_dir / "naive_linear_model_result.csv", index=False)
    age_bin_samples.to_csv(output_dir / "naive_age_bin_samples.csv", index=False)
    age_bin_model_result.to_csv(output_dir / "naive_age_bin_model_result.csv", index=False)
    age_bin_predictions.to_csv(output_dir / "naive_age_bin_predictions.csv", index=False)
    study_slope_model_result.to_csv(
        output_dir / "naive_study_slope_model_result.csv", index=False,
    )
    study_slope_result.to_csv(output_dir / "naive_study_slopes.csv", index=False)
    study_slope_predictions.to_csv(
        output_dir / "naive_study_slope_predictions.csv", index=False,
    )
    fit_diagnostics = pd.DataFrame(fit_diagnostics_rows)
    fit_diagnostics.to_csv(output_dir / "model_fit_diagnostics.csv", index=False)
    mixed_model_result.to_csv(output_dir / "naive_mixed_model_fixed_effects.csv", index=False)
    mixed_study_slopes.to_csv(output_dir / "naive_mixed_model_study_slopes.csv", index=False)
    mixed_study_predictions.to_csv(
        output_dir / "naive_mixed_model_predictions.csv", index=False,
    )
    correlation_table.to_csv(output_dir / "gene_residual_correlations.csv", index=False)
    top_correlated_residuals.to_csv(
        output_dir / "top_correlated_gene_residuals.csv", index=False,
    )
    full_vst = pd.DataFrame(vst_counts, index=metadata.index.astype(str), columns=fit_genes)
    full_vst.index.name = "pseudobulk_id"
    full_vst.to_csv(output_dir / "vst_normalized_counts.csv.gz", compression="gzip")
    model_metadata.assign(pseudobulk_id=model_metadata.index.astype(str)).to_csv(
        output_dir / "fit_sample_metadata.csv", index=False,
    )

    variance_table = variance.copy()
    if "_normed_means" not in variance_table:
        raise ValueError("PyDESeq2 did not return normalized means for dispersion diagnostics")
    variance_table["mean_normalized_count"] = variance_table["_normed_means"]
    variance_table["gene"] = variance_table.index.astype(str)
    variance_table["is_target_gene"] = variance_table["gene"] == gene
    variance_table.to_csv(output_dir / "mean_variance.csv", index=False)

    top_gene_correlations = correlation_table.head(10).merge(
        variance_table[["gene", "mean_normalized_count", "dispersions"]].reset_index(
            drop=True,
        ),
        # `variance_table` keeps a gene-named index as well as a gene column.
        # Drop the index name so pandas treats the merge key unambiguously.
        on="gene",
        validate="one_to_one",
    )
    target_variance = variance_table.loc[variance_table["gene"] == gene].iloc[0]
    target_correlation_row = pd.DataFrame([{
        "gene": gene,
        "pearson_r": 1.0,
        "absolute_pearson_r": 1.0,
        "mean_normalized_count": target_variance["mean_normalized_count"],
        "dispersions": target_variance["dispersions"],
    }])
    top_gene_correlations = pd.concat(
        [target_correlation_row, top_gene_correlations], ignore_index=True,
    )
    top_gene_correlations.to_csv(
        output_dir / "top_gene_residual_correlations.csv", index=False,
    )

    result = gene_result.iloc[0].to_dict()
    summary = {
        "cell_type": cell_type,
        "gene": gene,
        "input": str(input_path),
        "design": design,
        "geometric_mean_total_counts": geometric_mean_total_counts,
        "age_bin_model": {
            "design": age_bin_design,
            "bin_width_years": bin_width,
            "bins": observed_age_bins,
            "reference_bin": observed_age_bins[0],
            "n_samples": len(age_bin_model_metadata),
            "geometric_mean_total_counts": age_bin_geometric_mean,
            "normalization_n_samples": int(valid_all_type_totals.sum()),
            "metadata_exclusions": age_bin_exclusions,
        },
        "study_slope_model": {
            "design": slope_design,
            "studies": study_levels,
            "n_samples": len(slope_model_metadata),
        },
        "mixed_model": mixed_model_diagnostics,
        "n_samples": len(model_metadata),
        "n_studies": len(studies),
        "studies": studies,
        "n_genes_fit": len(fit_genes),
        "metadata_exclusions": exclusions,
        "age_result": {
            key: (None if pd.isna(value) else float(value) if isinstance(value, np.number) else str(value))
            for key, value in result.items()
        },
        "outputs": [
            "gene_sample_counts.csv", "gene_age_result.csv", "mean_variance.csv",
            "model_fit_diagnostics.csv",
            "naive_linear_model_result.csv", "naive_age_bin_samples.csv",
            "naive_age_bin_model_result.csv", "naive_age_bin_predictions.csv",
            "naive_study_slope_model_result.csv", "naive_study_slopes.csv",
            "naive_study_slope_predictions.csv",
            "naive_mixed_model_fixed_effects.csv", "naive_mixed_model_study_slopes.csv",
            "naive_mixed_model_predictions.csv",
            "gene_residual_correlations.csv", "top_correlated_gene_residuals.csv",
            "top_gene_residual_correlations.csv",
            "fit_sample_metadata.csv", "vst_normalized_counts.csv.gz",
        ],
    }
    (output_dir / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    template_path = Path(__file__).resolve().parents[1] / "reports" / "gene_profile_report.py"
    project_root = Path(__file__).resolve().parents[1]
    return generate_gene_profile_report(
        profile_dir=output_dir,
        template_path=template_path,
        project_root=project_root,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Refit the merged DESeq2 model for one cell type and profile one gene"
    )
    parser.add_argument(
        "--cell-type", required=True, help="Configured cell-type label or unique slug"
    )
    parser.add_argument("--gene", required=True, help="Exact gene identifier in the pseudobulk object")
    parser.add_argument("--input", type=Path, default=Path("output/merged/pseudobulk_merged.h5ad"))
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/gene_profiles"))
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()
    report = _fit_profile(
        input_path=args.input, config_path=args.config, cell_type=args.cell_type,
        gene=args.gene, output_dir=args.output_dir, cpus=args.cpus,
    )
    print(report)


if __name__ == "__main__":
    main()

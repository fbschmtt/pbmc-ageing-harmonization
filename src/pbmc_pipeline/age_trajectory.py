"""Age-bin support filtering and trajectory-specific DE result calculations."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .config import AgeTrajectorySettings

MANUALLY_EXCLUDED_AGE_BINS = ("90-100",)


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Adjust finite p-values with Benjamini-Hochberg, preserving missing entries."""
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(pvalues))
    if not len(finite_indices):
        return adjusted
    order = finite_indices[np.argsort(pvalues[finite_indices])]
    ranked = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def joint_wald_test(
    design: np.ndarray,
    dispersions: np.ndarray,
    coefficients: np.ndarray,
    means: np.ndarray,
    coefficient_indices: list[int],
    *,
    chunk_size: int = 512,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Test a group of GLM coefficients jointly using PyDESeq2's Wald covariance."""
    n_genes = coefficients.shape[0]
    statistics = np.full(n_genes, np.nan, dtype=float)
    pvalues = np.full(n_genes, np.nan, dtype=float)
    n_parameters = design.shape[1]
    ridge = np.eye(n_parameters) * 1e-6
    degrees_of_freedom = len(coefficient_indices)
    if degrees_of_freedom == 0:
        return statistics, pvalues, 0

    for start in range(0, n_genes, chunk_size):
        stop = min(start + chunk_size, n_genes)
        dispersion = dispersions[start:stop]
        gene_means = means[:, start:stop]
        gene_coefficients = coefficients[start:stop]
        valid = (
            np.isfinite(dispersion)
            & (dispersion >= 0)
            & np.isfinite(gene_means).all(axis=0)
            & np.isfinite(gene_coefficients).all(axis=1)
        )
        if not valid.any():
            continue
        local_means = gene_means[:, valid]
        weights = local_means / (1 + local_means * dispersion[valid][None, :])
        information = np.einsum(
            "ng,ni,nj->gij", weights, design, design, optimize=True
        )
        inverse = np.linalg.inv(information + ridge)
        covariance = inverse @ information @ inverse
        selected_covariance = covariance[:, coefficient_indices, :][:, :, coefficient_indices]
        selected_inverse = np.linalg.inv(selected_covariance)
        selected_coefficients = gene_coefficients[valid][:, coefficient_indices]
        local_statistics = np.einsum(
            "gi,gij,gj->g", selected_coefficients, selected_inverse, selected_coefficients
        )
        target = np.arange(start, stop)[valid]
        statistics[target] = local_statistics
        pvalues[target] = chi2.sf(local_statistics, degrees_of_freedom)
    return statistics, pvalues, degrees_of_freedom


def prepare_age_trajectory_metadata(
    metadata: pd.DataFrame,
    settings: AgeTrajectorySettings,
) -> tuple[pd.DataFrame | None, dict]:
    """Apply the strict age ceiling, then retain sufficiently supported bins."""
    age = pd.to_numeric(metadata["age"], errors="coerce")
    raw_starts = (
        np.floor(age / settings.bin_width_years) * settings.bin_width_years
    ).astype("Int64")
    raw_age_bins = raw_starts.map(
        lambda start: f"{int(start)}-{int(start) + settings.bin_width_years}"
        if not pd.isna(start) else pd.NA
    ).astype("string")
    raw_observed_starts = raw_starts.dropna()
    if raw_observed_starts.empty:
        prefilter_bin_counts = {}
    else:
        raw_maximum_start = int(raw_observed_starts.max())
        raw_candidate_bins = [
            f"{start}-{start + settings.bin_width_years}"
            for start in range(
                settings.reference_bin_start_age,
                raw_maximum_start + 1,
                settings.bin_width_years,
            )
        ]
        raw_counts = raw_age_bins.value_counts().to_dict()
        prefilter_bin_counts = {
            age_bin: int(raw_counts.get(age_bin, 0)) for age_bin in raw_candidate_bins
        }
    manually_excluded_count = int(raw_age_bins.isin(MANUALLY_EXCLUDED_AGE_BINS).sum())
    age_cutoff_keep = age < settings.strict_age_cutoff_exclusive
    excluded_by_age_cutoff = int((age.notna() & ~age_cutoff_keep).sum())
    keep = age_cutoff_keep & ~raw_age_bins.isin(MANUALLY_EXCLUDED_AGE_BINS)
    metadata = metadata.loc[keep].copy()
    age = age.loc[keep]
    starts = (np.floor(age / settings.bin_width_years) * settings.bin_width_years).astype("Int64")
    metadata["age_bin"] = starts.map(
        lambda start: f"{int(start)}-{int(start) + settings.bin_width_years}"
        if not pd.isna(start) else pd.NA
    ).astype("string")
    reference_bin = (
        f"{settings.reference_bin_start_age}-"
        f"{settings.reference_bin_start_age + settings.bin_width_years}"
    )
    observed_starts = starts.dropna()
    if observed_starts.empty:
        return None, {
            "status": "skipped",
            "reason": "no samples have a valid age for decade binning",
            "reference_bin": reference_bin,
            "strict_age_cutoff_exclusive": settings.strict_age_cutoff_exclusive,
            "samples_excluded_by_age_cutoff": excluded_by_age_cutoff,
            "n_samples_after_age_cutoff": len(metadata),
            "minimum_samples_per_bin": settings.minimum_samples_per_bin,
            "manually_excluded_age_bins": list(MANUALLY_EXCLUDED_AGE_BINS),
            "samples_in_manually_excluded_bins": manually_excluded_count,
            "age_bin_counts_before_filtering": prefilter_bin_counts,
            "age_bin_counts": {},
            "retained_bins": [],
            "dropped_bins": [],
        }

    maximum_start = int(observed_starts.max())
    candidate_bins = [
        f"{start}-{start + settings.bin_width_years}"
        for start in range(
            settings.reference_bin_start_age, maximum_start + 1, settings.bin_width_years
        )
    ]
    counts = metadata["age_bin"].value_counts().to_dict()
    age_bin_counts = {age_bin: int(counts.get(age_bin, 0)) for age_bin in candidate_bins}
    retained_bins = [
        age_bin for age_bin in candidate_bins
        if age_bin_counts[age_bin] >= settings.minimum_samples_per_bin
    ]
    dropped_bins = [age_bin for age_bin in candidate_bins if age_bin not in retained_bins]
    status = {
        "reference_bin": reference_bin,
        "strict_age_cutoff_exclusive": settings.strict_age_cutoff_exclusive,
        "samples_excluded_by_age_cutoff": excluded_by_age_cutoff,
        "n_samples_after_age_cutoff": len(metadata),
        "manually_excluded_age_bins": list(MANUALLY_EXCLUDED_AGE_BINS),
        "samples_in_manually_excluded_bins": manually_excluded_count,
        "age_bin_counts_before_filtering": prefilter_bin_counts,
        "bin_width_years": settings.bin_width_years,
        "minimum_samples_per_bin": settings.minimum_samples_per_bin,
        "minimum_bins": settings.minimum_bins,
        "age_bin_counts": age_bin_counts,
        "retained_bins": retained_bins,
        "dropped_bins": dropped_bins,
        "samples_in_dropped_bins": int(metadata["age_bin"].isin(dropped_bins).sum()),
        "n_samples_retained": int(metadata["age_bin"].isin(retained_bins).sum()),
        "reference_samples_by_study": {
            str(study): int(count)
            for study, count in metadata.loc[metadata["age_bin"] == reference_bin]
            .groupby("study", observed=True).size().items()
        },
        "samples_by_study_and_bin": {
            f"{study}|{age_bin}": int(count)
            for (study, age_bin), count in metadata.loc[
                metadata["age_bin"].isin(retained_bins)
            ].groupby(["study", "age_bin"], observed=True).size().items()
        },
    }
    if reference_bin not in retained_bins:
        return None, {
            **status,
            "status": "skipped",
            "reason": (
                f"reference age bin {reference_bin!r} has fewer than "
                f"{settings.minimum_samples_per_bin} eligible samples"
            ),
        }
    if len(retained_bins) < settings.minimum_bins:
        return None, {
            **status,
            "status": "skipped",
            "reason": (
                f"only {len(retained_bins)} age bins meet the minimum sample count; "
                f"at least {settings.minimum_bins} retained bins are required"
            ),
        }
    selected = metadata.loc[metadata["age_bin"].isin(retained_bins)].copy()
    selected["age_bin"] = pd.Categorical(
        selected["age_bin"], categories=retained_bins, ordered=True
    )
    return selected, {**status, "status": "ready"}


def age_trajectory_design(spec: dict) -> str:
    """Replace the configured continuous-age term with an ordered age-bin factor."""
    age_term = spec["age_term"]
    terms = [term.strip() for term in spec["merged_design"].removeprefix("~").split("+")]
    if age_term not in terms:
        raise ValueError(f"The merged DE design does not contain the configured age term {age_term!r}")
    return "~ " + " + ".join("age_bin" if term == age_term else term for term in terms)


def attach_trajectory_statistics(
    combined: pd.DataFrame,
    dds,
    design_matrix: pd.DataFrame,
    levels: list[str],
) -> pd.DataFrame:
    """Add joint age-bin inference, reference rows, and standardized profiles."""
    design_columns = list(design_matrix.columns)
    bin_terms = [term for term in design_columns if term.startswith("age_bin[T.")]
    bin_indices = [design_columns.index(term) for term in bin_terms]
    coefficients = dds.varm["LFC"].loc[:, design_columns].to_numpy(dtype=float)
    numeric_design = design_matrix.to_numpy(dtype=float)
    log_means = (
        numeric_design @ coefficients.T
        + np.log(dds.obs["size_factors"].to_numpy(dtype=float))[:, None]
    )
    means = np.exp(np.clip(log_means, -700, 700))
    omnibus_statistic, omnibus_pvalue, omnibus_df = joint_wald_test(
        numeric_design,
        dds.var["dispersions"].to_numpy(dtype=float),
        coefficients,
        means,
        bin_indices,
    )
    omnibus = pd.DataFrame({
        "gene": dds.var_names.astype(str),
        "omnibus_statistic": omnibus_statistic,
        "omnibus_df": omnibus_df,
        "omnibus_pvalue": omnibus_pvalue,
        "omnibus_padj": benjamini_hochberg(omnibus_pvalue),
    })
    combined = combined.merge(omnibus, on="gene", how="left", validate="many_to_one")

    reference_bin = levels[0]
    reference_rows = pd.DataFrame({
        "gene": dds.var_names.astype(str),
        "baseMean": dds.var["_normed_means"].to_numpy(dtype=float),
        "log2FoldChange": 0.0,
        "lfcSE": np.nan,
        "stat": np.nan,
        "pvalue": np.nan,
        "padj": np.nan,
        "covariate": "age_bin",
        "contrast": f"{reference_bin} reference",
        "age_bin": reference_bin,
    }).merge(omnibus, on="gene", how="left", validate="one_to_one")
    combined = pd.concat([reference_rows, combined], ignore_index=True)

    trajectory = combined.pivot(
        index="gene", columns="age_bin", values="log2FoldChange"
    ).reindex(columns=levels)
    complete = trajectory.notna().all(axis=1)
    trajectory_sd = trajectory.std(axis=1, ddof=0)
    scalable = complete & np.isfinite(trajectory_sd) & (trajectory_sd > 0)
    scaled = trajectory.div(trajectory_sd.where(scalable), axis=0)
    scaled_long = (
        scaled.rename_axis(columns="age_bin")
        .stack(future_stack=True)
        .rename("trajectory_scaled")
        .reset_index()
    )
    combined = combined.merge(
        scaled_long, on=["gene", "age_bin"], how="left", validate="one_to_one"
    )
    combined["trajectory_sd"] = combined["gene"].map(trajectory_sd.to_dict())
    combined["trajectory_scale_bins"] = len(levels)
    combined["age_bin"] = pd.Categorical(combined["age_bin"], categories=levels, ordered=True)
    return combined.sort_values(["gene", "age_bin"], kind="stable")

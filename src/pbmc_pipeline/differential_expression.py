"""Run sample-level covariate differential expression on pseudobulk counts."""
from __future__ import annotations

import argparse
import json
import logging
import re
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .age_trajectory import (
    age_trajectory_design,
    attach_trajectory_statistics,
    prepare_age_trajectory_metadata,
)
from .config import AgeTrajectorySettings, config_digest, read_json
from .covariates import categorical_levels, normalize_cmv_status

LOGGER = logging.getLogger(__name__)
_UNKNOWN = {"", "not_provided", "heterogeneous during bulk", "nan", "none", "unknown"}
_TEST_MODE_WARNING = (
    "TEST MODE: the configured minimum cells per pseudobulk was bypassed so the "
    "small fixture can exercise model fitting. These results are for software "
    "validation only and must not be interpreted as biological evidence."
)
_SYNTHETIC_DATA_WARNING = (
    "SYNTHETIC TEST DATA: pseudobulk counts and sample metadata were generated "
    "deterministically to exercise the DE workflow. These results are not biological evidence."
)


class _NonEstimableModel(ValueError):
    """Expected design or sample-size condition that prevents a DE fit."""


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "unnamed"


def _matrix_to_counts(matrix, *, context: str) -> np.ndarray:
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{context}: counts must be finite and non-negative")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError(f"{context}: PyDESeq2 requires integer raw counts")
    return values.astype(np.int64, copy=False)


def _prepare_metadata(
    obs: pd.DataFrame, inclusion: dict, *, test_mode: bool = False,
) -> tuple[pd.DataFrame, dict[str, int], list[dict]]:
    optional_covariates = [column for column in ("bmi", "cmv") if column in obs]
    metadata = obs[[
        "study", "study_site", "sample", "age", "sex", "n_cells", "total_counts",
        *optional_covariates,
    ]].copy()
    metadata["study"] = metadata["study"].astype("string").str.strip()
    metadata["study_site"] = metadata["study_site"].astype("string").str.strip()
    metadata["sample"] = metadata["sample"].astype("string").str.strip()
    metadata["sex"] = metadata["sex"].astype("string").str.strip().str.lower()
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce")
    metadata["n_cells"] = pd.to_numeric(metadata["n_cells"], errors="coerce")
    metadata["total_counts"] = pd.to_numeric(metadata["total_counts"], errors="coerce")
    positive_total_counts = metadata["total_counts"].where(metadata["total_counts"] > 0)
    metadata["log10_total_counts"] = np.log10(positive_total_counts)
    if "bmi" in metadata:
        metadata["bmi"] = pd.to_numeric(metadata["bmi"], errors="coerce")
    if "cmv" in metadata:
        metadata["cmv"] = normalize_cmv_status(metadata["cmv"])
    minimum_age = inclusion["minimum_age_years_inclusive"]
    minimum_cells = inclusion["minimum_cells_per_pseudobulk_inclusive"]
    valid_study = metadata["study"].notna() & ~metadata["study"].str.lower().isin(_UNKNOWN)
    valid_study_site = (
        metadata["study_site"].notna()
        & ~metadata["study_site"].str.lower().isin(_UNKNOWN)
    )
    valid_sample = metadata["sample"].notna() & ~metadata["sample"].str.lower().isin(_UNKNOWN)
    valid_age = metadata["age"].notna() & np.isfinite(metadata["age"])
    adult = valid_age & (metadata["age"] >= minimum_age)
    valid_sex = metadata["sex"].notna() & ~metadata["sex"].str.lower().isin(_UNKNOWN)
    valid_cells = metadata["n_cells"].notna() & np.isfinite(metadata["n_cells"])
    valid_total_counts = (
        metadata["total_counts"].notna()
        & np.isfinite(metadata["total_counts"])
        & (metadata["total_counts"] > 0)
    )
    enough_cells = valid_cells & (metadata["n_cells"] >= minimum_cells)
    valid_other_metadata = (
        valid_study & valid_study_site & valid_sample & adult & valid_sex & valid_cells
        & valid_total_counts
    )
    valid = valid_other_metadata & (enough_cells | test_mode)
    excluded = []
    for position in np.flatnonzero(~valid.to_numpy(dtype=bool)):
        row = metadata.iloc[position]
        reasons = []
        if not valid_study.iloc[position]:
            reasons.append("missing_or_unknown_study")
        if not valid_study_site.iloc[position]:
            reasons.append("missing_or_unknown_study_site")
        if not valid_sample.iloc[position]:
            reasons.append("missing_or_unknown_sample")
        if not valid_age.iloc[position]:
            reasons.append("missing_or_non_numeric_age")
        elif not adult.iloc[position]:
            reasons.append("age_below_minimum")
        if not valid_sex.iloc[position]:
            reasons.append("missing_or_unknown_sex")
        if not valid_cells.iloc[position]:
            reasons.append("missing_or_non_numeric_n_cells")
        elif not enough_cells.iloc[position] and not test_mode:
            reasons.append("n_cells_below_minimum")
        if not valid_total_counts.iloc[position]:
            reasons.append("missing_or_nonpositive_total_counts")
        excluded.append(
            {
                "pseudobulk_id": str(metadata.index[position]),
                "study": None if pd.isna(row["study"]) else str(row["study"]),
                "sample": None if pd.isna(row["sample"]) else str(row["sample"]),
                "reasons": reasons,
            }
        )
    exclusion_counts = {
        "missing_or_unknown_study": int((~valid_study).sum()),
        "missing_or_unknown_study_site": int((~valid_study_site).sum()),
        "missing_or_unknown_sample": int((~valid_sample).sum()),
        "missing_or_non_numeric_age": int((~valid_age).sum()),
        "age_below_minimum": int((valid_age & ~adult).sum()),
        "missing_or_unknown_sex": int((~valid_sex).sum()),
        "missing_or_non_numeric_n_cells": int((~valid_cells).sum()),
        "missing_or_nonpositive_total_counts": int((~valid_total_counts).sum()),
        "n_cells_below_minimum": int((valid_cells & ~enough_cells).sum()),
        "n_cells_filter_bypassed": int(
            (valid_other_metadata & ~enough_cells).sum() if test_mode else 0
        ),
        "eligible_before_duplicate_check": int(valid.sum()),
    }
    metadata = metadata.loc[valid].copy()
    duplicate = metadata.duplicated(["study", "sample"], keep=False)
    if duplicate.any():
        keys = metadata.loc[duplicate, ["study", "sample"]].astype(str).to_dict("records")
        raise ValueError(f"Pseudobulk must have one row per study/sample: {keys[:5]}")
    return metadata, exclusion_counts, excluded


def _availability_mask(var: pd.DataFrame, studies: list[str], *, context: str) -> np.ndarray:
    columns = [f"available_in_{study}" for study in studies]
    missing = sorted(set(columns) - set(var.columns))
    if missing:
        raise ValueError(f"{context}: gene-availability flags are missing: {missing}")
    availability = var[columns].fillna(False).astype(bool)
    return availability.all(axis=1).to_numpy()


def _has_covariate_variation(metadata: pd.DataFrame, covariate: str) -> bool:
    values = metadata[covariate]
    if covariate in {"age", "bmi", "log10_total_counts"}:
        values = pd.to_numeric(values, errors="coerce")
    return values.nunique() >= 2


def _maximal_per_study_covariates(metadata: pd.DataFrame, optional_covariates: list[str]) -> list[str]:
    """Select every varying covariate for one study's single DE design."""
    candidates = ["age", "sex", "log10_total_counts", *optional_covariates]
    return [
        covariate for covariate in candidates
        if covariate in metadata and _has_covariate_variation(metadata, covariate)
    ]


def _design_with_covariates(base_design: str, covariates: list[str]) -> str:
    """Append selected optional covariates once, preserving configured term order."""
    configured_terms = [
        term.strip() for term in base_design.removeprefix("~").split("+") if term.strip()
    ]
    optional_terms = [term for term in covariates if term not in configured_terms]
    return "~ " + " + ".join([*configured_terms, *optional_terms])


def _fit_covariate_model(
    counts: np.ndarray,
    genes: pd.Index,
    metadata: pd.DataFrame,
    *,
    design: str,
    covariates: list[str],
    alpha: float,
    cpus: int,
) -> dict[str, pd.DataFrame]:
    from formulaic_contrasts import FormulaicContrasts
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_df = pd.DataFrame(counts, index=metadata.index, columns=genes)
    model_metadata = metadata.copy()
    for column in ("study_site", "sex", "cmv", "age_bin"):
        if column in design and column in model_metadata:
            model_metadata[column] = pd.Categorical(
                model_metadata[column],
                categories=categorical_levels(column, model_metadata[column]),
            )

    design_matrix = FormulaicContrasts(model_metadata, design).design_matrix
    numeric_design = design_matrix.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(numeric_design))
    if rank < numeric_design.shape[1]:
        raise _NonEstimableModel("design matrix is rank deficient for the available samples")
    if numeric_design.shape[0] <= rank:
        raise _NonEstimableModel(
            f"{numeric_design.shape[0]} samples for {rank} design coefficients; "
            "positive residual degrees of freedom are required"
        )

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=model_metadata,
        design=design_matrix,
        n_cpus=cpus,
        quiet=True,
    )
    dds.deseq2()
    design_matrix = dds.obsm["design_matrix"]
    coefficient_results = {}
    for covariate in covariates:
        if covariate in {"age", "bmi", "log10_total_counts"}:
            terms = [covariate] if covariate in design_matrix.columns else []
            levels = []
        else:
            terms = [
                term for term in design_matrix.columns
                if term.startswith(f"{covariate}[T.")
            ]
            levels = categorical_levels(covariate, model_metadata[covariate])
        if not terms:
            raise _NonEstimableModel(
                f"no estimable coefficient was found for {covariate!r} in "
                f"design columns {list(design_matrix.columns)}"
            )
        pieces = []
        for term in terms:
            contrast = np.zeros(design_matrix.shape[1], dtype=float)
            contrast[list(design_matrix.columns).index(term)] = 1.0
            stats = DeseqStats(
                dds,
                contrast=contrast,
                alpha=alpha,
                n_cpus=cpus,
                quiet=True,
            )
            stats.summary()
            result = stats.results_df.copy()
            if covariate in {"age", "bmi", "log10_total_counts"}:
                contrast_label = {
                    "age": "per year",
                    "bmi": "per BMI unit",
                    "log10_total_counts": "per log10(total_counts) unit",
                }[covariate]
            else:
                level = term.removeprefix(f"{covariate}[T.").removesuffix("]")
                contrast_label = f"{level} vs {levels[0]}"
            result["covariate"] = covariate
            result["contrast"] = contrast_label
            if covariate == "age_bin":
                result["age_bin"] = level
            pieces.append(result)
        combined = pd.concat(pieces).rename_axis("gene").reset_index()
        if covariate == "age_bin":
            combined = attach_trajectory_statistics(combined, dds, design_matrix, levels)
        coefficient_results[covariate] = combined
    return coefficient_results


def _run_model(
    *,
    adata,
    obs: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_type: str,
    model_name: str,
    study: str | None,
    design: str,
    covariates: list[str],
    alpha: float,
    cpus: int,
    output_dir: Path,
    test_mode: bool,
    synthetic_test_data: bool,
) -> dict:
    details: dict = {
        "model": model_name,
        "study": study,
        "design": design,
        "covariates": covariates,
        "n_samples": len(metadata),
        "status": "skipped",
    }
    if metadata.empty:
        details["reason"] = "no samples pass the configured inclusion filters"
        return details

    design_covariates = [
        column.strip()
        for column in design.removeprefix("~").split("+")
        if column.strip() and column.strip() != "study_site"
    ]
    has_site_term = "study_site" in design
    n_samples_before_complete_case_filter = len(metadata)
    valid = pd.Series(True, index=metadata.index)
    missing_by_covariate = {}
    if has_site_term:
        site_values = metadata["study_site"].astype("string").str.strip().str.lower()
        site_observed = site_values.notna() & ~site_values.isin(_UNKNOWN)
        missing_by_covariate["study_site"] = int((~site_observed).sum())
        valid &= site_observed
    for column in design_covariates:
        if column in {"age", "bmi", "log10_total_counts"}:
            values = pd.to_numeric(metadata[column], errors="coerce")
            observed = values.notna() & np.isfinite(values)
        else:
            values = metadata[column].astype("string").str.strip().str.lower()
            observed = values.notna() & ~values.isin(_UNKNOWN)
        missing_by_covariate[column] = int((~observed).sum())
        valid &= observed
    metadata = metadata.loc[valid].copy()
    details.update({
        "n_samples_before_complete_case_filter": n_samples_before_complete_case_filter,
        "n_samples_excluded_missing_design_covariates": int((~valid).sum()),
        "missing_values_by_design_covariate": missing_by_covariate,
    })
    for covariate in covariates:
        if not _has_covariate_variation(metadata, covariate):
            details["reason"] = f"{covariate} has fewer than two observed values"
            return details
    if metadata.empty:
        details["reason"] = "no samples have complete values for this model"
        return details
    if has_site_term and metadata["study_site"].nunique() < 2:
        # A single-site fit remains estimable without a redundant batch term.
        has_site_term = False
    design_covariates = [
        covariate for covariate in design_covariates
        if covariate in covariates or _has_covariate_variation(metadata, covariate)
    ]
    design_terms = (["study_site"] if has_site_term else []) + design_covariates
    design = "~ " + " + ".join(design_terms)
    details["design"] = design

    zero_library_samples = []
    while True:
        studies = sorted(metadata["study"].astype(str).unique())
        study_sample_counts = {
            str(key): int(value)
            for key, value in metadata.groupby("study", observed=True).size().items()
        }
        study_site_sample_counts = {
            str(key): int(value)
            for key, value in metadata.groupby("study_site", observed=True).size().items()
        }
        details.update({
            "studies": studies,
            "study_sample_counts": study_sample_counts,
            "study_sites": sorted(metadata["study_site"].astype(str).unique()),
            "study_site_sample_counts": study_site_sample_counts,
            "n_samples": len(metadata),
        })
        var_mask = _availability_mask(adata.var, studies, context=f"{cell_type}/{model_name}")
        genes = adata.var_names[var_mask]
        if len(genes) == 0:
            details["reason"] = "no genes are available in every study included in this fit"
            return details

        # The order of pseudobulk rows in the AnnData object matches obs/metadata.
        row_positions = obs.index.get_indexer(metadata.index)
        col_positions = adata.var_names.get_indexer(genes)
        count_values = _matrix_to_counts(
            adata.X[row_positions, :][:, col_positions], context=f"{cell_type}/{model_name}"
        )
        nonzero_gene = count_values.sum(axis=0) > 0
        count_values = count_values[:, nonzero_gene]
        genes = genes[nonzero_gene]
        if len(genes) == 0:
            details["reason"] = "all genes in the analysis universe have zero counts"
            return details
        nonzero_library = count_values.sum(axis=1) > 0
        if nonzero_library.all():
            details["n_all_zero_genes_excluded"] = int((~nonzero_gene).sum())
            break
        zero_library_metadata = metadata.iloc[np.flatnonzero(~nonzero_library)]
        zero_library_samples.extend(
            {
                "pseudobulk_id": str(index),
                "study": str(row["study"]),
                "sample": str(row["sample"]),
            }
            for index, row in zero_library_metadata.iterrows()
        )
        metadata = metadata.iloc[np.flatnonzero(nonzero_library)].copy()
        if metadata.empty:
            details["reason"] = "all samples have zero counts over the fit's gene universe"
            return details
    if zero_library_samples:
        details["n_zero_library_samples_excluded"] = len(zero_library_samples)
        details["zero_library_samples_excluded"] = zero_library_samples
    if has_site_term and metadata["study_site"].nunique() < 2:
        design = "~ " + " + ".join(design_covariates)
        details["design"] = design
    try:
        results_by_covariate = _fit_covariate_model(
            count_values,
            genes,
            metadata,
            design=design,
            covariates=covariates,
            alpha=alpha,
            cpus=cpus,
        )
    except _NonEstimableModel as exc:
        details.update({"status": "skipped", "reason": str(exc)})
        return details
    except ValueError as exc:
        # Formula construction and DESeq validation use ValueError for models
        # that cannot be estimated from this subset; keep other fits runnable.
        details.update({"status": "skipped", "reason": str(exc)})
        return details
    except Exception as exc:  # noqa: BLE001 - report one failed fit and continue other models
        details.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        return details

    cell_type_dir = output_dir / _slug(cell_type)
    paths = {}
    test_only = test_mode or synthetic_test_data
    for covariate, result in results_by_covariate.items():
        result.insert(0, "cell_type", cell_type)
        result.insert(1, "model", model_name)
        result.insert(2, "study", study if study is not None else "all")
        result.insert(3, "analysis_mode", "test_only" if test_only else "standard")
        result.insert(
            4,
            "interpretation_warning",
            _SYNTHETIC_DATA_WARNING if synthetic_test_data else (
                _TEST_MODE_WARNING if test_mode else ""
            ),
        )
        result["studies_included"] = ";".join(studies)
        result["study_sample_counts"] = json.dumps(study_sample_counts, sort_keys=True)
        result["study_sites_included"] = ";".join(
            sorted(metadata["study_site"].astype(str).unique())
        )
        result["study_site_sample_counts"] = json.dumps(
            study_site_sample_counts, sort_keys=True
        )
        result["design"] = design
        result["n_samples_model"] = len(metadata)
        result["n_samples_before_complete_case_filter"] = n_samples_before_complete_case_filter
        result["n_samples_excluded_missing_design_covariates"] = int((~valid).sum())
        if model_name == "per_study":
            if covariate == "age":
                destination = cell_type_dir / "per_study" / f"{_slug(study or 'unknown-study')}.csv"
            else:
                destination = cell_type_dir / "per_study" / covariate / f"{_slug(study or 'unknown-study')}.csv"
        elif covariate == "age":
            destination = cell_type_dir / "merged.csv"
        else:
            destination = cell_type_dir / "combined" / f"{covariate}.csv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(destination, index=False)
        paths[covariate] = str(destination.relative_to(output_dir))

    # Keep the established age-result path in the summary fields while
    # indexing every fitted contrast explicitly for reports and downstream tools.
    primary_covariate = "age" if "age" in paths else next(iter(paths), None)
    details["results_by_covariate"] = paths
    details["results"] = paths.get(primary_covariate) if primary_covariate else None
    if not paths:
        return details
    significant = {}
    for covariate, result in results_by_covariate.items():
        significant[covariate] = int((result["padj"] < alpha).sum())
    details.update(
        {
            "status": "complete",
            "n_genes": len(genes),
            "n_significant_adjusted_p_lt_alpha": significant.get(primary_covariate, 0),
            "n_significant_by_covariate": significant,
        }
    )
    return details


def _write_age_model_diagnostics(
    *, adata, obs: pd.DataFrame, metadata: pd.DataFrame, cell_type: str, output_dir: Path
) -> Path | None:
    """Write the sample-level inputs needed for one shared-age-fit diagnostic figure."""
    if metadata.empty:
        return None
    studies = sorted(metadata["study"].astype(str).unique())
    var_mask = _availability_mask(adata.var, studies, context=f"{cell_type}/merged diagnostics")
    genes = adata.var_names[var_mask]
    if len(genes) == 0:
        return None
    row_positions = obs.index.get_indexer(metadata.index)
    col_positions = adata.var_names.get_indexer(genes)
    counts = _matrix_to_counts(
        adata.X[row_positions, :][:, col_positions], context=f"{cell_type}/merged diagnostics"
    )
    diagnostics = metadata[[
        "study", "study_site", "sample", "age", "sex", "n_cells", "total_counts",
        "log10_total_counts",
    ]].copy()
    diagnostics.insert(0, "pseudobulk_id", diagnostics.index.astype(str))
    diagnostics["counts_in_age_gene_intersection"] = counts.sum(axis=1, dtype=np.int64)
    diagnostics["genes_in_age_intersection"] = len(genes)
    destination = output_dir / _slug(cell_type) / "age_model_diagnostics.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(destination, index=False)
    return destination


def _load_pseudobulk(input_path: Path, pipeline: dict):
    """Load and validate the merged pseudobulk input and its configured labels."""
    import anndata as ad

    spec = pipeline["differential_expression"]
    if spec["method"] != "pydeseq2":
        raise ValueError("Only the configured PyDESeq2 method is supported")
    split_by = spec["split_by"]
    adata = ad.read_h5ad(input_path)
    required = {
        "study", "study_site", "sample", "age", "sex", "n_cells", "total_counts", split_by,
    }
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"{input_path}: pseudobulk metadata columns are missing: {missing}")
    if not adata.obs_names.is_unique:
        raise ValueError(f"{input_path}: pseudobulk identifiers must be unique")
    if not adata.var_names.is_unique:
        raise ValueError(f"{input_path}: gene identifiers must be unique")
    labels = adata.obs[split_by].astype("string")
    if labels.isna().any() or (labels.str.strip() == "").any():
        raise ValueError(f"{input_path}: split column {split_by!r} contains missing labels")
    return adata, spec, labels


def _cell_type_index(labels: pd.Series) -> list[tuple[str, str]]:
    """Return sorted cell-type labels with unique, stable output slugs."""
    used_slugs: dict[str, str] = {}
    cell_types = []
    for cell_type in sorted(labels.unique().tolist()):
        label = str(cell_type)
        slug = _slug(label)
        if slug in used_slugs and used_slugs[slug] != label:
            raise ValueError(
                f"AIFI L2 labels collide as output names: {used_slugs[slug]!r} and {label!r}"
            )
        used_slugs[slug] = label
        cell_types.append((slug, label))
    return cell_types


def list_pseudobulk_cell_types(input_path: Path, pipeline: dict) -> list[tuple[str, str]]:
    """List configured split labels without creating any per-type pseudobulk files."""
    _, _, labels = _load_pseudobulk(input_path, pipeline)
    return _cell_type_index(labels)


def _run_cell_type_differential_expression(
    adata,
    spec: dict,
    labels: pd.Series,
    output_dir: Path,
    *,
    cell_type: str,
    cpus: int,
    test_mode: bool,
    synthetic_test_data: bool,
) -> dict:
    """Fit one maximal per-study model and covariate-specific combined models."""
    slug = _slug(cell_type)
    selected = labels == cell_type
    type_adata = adata[selected.to_numpy(), :]
    type_obs = type_adata.obs.copy()
    (output_dir / slug).mkdir(parents=True, exist_ok=True)
    metadata, exclusions, excluded_samples = _prepare_metadata(
        type_obs, spec["sample_inclusion"], test_mode=test_mode
    )
    item = {
        "cell_type": cell_type,
        "slug": slug,
        "n_pseudobulks": int(type_adata.n_obs),
        "metadata_exclusions": exclusions,
        "excluded_samples": excluded_samples,
        "models": [],
    }
    diagnostics_path = _write_age_model_diagnostics(
        adata=type_adata,
        obs=type_obs,
        metadata=metadata,
        cell_type=cell_type,
        output_dir=output_dir,
    )
    if diagnostics_path is not None:
        item["age_model_diagnostics"] = str(diagnostics_path.relative_to(output_dir))
    optional_covariates = spec.get("optional_covariates", ["bmi", "cmv"])
    for study in sorted(metadata["study"].astype(str).unique()):
        study_metadata = metadata.loc[metadata["study"].astype(str) == study].copy()
        covariates = _maximal_per_study_covariates(study_metadata, optional_covariates)
        item["models"].append(
            _run_model(
                adata=type_adata,
                obs=type_obs,
                metadata=study_metadata,
                cell_type=cell_type,
                model_name="per_study",
                study=study,
                design=_design_with_covariates(spec["per_study_design"], covariates),
                covariates=covariates,
                alpha=spec["alpha"],
                cpus=cpus,
                output_dir=output_dir,
                test_mode=test_mode,
                synthetic_test_data=synthetic_test_data,
            )
        )
    core_covariates = [
        covariate for covariate in ("age", "sex", "log10_total_counts")
        if _has_covariate_variation(metadata, covariate)
    ] or ["age", "sex"]
    item["models"].append(
        _run_model(
            adata=type_adata,
            obs=type_obs,
            metadata=metadata,
            cell_type=cell_type,
            model_name="merged",
            study=None,
            design=spec["merged_design"],
            covariates=core_covariates,
            alpha=spec["alpha"],
            cpus=cpus,
            output_dir=output_dir,
            test_mode=test_mode,
            synthetic_test_data=synthetic_test_data,
        )
    )
    for covariate in optional_covariates:
        if covariate not in metadata:
            continue
        item["models"].append(
            _run_model(
                adata=type_adata,
                obs=type_obs,
                metadata=metadata,
                cell_type=cell_type,
                model_name="combined",
                study=None,
                design=f"{spec['merged_design']} + {covariate}",
                covariates=[covariate],
                alpha=spec["alpha"],
                cpus=cpus,
                output_dir=output_dir,
                test_mode=test_mode,
                synthetic_test_data=synthetic_test_data,
            )
        )
    trajectory_settings = AgeTrajectorySettings.from_mapping(spec["age_trajectory"])
    trajectory_metadata, trajectory_support = prepare_age_trajectory_metadata(
        metadata, trajectory_settings
    )
    if trajectory_metadata is None:
        trajectory_fit = {
            "model": "combined",
            "purpose": "age_trajectory",
            "study": None,
            "design": age_trajectory_design(spec),
            "covariates": ["age_bin"],
            "n_samples": trajectory_support["n_samples_retained"],
            "status": "skipped",
            "reason": trajectory_support["reason"],
        }
    else:
        trajectory_fit = _run_model(
            adata=type_adata,
            obs=type_obs,
            metadata=trajectory_metadata,
            cell_type=cell_type,
            model_name="combined",
            study=None,
            design=age_trajectory_design(spec),
            covariates=["age_bin"],
            alpha=spec["alpha"],
            cpus=cpus,
            output_dir=output_dir,
            test_mode=test_mode,
            synthetic_test_data=synthetic_test_data,
        )
        trajectory_fit["purpose"] = "age_trajectory"
        if trajectory_fit["status"] == "complete":
            trajectory_support["status"] = "complete"
        else:
            trajectory_support["status"] = "skipped"
            trajectory_support["reason"] = trajectory_fit.get("reason", "age-bin model did not fit")
    trajectory_fit.update({
        key: value for key, value in trajectory_support.items()
        if key not in {"status", "reason"}
    })
    item["models"].append(trajectory_fit)
    item["age_trajectory"] = {
        key: trajectory_support.get(key)
        for key in (
            "status", "reason", "reference_bin", "bin_width_years",
            "minimum_samples_per_bin", "minimum_bins", "age_bin_counts",
            "retained_bins", "dropped_bins", "samples_in_dropped_bins",
            "n_samples_retained", "reference_samples_by_study",
            "samples_by_study_and_bin",
        )
        if key in trajectory_support
    }
    type_metadata = {
        "kind": "pseudobulk_differential_expression_cell_type_metadata",
        "cell_type": cell_type,
        "analysis_mode": "test_only" if test_mode or synthetic_test_data else "standard",
        "test_mode": test_mode,
        "synthetic_test_data": synthetic_test_data,
        "interpretation_warning": (
            _SYNTHETIC_DATA_WARNING if synthetic_test_data else (
                _TEST_MODE_WARNING if test_mode else None
            )
        ),
        "sample_inclusion_configured": spec["sample_inclusion"],
        "sample_inclusion_applied": {
            **spec["sample_inclusion"],
            "minimum_cells_filter_enabled": not test_mode,
        },
        "age_trajectory": item.get("age_trajectory"),
        "metadata_exclusions": exclusions,
        "models": [
            {
                **model,
                "results_by_covariate": {
                    covariate: str(Path(path).relative_to(slug))
                    for covariate, path in model.get("results_by_covariate", {}).items()
                },
            }
            for model in item["models"]
        ],
    }
    cell_type_dir = output_dir / slug
    (cell_type_dir / "run_metadata.json").write_text(json.dumps(type_metadata, indent=2) + "\n")
    (cell_type_dir / "cell_type_result.json").write_text(json.dumps(item, indent=2) + "\n")
    return item


def _write_differential_expression_manifest(
    input_path: Path,
    output_dir: Path,
    pipeline: dict,
    records: list[dict],
    *,
    test_mode: bool,
    synthetic_test_data: bool,
) -> dict:
    """Write the stable top-level DE index after all cell-type tasks finish."""
    spec = pipeline["differential_expression"]
    split_by = spec["split_by"]
    output_dir.mkdir(parents=True, exist_ok=True)
    records = sorted(records, key=lambda item: item["cell_type"])

    report = {
        "schema_version": 1,
        "kind": "pseudobulk_differential_expression",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "configuration_sha256": config_digest(pipeline),
        "method": "PyDESeq2",
        "method_version": version("pydeseq2"),
        "split_by": split_by,
        "alpha": spec["alpha"],
        "analysis_mode": "test_only" if test_mode or synthetic_test_data else "standard",
        "test_mode": test_mode,
        "synthetic_test_data": synthetic_test_data,
        "interpretation_warning": (
            _SYNTHETIC_DATA_WARNING if synthetic_test_data else (
                _TEST_MODE_WARNING if test_mode else None
            )
        ),
        "gene_universe": {
            "per_study": "genes available in that study",
            "merged": "intersection of genes available in all studies included in each cell-type fit",
        },
        "optional_covariates": spec.get("optional_covariates", ["bmi", "cmv"]),
        "sample_inclusion": spec["sample_inclusion"],
        "sample_inclusion_applied": {
            **spec["sample_inclusion"],
            "minimum_cells_filter_enabled": not test_mode,
        },
        "models": {
            "per_study": (
                "One maximal complete-case design per study: study_site when it varies, "
                "the configured age/sex terms, log10_total_counts, and every configured "
                "optional covariate with at least two observed values."
            ),
            "merged": spec["merged_design"],
            "combined": (
                "One combined model per target covariate, adjusted for study_site, age, sex, "
                "and log10_total_counts; it uses complete cases from studies with recorded "
                "target values."
            ),
            "missing_values": (
                "A model excludes only samples missing a term in that model's design; exclusion "
                "counts are recorded with each fit."
            ),
        },
        "effect_scales": {
            "age": "log2 fold change per one-year increase",
            "bmi": "log2 fold change per one-unit increase in BMI",
            "log10_total_counts": "log2 fold change per one-unit increase in log10 pseudobulk total counts",
            "sex": "log2 fold change for the listed level versus female",
            "cmv": "log2 fold change for yes versus no CMV status",
        },
        "cell_types": records,
        "results_by_cell_type": {
            item["slug"]: {
                "cell_type": item["cell_type"],
                "analysis_mode": "test_only" if test_mode or synthetic_test_data else "standard",
                "test_mode": test_mode,
                "synthetic_test_data": synthetic_test_data,
                "interpretation_warning": (
                    _SYNTHETIC_DATA_WARNING if synthetic_test_data else (
                        _TEST_MODE_WARNING if test_mode else None
                    )
                ),
                "per_study": [
                    {
                        "study": study,
                        "paths": {
                            covariate: path
                            for model in item["models"]
                            if model["model"] == "per_study"
                            and model.get("study") == study
                            and model.get("status") == "complete"
                            for covariate, path in model.get("results_by_covariate", {}).items()
                        },
                    }
                    for study in sorted({
                        model["study"] for model in item["models"]
                        if model["model"] == "per_study"
                    })
                ],
                "combined": {
                    covariate: path
                    for model in item["models"]
                    if model["model"] in {"merged", "combined"}
                    and model.get("status") == "complete"
                    for covariate, path in model.get("results_by_covariate", {}).items()
                },
                "merged": next(
                    (
                        model.get("results_by_covariate", {}).get("age")
                        for model in item["models"]
                        if model["model"] == "merged" and model.get("status") == "complete"
                    ), None,
                ),
                "fit_status": (
                    "complete"
                    if any(model.get("status") == "complete" for model in item["models"])
                    else "failed"
                    if any(model.get("status") == "failed" for model in item["models"])
                    else "no_estimable_models"
                ),
            }
            for item in records
        },
    }
    completed = sum(
        model["status"] == "complete"
        for cell_type in records
        for model in cell_type["models"]
    )
    failed = sum(
        model["status"] == "failed"
        for cell_type in records
        for model in cell_type["models"]
    )
    if completed:
        status = "partial" if failed else "complete"
    elif failed:
        status = "failed"
    else:
        status = "no_estimable_models"
    report.update(
        {"status": status, "models_completed": int(completed), "models_failed": int(failed)}
    )
    (output_dir / "differential_expression.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    if failed and not completed:
        raise RuntimeError("No per-study or merged differential-expression model completed")
    return report


def run_differential_expression(
    input_path: Path,
    output_dir: Path,
    pipeline: dict,
    *,
    cpus: int = 1,
    test_mode: bool = False,
    synthetic_test_data: bool = False,
) -> dict:
    """Fit all cell types serially; retained for direct CLI/API use."""
    adata, spec, labels = _load_pseudobulk(input_path, pipeline)
    records = [
        _run_cell_type_differential_expression(
            adata, spec, labels, output_dir,
            cell_type=cell_type,
            cpus=cpus,
            test_mode=test_mode,
            synthetic_test_data=synthetic_test_data,
        )
        for _, cell_type in _cell_type_index(labels)
    ]
    return _write_differential_expression_manifest(
        input_path, output_dir, pipeline, records,
        test_mode=test_mode,
        synthetic_test_data=synthetic_test_data,
    )


def run_cell_type_differential_expression(
    input_path: Path,
    output_dir: Path,
    pipeline: dict,
    *,
    cell_type: str,
    cpus: int = 1,
    test_mode: bool = False,
    synthetic_test_data: bool = False,
) -> dict:
    """Load the merged pseudobulk H5AD and fit one cell type in this task."""
    adata, spec, labels = _load_pseudobulk(input_path, pipeline)
    available = dict(_cell_type_index(labels))
    if cell_type not in available.values():
        raise ValueError(f"{input_path}: requested cell type is absent: {cell_type!r}")
    return _run_cell_type_differential_expression(
        adata, spec, labels, output_dir,
        cell_type=cell_type,
        cpus=cpus,
        test_mode=test_mode,
        synthetic_test_data=synthetic_test_data,
    )


def combine_cell_type_differential_expression(
    input_path: Path,
    output_dir: Path,
    pipeline: dict,
    result_dirs: list[Path],
    *,
    test_mode: bool = False,
    synthetic_test_data: bool = False,
) -> dict:
    """Collect parallel task results into the established DE output contract."""
    import shutil

    expected = dict(list_pseudobulk_cell_types(input_path, pipeline))
    records = []
    seen_slugs = set()
    output_dir.mkdir(parents=True, exist_ok=True)
    for result_dir in result_dirs:
        record_path = result_dir / "cell_type_result.json"
        if not record_path.is_file():
            raise FileNotFoundError(f"{result_dir}: missing cell_type_result.json")
        item = json.loads(record_path.read_text())
        slug = item.get("slug")
        cell_type = item.get("cell_type")
        if slug not in expected or expected[slug] != cell_type:
            raise ValueError(f"{result_dir}: unexpected cell-type result {cell_type!r}/{slug!r}")
        if slug in seen_slugs:
            raise ValueError(f"Duplicate cell-type result for {slug!r}")
        seen_slugs.add(slug)
        shutil.copytree(result_dir, output_dir / slug)
        records.append(item)
    missing = sorted(set(expected) - seen_slugs)
    if missing:
        raise ValueError(f"Missing cell-type DE results: {missing}")
    return _write_differential_expression_manifest(
        input_path, output_dir, pipeline, records,
        test_mode=test_mode,
        synthetic_test_data=synthetic_test_data,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run per-study and combined PyDESeq2 covariate models on merged pseudobulk"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=1)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--list-cell-types",
        action="store_true",
        help="write the configured pseudobulk cell-type labels and slugs as TSV",
    )
    mode.add_argument(
        "--cell-type",
        help="fit only this cell type after loading and subsetting the merged pseudobulk H5AD",
    )
    mode.add_argument(
        "--combine-cell-type-results",
        type=Path,
        nargs="+",
        metavar="DIRECTORY",
        help="combine directories emitted by parallel --cell-type tasks",
    )
    parser.add_argument(
        "--cell-types-output",
        type=Path,
        help="TSV destination required with --list-cell-types",
    )
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help=(
            "bypass only the minimum-cells-per-pseudobulk filter to exercise "
            "small test fixtures; all outputs are marked test-only"
        ),
    )
    parser.add_argument(
        "--synthetic-test-data",
        action="store_true",
        help="label results as synthetic test data without bypassing production inclusion filters",
    )
    args = parser.parse_args()
    if args.cpus < 1:
        parser.error("--cpus must be positive")
    pipeline = read_json(args.config)
    if args.list_cell_types:
        if args.cell_types_output is None:
            parser.error("--list-cell-types requires --cell-types-output")
        records = list_pseudobulk_cell_types(args.input, pipeline)
        args.cell_types_output.write_text(
            "".join(f"{slug}\t{cell_type}\n" for slug, cell_type in records)
        )
        LOGGER.info("Wrote %s pseudobulk cell types to %s", len(records), args.cell_types_output)
        return
    if args.output_dir is None:
        parser.error("--output-dir is required when fitting or combining DE results")
    if args.combine_cell_type_results is not None:
        combine_cell_type_differential_expression(
            args.input,
            args.output_dir,
            pipeline,
            args.combine_cell_type_results,
            test_mode=args.test_mode,
            synthetic_test_data=args.synthetic_test_data,
        )
    elif args.cell_type:
        run_cell_type_differential_expression(
            args.input,
            args.output_dir,
            pipeline,
            cell_type=args.cell_type,
            cpus=args.cpus,
            test_mode=args.test_mode,
            synthetic_test_data=args.synthetic_test_data,
        )
    else:
        run_differential_expression(
            args.input,
            args.output_dir,
            pipeline,
            cpus=args.cpus,
            test_mode=args.test_mode,
            synthetic_test_data=args.synthetic_test_data,
        )
    LOGGER.info("Differential-expression results written to %s", args.output_dir)


if __name__ == "__main__":
    main()

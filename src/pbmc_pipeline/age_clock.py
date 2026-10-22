"""Train and apply site-held-out, cell-type pseudobulk age clocks."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import traceback
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import ParameterGrid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from pbmc_pipeline.age_clock_features import (
    WithinSampleRankTransformer,
    within_study_gene_zscore,
)
from pbmc_pipeline.config import config_digest
from pbmc_pipeline.utils import slugify

LOGGER = logging.getLogger(__name__)
MODEL_FORMAT_VERSION = 3
MISSING_COEFFICIENT_WEIGHT_WARNING_FRACTION = 0.05


class ClockEligibilityError(ValueError):
    """A cell type cannot support the configured study-level clock split."""


class ClockFitFailuresError(RuntimeError):
    """One or more all-cell-type fits failed after other types completed."""


@dataclass
class ClockFit:
    """In-memory fit result, deliberately separate from artifact publication."""

    manifest: dict
    inclusion_counts: dict
    fitted_methods: list[dict]


@dataclass(frozen=True)
class ClockMethod:
    """Small declarative method specification consumed by the shared evaluator."""

    name: str
    display_name: str
    description: str
    estimator_factory: Callable[[], Pipeline]
    parameter_grid: dict[str, list]
    default_parameters: dict[str, object]
    estimator_step: str = "regression"
    coefficient_step: str | None = "regression"
    feature_representation: str = "log-CPM expression"
    scaling_description: str = (
        "training-only standard deviation; sparse; model intercept fitted"
    )
    equal_study_weight: bool = False
    within_study_zscore: bool = False
    weighted_steps: tuple[str, ...] = ("scale", "regression")
    gene_selector: Callable[[object, list[str]], list[str]] | None = None


def _elastic_net_pipeline() -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler(with_mean=False)),
        ("regression", ElasticNet(max_iter=10000, selection="cyclic")),
    ])


def _rank_elastic_net_pipeline() -> Pipeline:
    return Pipeline([
        ("rank", WithinSampleRankTransformer()),
        ("scale", StandardScaler()),
        ("regression", ElasticNet(max_iter=10000, selection="cyclic")),
    ])


def _clock_methods(alphas: list[float], l1_ratios: list[float]) -> list[ClockMethod]:
    common_grid = {
        "regression__alpha": alphas,
        "regression__l1_ratio": l1_ratios,
    }
    return [
        ClockMethod(
            name="elastic_net",
            display_name="Elastic net",
            description="Each training pseudobulk has equal weight.",
            estimator_factory=_elastic_net_pipeline,
            parameter_grid=common_grid,
            default_parameters={"regression__alpha": 0.1, "regression__l1_ratio": 0.5},
        ),
        ClockMethod(
            name="equal_study_elastic_net",
            display_name="Study-balanced elastic net",
            description=(
                "Each training study has equal total weight; samples in larger "
                "studies receive proportionally less weight."
            ),
            estimator_factory=_elastic_net_pipeline,
            parameter_grid=common_grid,
            default_parameters={"regression__alpha": 0.1, "regression__l1_ratio": 0.5},
            equal_study_weight=True,
        ),
        ClockMethod(
            name="rank_elastic_net",
            display_name="Within-sample rank elastic net",
            description=(
                "Each sample is represented by percentile ranks of genes within "
                "the training gene universe; each training sample has equal weight."
            ),
            estimator_factory=_rank_elastic_net_pipeline,
            parameter_grid=common_grid,
            default_parameters={"regression__alpha": 0.1, "regression__l1_ratio": 0.5},
            feature_representation="within-sample gene percentile ranks",
            scaling_description=(
                "within-sample percentile ranks, then training-only centered "
                "standard deviation; model intercept fitted"
            ),
        ),
        ClockMethod(
            name="study_zscore_elastic_net",
            display_name="Within-study gene z-score elastic net",
            description=(
                "Each gene is centered and scaled across samples independently "
                "within each study before fitting and prediction."
            ),
            estimator_factory=_elastic_net_pipeline,
            parameter_grid=common_grid,
            default_parameters={"regression__alpha": 0.1, "regression__l1_ratio": 0.5},
            feature_representation="within-study per-gene z-scores of log-CPM expression",
            scaling_description=(
                "per-gene study-cohort z-score, then training-only standard deviation; "
                "target-study cohort expression is required at prediction"
            ),
            within_study_zscore=True,
        ),
    ]


def _load_config(path: Path) -> dict:
    return json.loads(path.read_text())


def _raw_library_sizes(adata, row_positions: np.ndarray) -> np.ndarray:
    """Validate raw counts and return their full-row totals for sample inclusion."""
    matrix = adata[row_positions, :].X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Age-clock expression must contain finite, non-negative counts")
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    if not np.isfinite(totals).all():
        raise ValueError("Age-clock raw-count library sizes must be finite")
    return totals


def _metadata_for_type(adata, pipeline: dict, cell_type: str) -> tuple[pd.DataFrame, dict]:
    spec = pipeline["differential_expression"]
    inclusion = spec["sample_inclusion"]
    split_by = spec["split_by"]
    required = {"study", "sample", "age", "n_cells", split_by}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"Pseudobulk metadata is missing required fields: {missing}")

    obs = adata.obs.loc[adata.obs[split_by].astype(str) == cell_type].copy()
    obs["age"] = pd.to_numeric(obs["age"], errors="coerce")
    obs["n_cells"] = pd.to_numeric(obs["n_cells"], errors="coerce")
    obs["total_counts"] = _raw_library_sizes(adata, adata.obs_names.get_indexer(obs.index))
    keep = pd.Series(True, index=obs.index)
    counts = {"n_cell_type_samples": len(obs)}
    excluded_masks = {}
    age_ok = np.isfinite(obs["age"]) & obs["age"].ge(inclusion["minimum_age_years_inclusive"])
    excluded_masks["excluded_invalid_or_below_minimum_age"] = ~age_ok
    keep &= age_ok
    cells_ok = np.isfinite(obs["n_cells"]) & obs["n_cells"].ge(
        inclusion["minimum_cells_per_pseudobulk_inclusive"]
    )
    excluded_masks["excluded_invalid_or_below_minimum_cells"] = keep & ~cells_ok
    keep &= cells_ok
    library_ok = np.isfinite(obs["total_counts"]) & obs["total_counts"].gt(0)
    excluded_masks["excluded_invalid_or_zero_library"] = keep & ~library_ok
    keep &= library_ok
    exclusions = inclusion.get("exclude_subject_prefixes_by_study", {})
    present_exclusion_studies = set(obs["study"].astype(str)) & set(exclusions)
    if present_exclusion_studies:
        if "subject" not in obs:
            raise ValueError("Configured subject-prefix exclusions require pseudobulk 'subject'")
        subject = obs["subject"].astype("string")
        excluded_subjects = pd.Series(False, index=obs.index)
        for study in sorted(present_exclusion_studies):
            prefixes = exclusions[study]
            in_study = obs["study"].astype(str) == study
            for prefix in prefixes:
                excluded = in_study & subject.str.startswith(prefix, na=False)
                excluded_subjects |= keep & excluded
                keep &= ~excluded
        excluded_masks["excluded_configured_subject_prefix"] = excluded_subjects
    counts.update({name: int(mask.sum()) for name, mask in excluded_masks.items()})
    counts["n_eligible_samples"] = int(keep.sum())
    counts["by_study"] = {}
    for study in sorted(adata.obs["study"].astype(str).unique()):
        in_study = obs["study"].astype(str) == study
        counts["by_study"][study] = {
            "n_cell_type_samples": int(in_study.sum()),
            **{name: int((mask & in_study).sum()) for name, mask in excluded_masks.items()},
            "n_eligible_samples": int((keep & in_study).sum()),
        }
    return obs.loc[keep].copy(), counts


def _training_gene_names(adata, studies: list[str]) -> list[str]:
    """Find the gene intersection using training studies only."""
    var = adata.var
    availability_columns = [f"available_in_{study}" for study in studies]
    present_columns = [column for column in availability_columns if column in var]
    if len(present_columns) != len(studies):
        raise ValueError(
            "Pseudobulk gene availability is incomplete for the requested training studies"
        )
    available = var[present_columns].fillna(False).astype(bool).all(axis=1)
    genes = var.index[available].astype(str).tolist()
    if not genes:
        raise ValueError("No genes are shared by all training studies")
    return genes


def _training_gene_intersection_diagnostics(
    adata, studies: list[str], model_genes: list[str],
) -> dict:
    """Describe the training-only gene loss caused by the final shared universe."""
    availability_columns = [f"available_in_{study}" for study in studies]
    if any(column not in adata.var for column in availability_columns):
        raise ValueError("Pseudobulk gene availability is incomplete for training studies")
    availability = adata.var[availability_columns].fillna(False).astype(bool)
    n_union = int(availability.any(axis=1).sum())
    n_model = len(model_genes)
    return {
        "n_genes_available_in_any_training_study": n_union,
        "n_genes_in_final_model_universe": n_model,
        "n_genes_dropped_from_training_union": n_union - n_model,
        "by_study": [
            {
                "study": study,
                "n_genes_available_before_intersection": int(availability[column].sum()),
                "n_genes_dropped_from_study_universe": int(availability[column].sum()) - n_model,
            }
            for study, column in zip(studies, availability_columns, strict=True)
        ],
    }


def _gene_coverage(adata, row_positions: np.ndarray, genes: list[str]) -> list[dict]:
    """Describe assay coverage without changing the training-defined gene universe."""
    if not genes:
        raise ValueError("Age-clock model gene universe is empty")
    present = [gene for gene in genes if gene in adata.var_names]
    obs = adata.obs.iloc[row_positions]
    studies = obs["study"].astype(str) if "study" in obs else pd.Series("input", index=obs.index)
    records = []
    for study in sorted(studies.unique()):
        column = f"available_in_{study}"
        n_available = (
            int(adata.var.loc[present, column].fillna(False).astype(bool).sum())
            if column in adata.var else len(present)
        )
        records.append({
            "study": study,
            "n_samples": int((studies == study).sum()),
            "n_model_genes": len(genes),
            "n_available_genes": n_available,
            "n_missing_genes": len(genes) - n_available,
            "fraction_available": n_available / len(genes),
        })
    return records


def _coefficient_coverage(
    adata, genes: list[str], coefficients: np.ndarray | None, study: str,
) -> dict:
    """Quantify fitted coefficient mass unavailable to one study's assay.

    Coefficients are on the estimator's scaled-feature scale.  Absolute mass is
    deliberately used as a transparent, linear proxy for how much of the fitted
    predictor is forced to zero; it does not require assumptions about feature
    covariance in the new study.
    """
    if coefficients is None:
        return {
            "missing_coefficient_absolute_weight": None,
            "total_coefficient_absolute_weight": None,
            "missing_coefficient_weight_fraction": None,
            "missing_coefficient_weight_warning": False,
        }
    coefficient_array = np.asarray(coefficients, dtype=float).ravel()
    if len(coefficient_array) != len(genes):
        raise ValueError("Age-clock coefficient count does not match model genes")
    present = np.asarray([gene in adata.var_names for gene in genes], dtype=bool)
    column = f"available_in_{study}"
    if column in adata.var:
        available = np.zeros(len(genes), dtype=bool)
        present_positions = np.flatnonzero(present)
        available[present_positions] = (
            adata.var.loc[np.asarray(genes, dtype=object)[present_positions], column]
            .fillna(False).astype(bool).to_numpy()
        )
    else:
        available = present
    absolute = np.abs(coefficient_array)
    total = float(absolute.sum())
    missing = float(absolute[~available].sum())
    fraction = missing / total if total > 0 else 0.0
    return {
        "missing_coefficient_absolute_weight": missing,
        "total_coefficient_absolute_weight": total,
        "missing_coefficient_weight_fraction": fraction,
        "missing_coefficient_weight_warning": (
            fraction > MISSING_COEFFICIENT_WEIGHT_WARNING_FRACTION
        ),
    }


def _log_cpm(adata, row_positions: np.ndarray, genes: list[str]):
    """Normalize over available model genes, then align and zero-fill missing genes."""
    for coverage in _gene_coverage(adata, row_positions, genes):
        if coverage["n_available_genes"] == 0:
            raise ValueError(
                f"No model genes are available for study {coverage['study']}; "
                "age prediction requires nonzero gene overlap"
            )
        if coverage["n_missing_genes"]:
            LOGGER.warning(
                "Study %s supplies %d/%d model genes (%.1f%%); missing genes are zero-filled",
                coverage["study"], coverage["n_available_genes"], coverage["n_model_genes"],
                100 * coverage["fraction_available"],
            )
    _raw_library_sizes(adata, row_positions)
    present = [gene for gene in genes if gene in adata.var_names]
    matrix = adata[row_positions, present].X
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr(copy=True)
    library_sizes = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    if not np.isfinite(library_sizes).all() or (library_sizes <= 0).any():
        raise ValueError(
            "Age-clock samples must have positive count totals across available model genes"
        )
    matrix = sparse.diags(1_000_000.0 / library_sizes) @ matrix
    matrix = matrix.tocsr()
    matrix.data = np.log1p(matrix.data)
    if present == genes:
        return matrix
    gene_positions = {gene: index for index, gene in enumerate(genes)}
    remap = np.asarray([gene_positions[gene] for gene in present], dtype=np.int32)
    matrix = matrix.tocoo()
    aligned = sparse.csr_matrix(
        (matrix.data, (matrix.row, remap[matrix.col])),
        shape=(matrix.shape[0], len(genes)),
    )
    return aligned


def predict_age_clock(model_path: Path, adata, row_positions: np.ndarray) -> np.ndarray:
    """Predict from raw pseudobulks, rejecting assays with zero model-gene overlap."""
    bundle = joblib.load(model_path)
    if bundle.get("model_format_version") != MODEL_FORMAT_VERSION:
        raise ValueError(
            "Age-clock model uses a deprecated normalization contract; refit the model"
        )
    genes = bundle["genes"]
    matrix = _log_cpm(adata, row_positions, genes)
    if bundle.get("within_study_zscore", False):
        studies = adata.obs.iloc[row_positions]["study"].astype(str).to_numpy()
        matrix, _ = within_study_gene_zscore(matrix, studies)
    return bundle["estimator"].predict(matrix)


def _metrics(observed: pd.Series, predicted: np.ndarray) -> dict:
    result = {
        "n_samples": len(observed),
        "mae_years": float(mean_absolute_error(observed, predicted)),
        "rmse_years": float(np.sqrt(mean_squared_error(observed, predicted))),
    }
    result["r_squared"] = (
        float(r2_score(observed, predicted)) if len(observed) > 1 else None
    )
    if len(observed) > 1 and observed.nunique() > 1 and np.std(predicted) > 0:
        result["calibration_slope"] = float(np.polyfit(predicted, observed, 1)[0])
        result["calibration_intercept_years"] = float(
            np.polyfit(predicted, observed, 1)[1]
        )
    else:
        result["calibration_slope"] = None
        result["calibration_intercept_years"] = None
    return result


def _sample_weights(studies: pd.Series, method: ClockMethod) -> np.ndarray | None:
    """Return mean-one weights giving every training study equal total weight."""
    if not method.equal_study_weight:
        return None
    study_ids = studies.astype(str)
    counts = study_ids.value_counts()
    weights = study_ids.map(1.0 / counts).to_numpy(dtype=float)
    return weights * (len(weights) / weights.sum())


def _fit_estimator(
    method: ClockMethod,
    parameters: dict[str, object],
    x_train,
    y_train: np.ndarray,
    sample_weights: np.ndarray | None,
) -> Pipeline:
    estimator = method.estimator_factory().set_params(**parameters)
    if sample_weights is None:
        estimator.fit(x_train, y_train)
    else:
        fit_params = {
            f"{step}__sample_weight": sample_weights
            for step in method.weighted_steps
        }
        estimator.fit(
            x_train,
            y_train,
            **fit_params,
        )
    return estimator


def _fit_diagnostics(estimator: Pipeline, method: ClockMethod) -> dict:
    model = estimator.named_steps[method.estimator_step]
    diagnostics = {}
    for attribute in ("n_iter_", "dual_gap_", "max_iter", "tol"):
        value = getattr(model, attribute, None)
        if value is not None:
            diagnostics[attribute.rstrip("_")] = np.asarray(value).tolist()
    if "n_iter" in diagnostics and "max_iter" in diagnostics:
        diagnostics["iteration_limit_reached"] = bool(
            np.any(np.asarray(diagnostics["n_iter"]) >= diagnostics["max_iter"])
        )
    return diagnostics


def _diagnostic_maximum(diagnostics: dict, name: str) -> float | None:
    """Return a scalar maximum from a sklearn scalar-or-array diagnostic."""
    value = diagnostics.get(name)
    if value is None:
        return None
    values = np.asarray(value, dtype=float)
    return float(np.max(values)) if values.size else None


def _optional_float(value: object) -> float | None:
    """Convert a persisted diagnostic value while retaining unavailable values."""
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _score_cv_candidate(
    method: ClockMethod,
    parameters: dict[str, object],
    x_fit,
    y_fit: np.ndarray,
    weights: np.ndarray | None,
    x_valid,
    y_valid: np.ndarray,
) -> tuple[float, dict]:
    """Fit and score one independent CV parameter candidate."""
    estimator = _fit_estimator(method, parameters, x_fit, y_fit, weights)
    score = float(mean_absolute_error(y_valid, estimator.predict(x_valid)))
    return score, _fit_diagnostics(estimator, method)


def _software_versions() -> dict:
    versions = {"python": platform.python_version()}
    for name in ("pbmc-ageing-pipeline", "anndata", "numpy", "scipy", "scikit-learn", "joblib"):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = "unknown"
    return versions


def _mean_study_mae(predictions: pd.DataFrame) -> float:
    per_study = predictions.groupby("study", observed=True)["age_error_years"].apply(
        lambda errors: float(np.abs(errors).mean())
    )
    return float(per_study.mean())


def _dense_features(matrix) -> np.ndarray:
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix, dtype=float)


def _fit_age_reference_expression(
    x, ages: np.ndarray, studies: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float, list[str]] | None:
    """Fit shared per-gene age trends with study-specific intercepts."""
    values = _dense_features(x)
    ages = np.asarray(ages, dtype=float)
    studies = np.asarray(studies, dtype=str)
    study_levels = sorted(np.unique(studies))
    if values.ndim != 2 or len(ages) != values.shape[0] or len(studies) != len(ages):
        raise ValueError("Age-reference features, ages, and studies must align")
    if len(study_levels) < 2 or len(np.unique(ages)) < 2:
        return None

    age_center = float(ages.mean())
    design = np.column_stack([
        np.ones(len(ages)),
        ages - age_center,
        *((studies == study).astype(float) for study in study_levels[1:]),
    ])
    study_counts = pd.Series(studies).value_counts().to_dict()
    row_weights = np.asarray([1.0 / study_counts[study] for study in studies])
    sqrt_weights = np.sqrt(row_weights)
    coefficients, _, rank, _ = np.linalg.lstsq(
        design * sqrt_weights[:, None],
        values * sqrt_weights[:, None],
        rcond=None,
    )
    if rank < design.shape[1]:
        return None
    study_intercept = (
        coefficients[0]
        + coefficients[2:].sum(axis=0) / len(study_levels)
    )
    return study_intercept, coefficients[1], age_center, study_levels


def _age_reference_values(
    reference: tuple[np.ndarray, np.ndarray, float, list[str]], ages: np.ndarray
) -> np.ndarray:
    """Evaluate the mean-across-studies gene expression trend at given ages."""
    intercept, age_slope, age_center, _ = reference
    return intercept[None, :] + (
        np.asarray(ages, dtype=float)[:, None] - age_center
    ) * age_slope[None, :]


def _target_expression_calibration(
    x_train,
    y_train: np.ndarray,
    train_studies: np.ndarray,
    x_external,
    y_external: np.ndarray,
    external_studies: np.ndarray,
    external_gene_available: np.ndarray,
    *,
    calibration_samples: int,
    repeats: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, dict]]:
    """Evaluate age-protected target-study mean shifts on repeated donor splits."""
    train_values = _dense_features(x_train)
    external_values = _dense_features(x_external)
    y_external = np.asarray(y_external, dtype=float)
    external_studies = np.asarray(external_studies, dtype=str)
    reference = _fit_age_reference_expression(
        train_values, y_train, train_studies
    )
    if reference is None:
        return pd.DataFrame(), {
            study: {"status": "non_estimable_training_age_reference"}
            for study in sorted(np.unique(external_studies))
        }

    records = []
    summaries = {}
    for study in sorted(np.unique(external_studies)):
        study_indices = np.flatnonzero(external_studies == study)
        n_study = len(study_indices)
        if n_study <= calibration_samples:
            summaries[study] = {
                "status": "too_few_samples_for_calibration_and_evaluation",
                "n_external_samples": n_study,
                "n_calibration_samples": calibration_samples,
                "n_evaluation_samples": max(n_study - calibration_samples, 0),
            }
            continue

        study_seed = seed + int.from_bytes(
            hashlib.sha256(study.encode("utf-8")).digest()[:4], "little"
        )
        rng = np.random.default_rng(study_seed)
        study_records = []
        for repeat in range(repeats):
            calibration_local = rng.choice(
                n_study, size=calibration_samples, replace=False
            )
            calibration_indices = study_indices[calibration_local]
            evaluation_mask = np.ones(n_study, dtype=bool)
            evaluation_mask[calibration_local] = False
            evaluation_indices = study_indices[evaluation_mask]

            calibration_reference = _age_reference_values(
                reference, y_external[calibration_indices]
            )
            residuals = external_values[calibration_indices] - calibration_reference
            site_shift = residuals.mean(axis=0)
            site_shift[~external_gene_available] = 0.0

            study_records.append({
                "study": study,
                "repeat": repeat + 1,
                "n_calibration_samples": calibration_samples,
                "n_evaluation_samples": len(evaluation_indices),
                "calibration_age_min": float(y_external[calibration_indices].min()),
                "calibration_age_max": float(y_external[calibration_indices].max()),
                "evaluation_age_min": float(y_external[evaluation_indices].min()),
                "evaluation_age_max": float(y_external[evaluation_indices].max()),
                "mean_absolute_expression_shift": float(
                    np.abs(site_shift[external_gene_available]).mean()
                ),
                "site_shift": site_shift,
                "evaluation_indices": evaluation_indices,
            })

        summaries[study] = {
            "status": "completed",
            "n_external_samples": n_study,
            "n_calibration_samples": calibration_samples,
            "n_evaluation_samples_per_repeat": n_study - calibration_samples,
            "n_repeats": repeats,
            "calibration_age_minimum_mean": float(np.mean([
                record["calibration_age_min"] for record in study_records
            ])),
            "calibration_age_maximum_mean": float(np.mean([
                record["calibration_age_max"] for record in study_records
            ])),
            "mean_absolute_expression_shift": float(np.mean([
                record["mean_absolute_expression_shift"] for record in study_records
            ])),
        }
        records.extend(study_records)

    # Feature matrices stay local to evaluation; publish only split summaries.
    return pd.DataFrame(records), summaries


def fit_age_clocks(
    input_path: Path,
    config_path: Path,
    external_studies: list[str],
    *,
    cell_type_filter: str | None = None,
    alphas: list[float] | None = None,
    l1_ratios: list[float] | None = None,
    method_names: list[str] | None = None,
    calibration_samples: int = 10,
    calibration_repeats: int = 50,
    calibration_seed: int = 42,
    cv: bool = True,
    n_jobs: int = 1,
) -> ClockFit:
    """Fit one cell type in memory without writing model artifacts."""
    import anndata as ad

    if not external_studies or any(not study for study in external_studies):
        raise ValueError("Specify one or more comma-separated external study IDs")
    if len(set(external_studies)) != len(external_studies):
        raise ValueError("External study IDs must be distinct")
    if n_jobs < 1:
        raise ValueError("n_jobs must be at least 1")
    if calibration_samples < 0:
        raise ValueError("calibration_samples must be zero or greater")
    if calibration_samples and calibration_repeats < 1:
        raise ValueError("calibration_repeats must be at least 1 when calibration is enabled")
    pipeline_config = _load_config(config_path)
    adata = ad.read_h5ad(input_path)
    available_studies = set(adata.obs["study"].astype(str).unique())
    unknown = sorted(set(external_studies) - available_studies)
    if unknown:
        raise ValueError(f"External studies are absent from the pseudobulk: {unknown}")
    all_studies = sorted(available_studies)
    training_studies = [study for study in all_studies if study not in external_studies]
    if not training_studies:
        raise ValueError("No training studies remain after external studies are held out")

    labels = pipeline_config["differential_expression"]["split_by"]
    available_cell_types = sorted(adata.obs[labels].astype(str).unique())
    requested = (cell_type_filter or "").strip()
    matches = [label for label in available_cell_types
               if label == requested or slugify(label) == requested]
    if not matches:
        choices = ", ".join(f"{label} ({slugify(label)})" for label in available_cell_types)
        raise ValueError(f"Unknown cell type {requested!r}. Available: {choices}")
    selected_cell_type = matches[0]
    if alphas is None:
        alphas = [0.001, 0.01, 0.1, 1.0, 10.0]
    if l1_ratios is None:
        l1_ratios = [0.1, 0.5, 0.9, 1.0]

    manifest = {
        "kind": "cell_type_transcriptomic_age_clock",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "configuration_sha256": config_digest(pipeline_config),
        "software_versions": _software_versions(),
        "normalization": (
            "CPM from available training-defined model genes, then log1p; "
            "missing model genes are zero-filled after normalization"
        ),
        "missing_gene_policy": "zero-fill missing model genes; reject zero overlap per study",
        "missing_coefficient_weight_warning_fraction": (
            MISSING_COEFFICIENT_WEIGHT_WARNING_FRACTION
        ),
        "synthetic_test_data": bool(adata.uns.get("synthetic_test_data", False)),
        "target": "continuous age in years",
        "methods": [],
        "preprocessing": (
            "per-sample CPM and log1p; scaling is fit on training samples only "
            "inside each fold and final development fit"
        ),
        "feature_selection": (
            "gene intersection is recalculated from inner-training studies in each "
            "CV fold and from all development studies for final fitting"
        ),
        "external_expression_use": "raw pseudobulk counts transformed only after training and tuning",
        "target_expression_calibration": {
            "method": "age-protected per-gene additive shift",
            "calibration_samples_per_target_study": calibration_samples,
            "repeated_splits": calibration_repeats if calibration_samples else 0,
            "seed": calibration_seed,
            "age_reference": (
                "per-gene linear age trend with study-specific intercepts and "
                "equal total study weights; reference intercepts averaged "
                "equally across training studies"
            ),
            "evaluation": "calibration donors excluded; compare raw and corrected predictions on the same held-out donors",
        } if calibration_samples else {"method": "disabled"},
        "training_studies": training_studies,
        "external_studies": external_studies,
        "cross_validation": "leave-one-training-study-out; fold-specific gene intersection" if cv else "disabled",
        "cv_parallel_jobs": n_jobs,
        "cell_type": selected_cell_type,
    }

    metadata, inclusion_counts = _metadata_for_type(
        adata, pipeline_config, selected_cell_type
    )
    train_meta = metadata.loc[metadata["study"].astype(str).isin(training_studies)].copy()
    external_meta = metadata.loc[metadata["study"].astype(str).isin(external_studies)].copy()
    missing_external = sorted(set(external_studies) - set(external_meta["study"].astype(str)))
    missing_training = sorted(set(training_studies) - set(train_meta["study"].astype(str)))
    if missing_external or missing_training:
        missing_studies = sorted(set(missing_external) | set(missing_training))
        details = {study: inclusion_counts["by_study"][study] for study in missing_studies}
        problems = []
        if missing_external:
            problems.append(
                f"External studies have no eligible samples for {selected_cell_type}: "
                f"{missing_external}"
            )
        if missing_training:
            problems.append(
                f"Training studies have no eligible samples for {selected_cell_type}: "
                f"{missing_training}"
            )
        raise ClockEligibilityError(
            "; ".join(problems)
            + f". Per-study inclusion counts: {json.dumps(details)}"
        )
    if train_meta.empty:
        raise ClockEligibilityError("Selected cell type has no eligible training samples")
    train_meta["study"] = train_meta["study"].astype(str)
    external_meta["study"] = external_meta["study"].astype(str)
    train_study_ids = sorted(train_meta["study"].unique())
    manifest["training_studies"] = train_study_ids
    manifest["sample_inclusion_counts"] = inclusion_counts
    if cv and len(train_study_ids) < 2:
        raise ValueError("Leave-one-study-out tuning needs at least two training studies")
    train_positions = adata.obs_names.get_indexer(train_meta.index)
    external_positions = adata.obs_names.get_indexer(external_meta.index)
    y_train = train_meta["age"].to_numpy(dtype=float)
    alphas = alphas or [0.001, 0.01, 0.1, 1.0, 10.0]
    l1_ratios = l1_ratios or [0.1, 0.5, 0.9, 1.0]
    available_methods = _clock_methods(alphas, l1_ratios)
    if method_names is None:
        methods = available_methods
    else:
        requested_methods = list(dict.fromkeys(name.strip() for name in method_names))
        unknown_methods = sorted(
            set(requested_methods) - {method.name for method in available_methods}
        )
        if unknown_methods:
            raise ValueError(
                f"Unknown age-clock method(s) {unknown_methods}; available methods: "
                f"{[method.name for method in available_methods]}"
            )
        if not requested_methods:
            raise ValueError("At least one age-clock method must be selected")
        methods = [method for method in available_methods if method.name in requested_methods]
    manifest["selected_methods"] = [method.name for method in methods]
    fitted = []

    for method in methods:
        LOGGER.info(
            "Starting clock %s (%s) for %s",
            method.display_name,
            method.name,
            selected_cell_type,
        )
        cv_rows = []
        if cv:
            fold_data = []
            fold_coverage = {}
            study_values = train_meta["study"].to_numpy(dtype=str)
            for held_out_study in train_study_ids:
                validation_mask = study_values == held_out_study
                fit_indices = np.flatnonzero(~validation_mask)
                validation_indices = np.flatnonzero(validation_mask)
                fold_train_meta = train_meta.iloc[fit_indices]
                fold_studies = sorted(fold_train_meta["study"].unique())
                gene_selector = method.gene_selector or _training_gene_names
                fold_genes = gene_selector(adata, fold_studies)
                x_fit = _log_cpm(adata, train_positions[fit_indices], fold_genes)
                x_validation = _log_cpm(
                    adata, train_positions[validation_indices], fold_genes
                )
                if method.within_study_zscore:
                    x_fit, _ = within_study_gene_zscore(
                        x_fit, fold_train_meta["study"].to_numpy()
                    )
                    x_validation, _ = within_study_gene_zscore(
                        x_validation,
                        train_meta.iloc[validation_indices]["study"].to_numpy(),
                    )
                fold_coverage[held_out_study] = _gene_coverage(
                    adata, train_positions[validation_indices], fold_genes
                )[0]
                fold_weights = _sample_weights(fold_train_meta["study"], method)
                fold_data.append((
                    held_out_study,
                    fold_genes,
                    x_fit,
                    y_train[fit_indices],
                    fold_weights,
                    x_validation,
                    y_train[validation_indices],
                ))

            parameter_candidates = list(ParameterGrid(method.parameter_grid))
            cv_rows = [
                {**parameters, "candidate_order": candidate_index}
                for candidate_index, parameters in enumerate(parameter_candidates)
            ]
            for fold_index, (study, _, x_fit, y_fit, weights, x_valid, y_valid) in enumerate(
                fold_data, start=1
            ):
                LOGGER.info(
                    "Fitting CV fold %d/%d with %s held out (%d validation, %d training samples; %d candidate workers)",
                    fold_index,
                    len(fold_data),
                    study,
                    len(y_valid),
                    len(y_fit),
                    min(n_jobs, len(parameter_candidates)),
                )
                candidate_results = joblib.Parallel(n_jobs=n_jobs, prefer="threads")(
                    joblib.delayed(_score_cv_candidate)(
                        method, parameters, x_fit, y_fit, weights, x_valid, y_valid,
                    )
                    for parameters in parameter_candidates
                )
                for candidate_index, (score, diagnostics) in enumerate(candidate_results):
                    cv_rows[candidate_index][f"mae_{study}"] = score
                    cv_rows[candidate_index][f"iteration_limit_reached_{study}"] = diagnostics.get(
                        "iteration_limit_reached", False
                    )
                    cv_rows[candidate_index][f"n_iterations_{study}"] = _diagnostic_maximum(
                        diagnostics, "n_iter"
                    )
                    cv_rows[candidate_index][f"dual_gap_{study}"] = _diagnostic_maximum(
                        diagnostics, "dual_gap"
                    )

            for row in cv_rows:
                scores = [row[f"mae_{study}"] for study in train_study_ids]
                row["mean_mae_years"] = float(np.mean(scores))
                row["sd_mae_years"] = float(np.std(scores, ddof=1)) if len(scores) > 1 else 0.0
                row["n_iteration_limit_folds"] = sum(
                    row[f"iteration_limit_reached_{study}"] for study in train_study_ids
                )
            cv_grid = pd.DataFrame(cv_rows).sort_values(
                ["mean_mae_years", "candidate_order"],
                kind="stable",
            ).reset_index(drop=True)
            best_row = cv_grid.iloc[0]
            best_parameters = parameter_candidates[int(best_row["candidate_order"])]
            cv_best_mae = float(best_row["mean_mae_years"])
            cv_by_study = [
                {
                    "held_out_study": study,
                    "mae_years": float(best_row[f"mae_{study}"]),
                    **fold_coverage[study],
                    "iteration_limit_reached": bool(best_row[f"iteration_limit_reached_{study}"]),
                    "n_iterations": _optional_float(best_row[f"n_iterations_{study}"]),
                    "dual_gap": _optional_float(best_row[f"dual_gap_{study}"]),
                }
                for study in train_study_ids
            ]
            cv_oof_frames = []
            fold_coefficient_coverage = {}
            for fold_index, (study, fold_genes, x_fit, y_fit, weights, x_valid, y_valid) in enumerate(
                fold_data, start=1
            ):
                estimator = _fit_estimator(method, best_parameters, x_fit, y_fit, weights)
                predicted = estimator.predict(x_valid)
                fold_coefficients = (
                    getattr(estimator.named_steps[method.coefficient_step], "coef_", None)
                    if method.coefficient_step else None
                )
                fold_coefficient_coverage[study] = _coefficient_coverage(
                    adata, fold_genes, fold_coefficients, study
                )
                held_out = train_meta.loc[train_meta["study"] == study].copy()
                if len(held_out) != len(predicted):
                    raise RuntimeError(f"CV prediction alignment failed for study {study}")
                held_out["predicted_age"] = predicted
                held_out["age_error_years"] = predicted - y_valid
                held_out["cv_fold"] = fold_index
                for name in ("n_model_genes", "n_available_genes", "fraction_available"):
                    held_out[name] = fold_coverage[study][name]
                cv_oof_frames.append(held_out)
            for record in cv_by_study:
                record.update(fold_coefficient_coverage[record["held_out_study"]])
            cv_oof_predictions = pd.concat(cv_oof_frames).sort_index()
        else:
            cv_by_study = []
            cv_oof_predictions = pd.DataFrame()
            best_parameters = method.default_parameters
            cv_best_mae = None
            cv_grid = pd.DataFrame()

        # Final gene universe uses all development studies; external sites play no role.
        gene_selector = method.gene_selector or _training_gene_names
        final_genes = gene_selector(adata, train_study_ids)
        final_gene_intersection = _training_gene_intersection_diagnostics(
            adata, train_study_ids, final_genes
        )
        x_final = _log_cpm(adata, train_positions, final_genes)
        x_final_raw = x_final
        final_zscore_diagnostics = []
        if method.within_study_zscore:
            x_final, final_zscore_diagnostics = within_study_gene_zscore(
                x_final, train_meta["study"].to_numpy()
            )
        final_weights = _sample_weights(train_meta["study"], method)
        estimator = _fit_estimator(
            method, best_parameters, x_final, y_train, final_weights
        )
        fitted.append({
            "method": method,
            "genes": final_genes,
            "estimator": estimator,
            "best_parameters": best_parameters,
            "cv_best_mae": cv_best_mae,
            "cv_by_study": cv_by_study,
            "cv_grid": cv_grid,
            "cv_oof_predictions": cv_oof_predictions,
            "training_gene_intersection": final_gene_intersection,
            "within_study_zscore_diagnostics": final_zscore_diagnostics,
            "training_features_raw": x_final_raw,
        })

    # Generate external features only after every method has completed tuning and fitting.
    for item in fitted:
        method = item["method"]
        estimator = item["estimator"]
        genes = item["genes"]
        model = estimator.named_steps[method.estimator_step]
        coefficients = (
            getattr(estimator.named_steps[method.coefficient_step], "coef_", None)
            if method.coefficient_step else None
        )
        x_external_raw = _log_cpm(adata, external_positions, genes)
        x_external = x_external_raw
        external_zscore_diagnostics = []
        if method.within_study_zscore:
            x_external, external_zscore_diagnostics = within_study_gene_zscore(
                x_external, external_meta["study"].to_numpy()
            )
        external_coverage = _gene_coverage(adata, external_positions, genes)
        for coverage in external_coverage:
            coverage.update(_coefficient_coverage(
                adata, genes, coefficients, coverage["study"]
            ))
            if coverage["missing_coefficient_weight_warning"]:
                LOGGER.warning(
                    "Study %s is missing %.1f%% of absolute coefficient weight; "
                    "predictions may be unreliable",
                    coverage["study"],
                    100 * coverage["missing_coefficient_weight_fraction"],
                )
        predicted = estimator.predict(x_external)
        calibration_splits = pd.DataFrame()
        calibration_summary = {}
        if calibration_samples:
            if method.within_study_zscore:
                calibration_summary = {
                    study: {
                        "status": "not_compatible_with_whole_cohort_zscore",
                        "n_external_samples": int(
                            (external_meta["study"].astype(str) == study).sum()
                        ),
                    }
                    for study in sorted(external_meta["study"].astype(str).unique())
                }
            else:
                calibration_splits, calibration_summary = _target_expression_calibration(
                    item["training_features_raw"],
                    y_train,
                    train_meta["study"].to_numpy(dtype=str),
                    x_external_raw,
                    external_meta["age"].to_numpy(dtype=float),
                    external_meta["study"].to_numpy(dtype=str),
                    np.isin(genes, np.asarray(adata.var_names, dtype=str)),
                    calibration_samples=calibration_samples,
                    repeats=calibration_repeats,
                    seed=calibration_seed,
                )
                comparison_rows = []
                if not calibration_splits.empty:
                    for split in calibration_splits.itertuples(index=False):
                        evaluation_indices = split.evaluation_indices
                        y_evaluation = external_meta.iloc[evaluation_indices][
                            "age"
                        ].to_numpy(dtype=float)
                        raw_evaluation = estimator.predict(
                            x_external_raw[evaluation_indices]
                        )
                        evaluation_features = _dense_features(
                            x_external_raw[evaluation_indices]
                        )
                        corrected_evaluation = estimator.predict(
                            evaluation_features - split.site_shift[None, :]
                        )
                        for calibration_name, values in (
                            ("none", raw_evaluation),
                            ("age_protected_mean_shift", corrected_evaluation),
                        ):
                            comparison_rows.append({
                                "study": split.study,
                                "repeat": split.repeat,
                                "calibration": calibration_name,
                                "n_calibration_samples": split.n_calibration_samples,
                                "n_evaluation_samples": len(y_evaluation),
                                "mae_years": float(mean_absolute_error(y_evaluation, values)),
                                "rmse_years": float(np.sqrt(mean_squared_error(y_evaluation, values))),
                                "mean_error_years": float(np.mean(values - y_evaluation)),
                                "calibration_age_min": split.calibration_age_min,
                                "calibration_age_max": split.calibration_age_max,
                                "evaluation_age_min": split.evaluation_age_min,
                                "evaluation_age_max": split.evaluation_age_max,
                                "mean_absolute_expression_shift": split.mean_absolute_expression_shift,
                            })
                calibration_splits = pd.DataFrame(comparison_rows)
                for study, summary in calibration_summary.items():
                    if summary.get("status") != "completed":
                        continue
                    study_rows = calibration_splits.loc[
                        calibration_splits["study"] == study
                    ]
                    summary["metrics"] = {
                        name: {
                            "mean_mae_years": float(group["mae_years"].mean()),
                            "sd_mae_across_repeats_years": float(
                                group["mae_years"].std(ddof=1)
                                if len(group) > 1 else 0.0
                            ),
                            "mean_rmse_years": float(group["rmse_years"].mean()),
                            "mean_bias_years": float(group["mean_error_years"].mean()),
                        }
                        for name, group in study_rows.groupby("calibration", sort=True)
                    }
        predictions = external_meta.copy()
        coverage_by_study = {record["study"]: record for record in external_coverage}
        for name in (
            "n_model_genes", "n_available_genes", "fraction_available",
            "missing_coefficient_weight_fraction", "missing_coefficient_weight_warning",
        ):
            predictions[name] = predictions["study"].map(
                {study: record[name] for study, record in coverage_by_study.items()}
            )
        predictions["predicted_age"] = predicted
        predictions["age_error_years"] = predicted - predictions["age"].to_numpy(float)
        metrics_by_study = {
            str(study): _metrics(group["age"], group["predicted_age"].to_numpy(float))
            for study, group in predictions.groupby("study", observed=True, sort=True)
        }
        training_age_coverage = [
            {
                "study": str(study),
                "n_samples": len(group),
                "minimum_age": float(group["age"].min()),
                "median_age": float(group["age"].median()),
                "maximum_age": float(group["age"].max()),
            }
            for study, group in train_meta.groupby("study", observed=True, sort=True)
        ]
        result = {
            "method": method.name,
            "method_name": method.display_name,
            "method_description": method.description,
            "feature_representation": method.feature_representation,
            "training_weighting": "equal total weight per study" if method.equal_study_weight else "equal weight per sample",
            "within_study_zscore": method.within_study_zscore,
            "cell_type": selected_cell_type,
            "n_genes": len(genes),
            "n_training_samples": len(train_meta),
            "n_external_samples": len(predictions),
            "n_nonzero_coefficients": (
                int(np.count_nonzero(coefficients)) if coefficients is not None else None
            ),
            "training_studies_present": train_study_ids,
            "training_gene_intersection": item["training_gene_intersection"],
            "within_study_zscore_diagnostics": item["within_study_zscore_diagnostics"],
            "external_within_study_zscore_diagnostics": external_zscore_diagnostics,
            "training_age_coverage_by_study": training_age_coverage,
            "best_parameters": item["best_parameters"],
            "best_training_cv_mae_years": item["cv_best_mae"],
            "cv_mae_by_held_out_study": item["cv_by_study"],
            "fit_diagnostics": _fit_diagnostics(estimator, method),
            "external_gene_coverage_by_study": external_coverage,
            "external_metrics_by_study": metrics_by_study,
            "external_metrics_pooled": _metrics(predictions["age"], predicted),
            "external_metrics_macro_study_mae_years": _mean_study_mae(predictions),
            "target_expression_calibration": calibration_summary,
        }
        item.update({
            "predictions": predictions,
            "calibration_splits": calibration_splits,
            "coefficients": coefficients,
            "model": model,
            "result": result,
        })
        manifest["methods"].append(result)

    return ClockFit(manifest, inclusion_counts, fitted)


def publish_age_clock_fit(fit: ClockFit, output_dir: Path) -> dict:
    """Publish a completed in-memory fit; publication failures remain fatal."""
    manifest = dict(fit.manifest)
    cell_type = str(manifest["cell_type"])
    output_type_dir = output_dir / slugify(cell_type) / "ageing_clock"
    output_type_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_type_dir / "age_clock.json"
    newly_fitted_methods = {item["method"].name for item in fit.fitted_methods}
    manifest["newly_fitted_methods"] = sorted(newly_fitted_methods)
    if manifest_path.is_file():
        previous = json.loads(manifest_path.read_text())
        identity_fields = (
            "input", "configuration_sha256", "cell_type",
            "training_studies", "external_studies",
            "target_expression_calibration",
        )
        compatible = all(previous.get(key) == manifest.get(key) for key in identity_fields)
        if compatible:
            previous_methods = {
                item["method"]: item for item in previous.get("methods", [])
                if item.get("method") not in newly_fitted_methods
            }
            current_methods = {item["method"]: item for item in manifest["methods"]}
            merged_methods = {**previous_methods, **current_methods}
            manifest["methods"] = list(merged_methods.values())
            manifest["reused_prior_methods"] = sorted(previous_methods)
    run_metadata = {key: value for key, value in manifest.items() if key != "methods"}

    for item in fit.fitted_methods:
        method = item["method"]
        method_dir = output_type_dir / method.name
        method_dir.mkdir(parents=True, exist_ok=True)
        if item["cv_best_mae"] is not None:
            item["cv_grid"].to_csv(method_dir / "cv_grid_search.csv", index=False)
            pd.DataFrame(item["cv_by_study"]).to_csv(method_dir / "cv_by_study.csv", index=False)
            item["cv_oof_predictions"].to_csv(
                method_dir / "cv_oof_predictions.csv", index=True, index_label="pseudobulk_id"
            )
        else:
            for filename in ("cv_grid_search.csv", "cv_by_study.csv", "cv_oof_predictions.csv"):
                (method_dir / filename).unlink(missing_ok=True)

        genes = item["genes"]
        estimator = item["estimator"]
        model = item["model"]
        result = item["result"]
        result["model"] = str(method_dir / "model.joblib")
        result["predictions"] = str(method_dir / "external_predictions.csv")
        item["predictions"].to_csv(method_dir / "external_predictions.csv", index=True)
        calibration_splits = item.get("calibration_splits", pd.DataFrame())
        calibration_path = method_dir / "external_expression_calibration_repeats.csv"
        if not calibration_splits.empty:
            calibration_splits.to_csv(calibration_path, index=False)
        else:
            calibration_path.unlink(missing_ok=True)
        coefficients = item["coefficients"]
        if method.coefficient_step and coefficients is not None:
            pd.DataFrame({
                "gene": genes,
                "coefficient": coefficients,
                "absolute_coefficient": np.abs(coefficients),
            }).sort_values("absolute_coefficient", ascending=False).to_csv(
                method_dir / "coefficients.csv", index=False
            )
        else:
            (method_dir / "coefficients.csv").unlink(missing_ok=True)
        bundle = {
            "model_format_version": MODEL_FORMAT_VERSION,
            "estimator": estimator,
            "genes": genes,
            "transform": "counts_per_million_then_log1p",
            "scaling": method.scaling_description,
            "cell_type": cell_type,
            "method": method.name,
            "training_weighting": result["training_weighting"],
            "target": "age_years",
            "normalization": manifest["normalization"],
            "within_study_zscore": method.within_study_zscore,
            "missing_gene_policy": manifest["missing_gene_policy"],
            "configuration_sha256": manifest["configuration_sha256"],
            "software_versions": manifest["software_versions"],
        }
        joblib.dump(bundle, method_dir / "model.joblib")
        method_metadata = {
            **run_metadata,
            **result,
            "genes": genes,
            "sample_inclusion_counts": fit.inclusion_counts,
            "model_intercept_years": (
                float(model.intercept_) if hasattr(model, "intercept_") else None
            ),
            "method_dir": str(method_dir),
        }
        (method_dir / "metadata.json").write_text(
            json.dumps(method_metadata, indent=2) + "\n"
        )
        method_manifest = {
            **run_metadata,
            "estimator": type(model).__name__,
            "cell_types": [result],
        }
        (method_dir / "age_clock.json").write_text(
            json.dumps(method_manifest, indent=2) + "\n"
        )

    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def train_age_clocks(
    input_path: Path,
    config_path: Path,
    output_dir: Path,
    external_studies: list[str],
    *,
    cell_type_filter: str | None = None,
    alphas: list[float] | None = None,
    l1_ratios: list[float] | None = None,
    method_names: list[str] | None = None,
    calibration_samples: int = 10,
    calibration_repeats: int = 50,
    calibration_seed: int = 42,
    cv: bool = True,
    n_jobs: int = 1,
) -> dict:
    """Fit and publish configured clock methods for one cell type."""
    fit = fit_age_clocks(
        input_path, config_path, external_studies,
        cell_type_filter=cell_type_filter, alphas=alphas, l1_ratios=l1_ratios,
        method_names=method_names, calibration_samples=calibration_samples,
        calibration_repeats=calibration_repeats, calibration_seed=calibration_seed,
        cv=cv,
        n_jobs=n_jobs,
    )
    return publish_age_clock_fit(fit, output_dir)


def train_all_age_clocks(
    input_path: Path,
    config_path: Path,
    output_dir: Path,
    external_studies: list[str],
    *,
    alphas: list[float] | None = None,
    l1_ratios: list[float] | None = None,
    method_names: list[str] | None = None,
    calibration_samples: int = 10,
    calibration_repeats: int = 50,
    calibration_seed: int = 42,
    cv: bool = True,
    n_jobs: int = 1,
) -> dict[str, dict]:
    """Fit every type, isolating expected and unexpected fit failures per type."""
    import anndata as ad

    pipeline_config = _load_config(config_path)
    split_by = pipeline_config["differential_expression"]["split_by"]
    adata = ad.read_h5ad(input_path, backed="r")
    if split_by not in adata.obs:
        raise ValueError(f"Pseudobulk metadata is missing configured cell-type field {split_by!r}")
    cell_types = sorted(adata.obs[split_by].astype(str).unique())
    adata.file.close()
    if not cell_types:
        raise ValueError("Pseudobulk input contains no cell types")
    results = {}
    summary = {
        "kind": "all_cell_type_age_clock_run",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "external_studies": external_studies,
        "successful_cell_types": [],
        "eligibility_skipped_cell_types": [],
        "failed_cell_types": [],
    }
    for cell_type in cell_types:
        LOGGER.info("Starting all-cell-type age-clock run for %s", cell_type)
        try:
            fit = fit_age_clocks(
                input_path, config_path, external_studies,
                cell_type_filter=cell_type, alphas=alphas, l1_ratios=l1_ratios,
                method_names=method_names, calibration_samples=calibration_samples,
                calibration_repeats=calibration_repeats,
                calibration_seed=calibration_seed, cv=cv,
                n_jobs=n_jobs,
            )
        except ClockEligibilityError as error:
            LOGGER.error("Skipping age clock for %s: %s", cell_type, error)
            summary["eligibility_skipped_cell_types"].append({
                "cell_type": cell_type,
                "exception": type(error).__name__,
                "message": str(error),
            })
            continue
        except Exception as error:
            formatted_traceback = traceback.format_exc()
            LOGGER.exception("Age-clock fit failed for %s; continuing with other cell types", cell_type)
            summary["failed_cell_types"].append({
                "cell_type": cell_type,
                "exception": type(error).__name__,
                "message": str(error),
                "traceback": formatted_traceback,
            })
            continue
        results[cell_type] = publish_age_clock_fit(fit, output_dir)
        summary["successful_cell_types"].append(cell_type)

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "age_clock_all_types_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    if summary["failed_cell_types"]:
        raise ClockFitFailuresError(
            f"Age-clock fitting failed for {len(summary['failed_cell_types'])} cell type(s); "
            f"successful types were published and details are in {summary_path}"
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train continuous-age elastic-net clocks from pseudobulk H5AD"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--external-studies", required=True,
        help="Comma-separated study IDs held out from preprocessing, tuning, and fitting",
    )
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--cell-type", help="AIFI L2 cell type label or slug to fit")
    selection.add_argument(
        "--all-cell-types", action="store_true",
        help="Fit clocks for every cell type present in the pseudobulk input",
    )
    parser.add_argument("--alphas", default="0.001,0.01,0.1,1,10")
    parser.add_argument("--l1-ratios", default="0.1,0.5,0.9,1")
    parser.add_argument(
        "--methods",
        help="Comma-separated method IDs (default: all methods)",
    )
    parser.add_argument(
        "--calibration-samples",
        type=int,
        default=10,
        help="Known-age target donors per study used for expression calibration; 0 disables it",
    )
    parser.add_argument(
        "--calibration-repeats",
        type=int,
        default=50,
        help="Repeated calibration/evaluation splits per target study (default: 50)",
    )
    parser.add_argument(
        "--calibration-seed",
        type=int,
        default=42,
        help="Random seed for target calibration donor splits (default: 42)",
    )
    parser.add_argument(
        "--n-jobs", type=int, default=1,
        help="parallel CV candidate fits per held-out-study fold (default: 1)",
    )
    parser.add_argument("--no-cv", action="store_true", help="Fit fixed alpha=0.1, l1_ratio=0.5")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fit_args = {
        "alphas": [float(value) for value in args.alphas.split(",")],
        "l1_ratios": [float(value) for value in args.l1_ratios.split(",")],
        "method_names": (
            [value.strip() for value in args.methods.split(",") if value.strip()]
            if args.methods else None
        ),
        "calibration_samples": args.calibration_samples,
        "calibration_repeats": args.calibration_repeats,
        "calibration_seed": args.calibration_seed,
        "cv": not args.no_cv,
        "n_jobs": args.n_jobs,
    }
    external_studies = [value.strip() for value in args.external_studies.split(",") if value.strip()]
    if args.all_cell_types:
        train_all_age_clocks(
            args.input, args.config, args.output_dir, external_studies, **fit_args
        )
    else:
        train_age_clocks(
            args.input, args.config, args.output_dir, external_studies,
            cell_type_filter=args.cell_type, **fit_args,
        )


if __name__ == "__main__":
    main()

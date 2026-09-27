"""Run sample-level age-associated differential expression on pseudobulk counts."""
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

from .config import config_digest, read_json

LOGGER = logging.getLogger(__name__)
_UNKNOWN = {"", "not_provided", "heterogeneous during bulk", "nan", "none"}
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
    metadata = obs[["study", "sample", "age", "sex", "n_cells"]].copy()
    metadata["study"] = metadata["study"].astype("string").str.strip()
    metadata["sample"] = metadata["sample"].astype("string").str.strip()
    metadata["sex"] = metadata["sex"].astype("string").str.strip().str.lower()
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce")
    metadata["n_cells"] = pd.to_numeric(metadata["n_cells"], errors="coerce")
    minimum_age = inclusion["minimum_age_years_inclusive"]
    minimum_cells = inclusion["minimum_cells_per_pseudobulk_inclusive"]
    valid_study = metadata["study"].notna() & ~metadata["study"].str.lower().isin(_UNKNOWN)
    valid_sample = metadata["sample"].notna() & ~metadata["sample"].str.lower().isin(_UNKNOWN)
    valid_age = metadata["age"].notna() & np.isfinite(metadata["age"])
    adult = valid_age & (metadata["age"] >= minimum_age)
    valid_sex = metadata["sex"].notna() & ~metadata["sex"].str.lower().isin(_UNKNOWN)
    valid_cells = metadata["n_cells"].notna() & np.isfinite(metadata["n_cells"])
    enough_cells = valid_cells & (metadata["n_cells"] >= minimum_cells)
    valid_other_metadata = valid_study & valid_sample & adult & valid_sex & valid_cells
    valid = valid_other_metadata & (enough_cells | test_mode)
    excluded = []
    for position in np.flatnonzero(~valid.to_numpy(dtype=bool)):
        row = metadata.iloc[position]
        reasons = []
        if not valid_study.iloc[position]:
            reasons.append("missing_or_unknown_study")
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
        "missing_or_unknown_sample": int((~valid_sample).sum()),
        "missing_or_non_numeric_age": int((~valid_age).sum()),
        "age_below_minimum": int((valid_age & ~adult).sum()),
        "missing_or_unknown_sex": int((~valid_sex).sum()),
        "missing_or_non_numeric_n_cells": int((~valid_cells).sum()),
        "n_cells_below_minimum": int((valid_cells & ~enough_cells).sum()),
        "n_cells_filter_bypassed": int(
            (valid_other_metadata & ~enough_cells).sum() if test_mode else 0
        ),
        "eligible_before_duplicate_check": int(valid.sum()),
    }
    metadata = metadata.loc[valid, ["study", "sample", "age", "sex", "n_cells"]].copy()
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


def _fit_age_model(
    counts: np.ndarray,
    genes: pd.Index,
    metadata: pd.DataFrame,
    *,
    design: str,
    age_term: str,
    alpha: float,
    cpus: int,
) -> pd.DataFrame:
    from formulaic_contrasts import FormulaicContrasts
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_df = pd.DataFrame(counts, index=metadata.index, columns=genes)
    model_metadata = metadata.copy()
    for column in ("study", "sex"):
        if column in design:
            model_metadata[column] = pd.Categorical(
                model_metadata[column], categories=sorted(model_metadata[column].unique())
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
    if age_term not in design_matrix.columns:
        raise ValueError(
            f"Age coefficient {age_term!r} not found in design matrix columns "
            f"{list(design_matrix.columns)}"
        )
    contrast = np.zeros(design_matrix.shape[1], dtype=float)
    contrast[list(design_matrix.columns).index(age_term)] = 1.0
    stats = DeseqStats(
        dds,
        contrast=contrast,
        alpha=alpha,
        n_cpus=cpus,
        quiet=True,
    )
    stats.summary()
    return stats.results_df


def _run_model(
    *,
    adata,
    obs: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_type: str,
    model_name: str,
    study: str | None,
    design: str,
    age_term: str,
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
        "contrast": f"log2 fold change per one-year increase in {age_term}",
        "n_samples": len(metadata),
        "status": "skipped",
    }
    if metadata.empty:
        details["reason"] = (
            "no samples pass the configured age, sex, and minimum-cell inclusion filters"
        )
        return details
    if metadata["age"].nunique() < 2:
        details["reason"] = "age has fewer than two distinct values"
        return details
    if metadata["sex"].nunique() < 2:
        details["reason"] = "sex has fewer than two observed levels"
        return details
    if "study" in design and metadata["study"].nunique() < 2:
        details["reason"] = "study has fewer than two observed levels"
        return details

    studies = sorted(metadata["study"].astype(str).unique())
    var_mask = _availability_mask(adata.var, studies, context=f"{cell_type}/{model_name}")
    genes = adata.var_names[var_mask]
    if len(genes) == 0:
        details["reason"] = "no genes are available in every study included in this fit"
        return details

    # The order of pseudobulk rows in the AnnData object matches obs/metadata.
    obs_indices = metadata.index
    row_positions = obs.index.get_indexer(obs_indices)
    col_positions = adata.var_names.get_indexer(genes)
    count_values = _matrix_to_counts(
        adata.X[row_positions, :][:, col_positions], context=f"{cell_type}/{model_name}"
    )
    nonzero_gene = count_values.sum(axis=0) > 0
    details["n_all_zero_genes_excluded"] = int((~nonzero_gene).sum())
    count_values = count_values[:, nonzero_gene]
    genes = genes[nonzero_gene]
    if len(genes) == 0:
        details["reason"] = "all genes in the analysis universe have zero counts"
        return details
    nonzero_library = count_values.sum(axis=1) > 0
    if not nonzero_library.all():
        zero_library_metadata = metadata.iloc[np.flatnonzero(~nonzero_library)]
        metadata = metadata.iloc[np.flatnonzero(nonzero_library)].copy()
        count_values = count_values[nonzero_library]
        details["n_zero_library_samples_excluded"] = int((~nonzero_library).sum())
        details["zero_library_samples_excluded"] = [
            {
                "pseudobulk_id": str(index),
                "study": str(row["study"]),
                "sample": str(row["sample"]),
            }
            for index, row in zero_library_metadata.iterrows()
        ]
    try:
        results = _fit_age_model(
            count_values,
            genes,
            metadata,
            design=design,
            age_term=age_term,
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

    result = results.rename_axis("gene").reset_index()
    result.insert(0, "cell_type", cell_type)
    result.insert(1, "model", model_name)
    result.insert(2, "study", study if study is not None else "all")
    test_only = test_mode or synthetic_test_data
    result.insert(3, "analysis_mode", "test_only" if test_only else "standard")
    result.insert(
        4,
        "interpretation_warning",
        _SYNTHETIC_DATA_WARNING if synthetic_test_data else (
            _TEST_MODE_WARNING if test_mode else ""
        ),
    )
    cell_type_dir = output_dir / _slug(cell_type)
    if model_name == "per_study":
        destination = cell_type_dir / "per_study" / f"{_slug(study or 'unknown-study')}.csv"
    else:
        destination = cell_type_dir / "merged.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, index=False)
    details.update(
        {
            "status": "complete",
            "n_genes": len(genes),
            "results": str(destination.relative_to(output_dir)),
            "n_significant_adjusted_p_lt_alpha": int((results["padj"] < alpha).sum()),
        }
    )
    return details


def run_differential_expression(
    input_path: Path,
    output_dir: Path,
    pipeline: dict,
    *,
    cpus: int = 1,
    test_mode: bool = False,
    synthetic_test_data: bool = False,
) -> dict:
    """Fit per-study and shared-slope PyDESeq2 models by AIFI L2 type."""
    import anndata as ad

    spec = pipeline["differential_expression"]
    if spec["method"] != "pydeseq2":
        raise ValueError("Only the configured PyDESeq2 method is supported")
    split_by = spec["split_by"]
    adata = ad.read_h5ad(input_path)
    required = {"study", "sample", "age", "sex", "n_cells", split_by}
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

    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    used_slugs: dict[str, str] = {}
    for cell_type in sorted(labels.unique().tolist()):
        slug = _slug(str(cell_type))
        if slug in used_slugs and used_slugs[slug] != cell_type:
            raise ValueError(
                f"AIFI L2 labels collide as output names: {used_slugs[slug]!r} and {cell_type!r}"
            )
        used_slugs[slug] = str(cell_type)
        selected = labels == cell_type
        type_adata = adata[selected.to_numpy(), :]
        type_obs = type_adata.obs.copy()
        (output_dir / slug).mkdir(parents=True, exist_ok=True)
        metadata, exclusions, excluded_samples = _prepare_metadata(
            type_obs, spec["sample_inclusion"], test_mode=test_mode
        )
        item = {
            "cell_type": str(cell_type),
            "slug": slug,
            "n_pseudobulks": int(type_adata.n_obs),
            "metadata_exclusions": exclusions,
            "excluded_samples": excluded_samples,
            "models": [],
        }
        for study in sorted(metadata["study"].astype(str).unique()):
            study_metadata = metadata.loc[metadata["study"].astype(str) == study].copy()
            model_metadata = study_metadata[["study", "sample", "age", "sex"]]
            item["models"].append(
                _run_model(
                    adata=type_adata,
                    obs=type_obs,
                    metadata=model_metadata,
                    cell_type=str(cell_type),
                    model_name="per_study",
                    study=study,
                    design=spec["per_study_design"],
                    age_term=spec["age_term"],
                    alpha=spec["alpha"],
                    cpus=cpus,
                    output_dir=output_dir,
                    test_mode=test_mode,
                    synthetic_test_data=synthetic_test_data,
                )
            )
        item["models"].append(
            _run_model(
                adata=type_adata,
                obs=type_obs,
                metadata=metadata[["study", "sample", "age", "sex"]],
                cell_type=str(cell_type),
                model_name="merged",
                study=None,
                design=spec["merged_design"],
                age_term=spec["age_term"],
                alpha=spec["alpha"],
                cpus=cpus,
                output_dir=output_dir,
                test_mode=test_mode,
                synthetic_test_data=synthetic_test_data,
            )
        )
        type_metadata = {
            "kind": "pseudobulk_differential_expression_cell_type_metadata",
            "cell_type": str(cell_type),
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
            "metadata_exclusions": exclusions,
        }
        (output_dir / slug / "run_metadata.json").write_text(
            json.dumps(type_metadata, indent=2) + "\n"
        )
        records.append(item)

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
        "sample_inclusion": spec["sample_inclusion"],
        "sample_inclusion_applied": {
            **spec["sample_inclusion"],
            "minimum_cells_filter_enabled": not test_mode,
        },
        "models": {
            "per_study": spec["per_study_design"],
            "merged": spec["merged_design"],
        },
        "age_effect": "log2 fold change per one-year increase in age",
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
                    {"study": model["study"], "path": model["results"]}
                    for model in item["models"]
                    if model["model"] == "per_study" and model.get("status") == "complete"
                ],
                "merged": next(
                    (
                        model["results"] for model in item["models"]
                        if model["model"] == "merged" and model.get("status") == "complete"
                    ),
                    None,
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run per-study and merged PyDESeq2 models on a merged pseudobulk H5AD"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=1)
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
    run_differential_expression(
        args.input,
        args.output_dir,
        read_json(args.config),
        cpus=args.cpus,
        test_mode=args.test_mode,
        synthetic_test_data=args.synthetic_test_data,
    )
    LOGGER.info("Differential-expression results written to %s", args.output_dir)


if __name__ == "__main__":
    main()

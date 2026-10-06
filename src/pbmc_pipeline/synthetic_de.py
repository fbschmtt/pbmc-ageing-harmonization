"""Create a compact, deterministic pseudobulk fixture for positive DE tests."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .config import AgeTrajectorySettings, read_json

TRAJECTORY_PROFILE_MULTIPLIERS = [
    [1.0, 1.8, 0.55, 2.6, 0.7, 2.0, 0.65],
    [1.0, 2.0, 2.0, 0.5, 0.5, 2.4, 0.75],
    [1.0, 0.6, 1.8, 2.2, 0.5, 2.0, 1.7],
    [1.0, 2.2, 0.6, 0.6, 2.5, 1.0, 2.1],
    [1.0, 0.5, 2.2, 0.5, 1.7, 2.5, 0.55],
]


def _default_samples_per_study(settings: AgeTrajectorySettings, study_count: int) -> int:
    """Derive balanced synthetic replication from the active trajectory support rule."""
    age_bin_count = (
        settings.strict_age_cutoff_exclusive - settings.reference_bin_start_age
    ) // settings.bin_width_years
    if age_bin_count != len(TRAJECTORY_PROFILE_MULTIPLIERS[0]):
        raise ValueError(
            f"The synthetic age-trajectory fixture plants {len(TRAJECTORY_PROFILE_MULTIPLIERS[0])} "
            "bins; pipeline settings must produce that many bins below the strict age cutoff"
        )
    if settings.minimum_bins > age_bin_count:
        raise ValueError(
            "The synthetic age-trajectory fixture cannot meet minimum_bins with its "
            f"{age_bin_count} available bins"
        )
    per_study_per_bin = int(np.ceil(settings.minimum_samples_per_bin / study_count))
    return age_bin_count * per_study_per_bin


def create_synthetic_pseudobulk(
    output_path: Path,
    *,
    trajectory_settings: AgeTrajectorySettings,
    seed: int = 413,
    samples_per_study: int | None = None,
    cells_per_pseudobulk: int = 15,
) -> Path:
    """Write a small pseudobulk H5AD with estimable, deliberately synthetic fits.

    Keep the planted non-linear profiles aligned with the production age-bin model;
    ``check_synthetic_de.py`` verifies that the fitted coefficients recover them.
    """
    studies = ["synthetic_study_a", "synthetic_study_b", "synthetic_study_c"]
    if samples_per_study is None:
        samples_per_study = _default_samples_per_study(trajectory_settings, len(studies))
    if samples_per_study < 6:
        raise ValueError("samples_per_study must be at least 6 to exercise both models")
    if cells_per_pseudobulk < 10:
        raise ValueError("cells_per_pseudobulk must meet the configured DE cell cutoff (10)")

    rng = np.random.default_rng(seed)
    cell_types = ["CD14 monocyte", "Naive CD4 T cell"]
    age_bin_starts = list(range(
        trajectory_settings.reference_bin_start_age,
        trajectory_settings.strict_age_cutoff_exclusive,
        trajectory_settings.bin_width_years,
    ))
    if (
        len(age_bin_starts) != len(TRAJECTORY_PROFILE_MULTIPLIERS[0])
        or trajectory_settings.minimum_bins > len(age_bin_starts)
    ):
        raise ValueError(
            "The synthetic age-trajectory fixture's planted bins must cover "
            "the configured minimum-bin threshold"
        )
    per_bin_sample_count = samples_per_study // len(age_bin_starts)
    if samples_per_study % len(age_bin_starts):
        raise ValueError(
            "samples_per_study must divide evenly across the configured synthetic age bins"
        )
    if per_bin_sample_count * len(studies) < trajectory_settings.minimum_samples_per_bin:
        raise ValueError(
            "samples_per_study is too small to meet minimum_samples_per_bin after "
            "combining the synthetic studies"
        )
    ages = np.array([
        start + offset
        for start in age_bin_starts
        for offset in np.linspace(
            trajectory_settings.bin_width_years * 0.1,
            trajectory_settings.bin_width_years * 0.9,
            per_bin_sample_count,
        )
    ], dtype=float)
    # Balance sex within each decade without aligning it with the alternating CMV labels.
    sex_pattern = np.resize(
        np.array(["female", "female", "male", "male"]), per_bin_sample_count
    )
    sexes = np.tile(sex_pattern, len(age_bin_starts))
    bmi_values = rng.permutation(np.linspace(21.0, 34.0, samples_per_study))
    cmv_values = np.resize(np.array(["negative", "positive", "negative", "positive"]), samples_per_study)

    base_genes = [f"SYNTH_GENE_{index:03d}" for index in range(80)]
    age_genes = ["SYNTH_AGE_MARKER_1", "SYNTH_AGE_MARKER_2", "SYNTH_AGE_MARKER_3"]
    trajectory_genes = [f"SYNTH_AGE_TRAJECTORY_MARKER_{index}" for index in range(1, 6)]
    if len(trajectory_genes) != len(TRAJECTORY_PROFILE_MULTIPLIERS):
        raise ValueError("Synthetic age-trajectory gene and profile counts must match")
    sex_genes = ["SYNTH_SEX_MARKER_1", "SYNTH_SEX_MARKER_2"]
    bmi_genes = ["SYNTH_BMI_MARKER_1", "SYNTH_BMI_MARKER_2"]
    cmv_genes = ["SYNTH_CMV_MARKER_1", "SYNTH_CMV_MARKER_2"]
    study_specific_genes = ["SYNTH_STUDY_A_ONLY", "SYNTH_STUDY_B_ONLY"]
    genes = (
        base_genes + age_genes + trajectory_genes + sex_genes + bmi_genes
        + cmv_genes + study_specific_genes
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    for study in studies:
        var[f"available_in_{study}"] = True
    var.loc["SYNTH_STUDY_A_ONLY", [f"available_in_{study}" for study in studies]] = [True, False, False]
    var.loc["SYNTH_STUDY_B_ONLY", [f"available_in_{study}" for study in studies]] = [False, True, False]

    rows: list[dict] = []
    count_rows: list[np.ndarray] = []
    baseline = np.exp(rng.normal(np.log(120), 0.45, size=len(genes)))
    alpha = 0.08  # negative-binomial dispersion, independent of the analysis alpha
    for study_index, study in enumerate(studies):
        for cell_type_index, cell_type in enumerate(cell_types):
            for sample_index, (age, sex) in enumerate(zip(ages, sexes, strict=True)):
                sample = f"SYN{sample_index + 1:02d}"
                site_index = int(
                    study == studies[0] and sample_index % 4 in {1, 2}
                )
                study_site = f"{study}_site_{site_index + 1}"
                bmi = float(bmi_values[sample_index]) if study in {studies[0], studies[2]} else np.nan
                # Leave one additional missing value within each supported study to
                # exercise complete-case selection, while retaining estimable fits.
                if study in {studies[0], studies[2]} and sample_index == 1:
                    bmi = np.nan
                cmv = str(cmv_values[sample_index]) if study == studies[0] else pd.NA
                if study == studies[0] and sample_index == 1:
                    cmv = pd.NA
                mean = baseline.copy() * (1.0 + 0.08 * cell_type_index)
                age_effect = np.array([gene in age_genes for gene in genes])
                sex_effect = np.array([gene in sex_genes for gene in genes])
                bmi_effect = np.array([gene in bmi_genes for gene in genes])
                cmv_effect = np.array([gene in cmv_genes for gene in genes])
                mean[age_effect] *= np.exp(0.025 * (age - 50))
                age_bin_start = (
                    int(age // trajectory_settings.bin_width_years)
                    * trajectory_settings.bin_width_years
                )
                bin_index = (
                    age_bin_start - trajectory_settings.reference_bin_start_age
                ) // trajectory_settings.bin_width_years
                for gene, profile in zip(
                    trajectory_genes, TRAJECTORY_PROFILE_MULTIPLIERS, strict=True
                ):
                    mean[genes.index(gene)] *= profile[bin_index]
                mean[sex_effect & (sex == "male")] *= 1.7
                if np.isfinite(bmi):
                    mean[bmi_effect] *= np.exp(0.15 * (bmi - 27.5))
                if not pd.isna(cmv) and cmv == "positive":
                    mean[cmv_effect] *= 5.0
                if study == "synthetic_study_a":
                    mean[genes.index("SYNTH_STUDY_A_ONLY")] *= 1.2
                elif study == "synthetic_study_b":
                    mean[genes.index("SYNTH_STUDY_B_ONLY")] *= 1.2
                unavailable = ~var[f"available_in_{study}"].to_numpy(dtype=bool)
                mean[unavailable] = 0
                size = 1.0 / alpha
                probability = size / (size + mean)
                counts = rng.negative_binomial(size, probability).astype(np.int32)
                rows.append(
                    {
                        "study": study,
                        "study_site": study_site,
                        "sample": sample,
                        "age": float(age),
                        "sex": sex,
                        "bmi": bmi,
                        "cmv": cmv,
                        "n_cells": int(cells_per_pseudobulk),
                        "total_counts": int(counts.sum()),
                        "aifi_l2_majority": cell_type,
                    }
                )
                count_rows.append(counts)

    obs = pd.DataFrame(rows, index=pd.Index(
        [f"{row['study']}::{row['sample']}::{row['aifi_l2_majority']}" for row in rows],
        name="pseudobulk_id",
    ))
    matrix = np.vstack(count_rows)
    adata = ad.AnnData(X=sparse.csr_matrix(matrix), obs=obs, var=var)
    adata.uns["synthetic_test_data"] = {
        "seed": int(seed),
        "generator": "pbmc_pipeline.synthetic_de.create_synthetic_pseudobulk",
        "not_biological_evidence": True,
        "expected_age_associated_genes": age_genes,
        "expected_age_trajectory_associated_genes": trajectory_genes,
        "expected_age_trajectory_profiles": {
            gene: profile for gene, profile in zip(
                trajectory_genes, TRAJECTORY_PROFILE_MULTIPLIERS, strict=True
            )
        },
        "expected_sex_associated_genes": sex_genes,
        "expected_bmi_associated_genes": bmi_genes,
        "expected_cmv_associated_genes": cmv_genes,
        "expected_covariate_studies": {
            "bmi": [studies[0], studies[2]],
            "cmv": [studies[0]],
        },
        "samples_per_study": int(samples_per_study),
        "age_trajectory_settings": trajectory_settings.to_mapping(),
        "cells_per_pseudobulk": int(cells_per_pseudobulk),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=413)
    parser.add_argument("--samples-per-study", type=int)
    parser.add_argument("--cells-per-pseudobulk", type=int, default=15)
    args = parser.parse_args()
    pipeline = read_json(args.config)
    trajectory_settings = AgeTrajectorySettings.from_mapping(
        pipeline["differential_expression"]["age_trajectory"]
    )
    output = create_synthetic_pseudobulk(
        args.output,
        trajectory_settings=trajectory_settings,
        seed=args.seed,
        samples_per_study=args.samples_per_study,
        cells_per_pseudobulk=args.cells_per_pseudobulk,
    )
    print(output)


if __name__ == "__main__":
    main()

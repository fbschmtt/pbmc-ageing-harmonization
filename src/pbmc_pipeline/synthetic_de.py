"""Create a compact, deterministic pseudobulk fixture for positive DE tests."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def create_synthetic_pseudobulk(
    output_path: Path,
    *,
    seed: int = 413,
    samples_per_study: int = 8,
    cells_per_pseudobulk: int = 15,
) -> Path:
    """Write a small pseudobulk H5AD with estimable, deliberately synthetic fits."""
    if samples_per_study < 6:
        raise ValueError("samples_per_study must be at least 6 to exercise both models")
    if cells_per_pseudobulk < 10:
        raise ValueError("cells_per_pseudobulk must meet the configured DE cell cutoff (10)")

    rng = np.random.default_rng(seed)
    studies = ["synthetic_study_a", "synthetic_study_b", "synthetic_study_c"]
    cell_types = ["CD14 monocyte", "Naive CD4 T cell"]
    ages = np.array([24, 31, 38, 45, 52, 60, 68, 76], dtype=float)
    sexes = np.array(["female", "male", "male", "female", "female", "male", "male", "female"])
    bmi_values = rng.permutation(np.linspace(21.0, 34.0, samples_per_study))
    cmv_values = np.resize(np.array(["negative", "positive", "negative", "positive"]), samples_per_study)
    if samples_per_study != len(ages):
        ages = np.linspace(22, 78, samples_per_study).round().astype(float)
        sexes = np.resize(np.array(["female", "male"]), samples_per_study)
        # Alternate in mirrored blocks to avoid making sex a proxy for age.
        sexes[np.arange(samples_per_study) % 4 >= 2] = np.where(
            sexes[np.arange(samples_per_study) % 4 >= 2] == "female", "male", "female"
        )

    base_genes = [f"SYNTH_GENE_{index:03d}" for index in range(80)]
    age_genes = ["SYNTH_AGE_MARKER_1", "SYNTH_AGE_MARKER_2", "SYNTH_AGE_MARKER_3"]
    sex_genes = ["SYNTH_SEX_MARKER_1", "SYNTH_SEX_MARKER_2"]
    bmi_genes = ["SYNTH_BMI_MARKER_1", "SYNTH_BMI_MARKER_2"]
    cmv_genes = ["SYNTH_CMV_MARKER_1", "SYNTH_CMV_MARKER_2"]
    study_specific_genes = ["SYNTH_STUDY_A_ONLY", "SYNTH_STUDY_B_ONLY"]
    genes = base_genes + age_genes + sex_genes + bmi_genes + cmv_genes + study_specific_genes
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
                    study == studies[0]
                    and (sample_index // 2) % 2 == 1
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
        "expected_sex_associated_genes": sex_genes,
        "expected_bmi_associated_genes": bmi_genes,
        "expected_cmv_associated_genes": cmv_genes,
        "expected_covariate_studies": {
            "bmi": [studies[0], studies[2]],
            "cmv": [studies[0]],
        },
        "samples_per_study": int(samples_per_study),
        "cells_per_pseudobulk": int(cells_per_pseudobulk),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=413)
    parser.add_argument("--samples-per-study", type=int, default=8)
    parser.add_argument("--cells-per-pseudobulk", type=int, default=15)
    args = parser.parse_args()
    output = create_synthetic_pseudobulk(
        args.output,
        seed=args.seed,
        samples_per_study=args.samples_per_study,
        cells_per_pseudobulk=args.cells_per_pseudobulk,
    )
    print(output)


if __name__ == "__main__":
    main()

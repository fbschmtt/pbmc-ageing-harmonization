from __future__ import annotations

import numpy as np
import pandas as pd

from pbmc_pipeline.config import ResidualWGCNASettings
from pbmc_pipeline.residual_wgcna import run_residual_wgcna


def _settings() -> ResidualWGCNASettings:
    return ResidualWGCNASettings(
        minimum_samples=20,
        maximum_genes=100,
        minimum_module_size=10,
        soft_threshold_powers=[1, 2, 3, 4, 5, 6],
        scale_free_fit_threshold=0.8,
        mean_connectivity_threshold=10,
        module_eigengene_dissimilarity_threshold=0.2,
        dynamic_tree_deep_split=2,
        cmv_fdr_threshold=0.05,
    )


def test_residual_wgcna_publishes_modules_eigengenes_and_cmv_checks() -> None:
    rng = np.random.default_rng(7)
    n_samples = 40
    cmv_signal = rng.normal(size=n_samples)
    other_signal = rng.normal(size=n_samples)
    residuals = pd.DataFrame({
        **{
            f"cmv_gene_{index}": cmv_signal + rng.normal(scale=0.15, size=n_samples)
            for index in range(15)
        },
        **{
            f"other_gene_{index}": other_signal + rng.normal(scale=0.15, size=n_samples)
            for index in range(15)
        },
        **{f"noise_gene_{index}": rng.normal(size=n_samples) for index in range(10)},
    }, index=[f"sample-{index}" for index in range(n_samples)])
    sample_metadata = pd.DataFrame({
        "sample_key": residuals.index,
        "cmv": np.where(cmv_signal > 0, "yes", "no"),
    })
    cmv_results = pd.DataFrame({
        "gene": [f"cmv_gene_{index}" for index in range(4)],
        "padj": [0.01, 0.02, 0.03, 0.04],
    })

    result = run_residual_wgcna(
        residuals, sample_metadata, settings=_settings(), cmv_results=cmv_results
    )

    assert result.metadata["status"] == "complete"
    assert result.metadata["n_genes_selected"] == residuals.shape[1]
    assert result.metadata["selected_soft_threshold_power"] in {1, 2, 3, 4, 5, 6}
    assert {"Power", "SFT.R.sq", "mean(k)"}.issubset(result.soft_threshold)
    assert set(result.module_assignments["gene"]) == set(residuals.columns)
    assert result.module_assignments["is_cmv_associated"].sum() == 4
    assert result.eigengenes["sample_key"].tolist() == residuals.index.tolist()
    assert {
        "n_cmv_associated_genes", "cmv_gene_overlap_fdr",
        "cmv_eigengene_correlation", "cmv_eigengene_fdr",
    }.issubset(result.module_summary)
    assert result.module_summary["n_cmv_samples"].eq(n_samples).all()


def test_residual_wgcna_records_when_too_few_samples_are_available() -> None:
    residuals = pd.DataFrame(
        np.ones((3, 10)), index=["sample-1", "sample-2", "sample-3"]
    )

    result = run_residual_wgcna(residuals, pd.DataFrame(index=residuals.index), settings=_settings())

    assert result.metadata["status"] == "skipped"
    assert "at least 20" in result.metadata["reason"]

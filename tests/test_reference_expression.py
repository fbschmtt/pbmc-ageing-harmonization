import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pbmc_pipeline.reference_expression import reference_expression


def test_reference_expression_balances_studies_and_excludes_unavailable_genes():
    obs = pd.DataFrame({
        "study": ["a", "a", "b", "b", "a"],
        "sample": ["1", "2", "3", "4", "5"],
        "age": [20, 29, 25, 30, 22],
        "n_cells": [10, 12, 20, 20, 2],
        "total_counts": [100, 100, 100, 100, 100],
        "sex": ["female", "not_provided", "male", "male", "male"],
    }, index=["x1", "x2", "x3", "x4", "x5"])
    var = pd.DataFrame({"available_in_a": [True, True], "available_in_b": [True, False]}, index=["G", "H"])
    counts = sparse.csr_matrix([[10, 90], [0, 100], [40, 0], [90, 10], [80, 20]])
    adata = ad.AnnData(counts, obs=obs, var=var)
    spec = {"age_trajectory": {"reference_bin_start_age": 20, "bin_width_years": 10},
            "sample_inclusion": {"minimum_cells_per_pseudobulk_inclusive": 10, "minimum_age_years_inclusive": 20}}
    per_study, combined = reference_expression(adata, spec, cell_type="T")
    study_a = np.expm1((np.log1p(100000) + np.log1p(0)) / 2)
    expected = np.expm1((np.log1p(study_a) + np.log1p(1000000)) / 2)
    gene = combined.set_index("gene").loc["G"]
    assert gene.reference_cpm == pytest.approx(expected)
    assert gene.n_studies == 2 and gene.n_samples == 3
    assert combined.set_index("gene").loc["H", "n_studies"] == 1
    assert len(per_study) == 3
    assert set(per_study.reference_bin) == {"20-30"}


def test_reference_expression_requires_study_gene_availability():
    obs = pd.DataFrame({"study": ["a"], "sample": ["x"], "age": [25], "n_cells": [10], "total_counts": [5]}, index=["x"])
    adata = ad.AnnData(np.array([[5]]), obs=obs, var=pd.DataFrame(index=["G"]))
    spec = {"age_trajectory": {"reference_bin_start_age": 20, "bin_width_years": 10}, "sample_inclusion": {"minimum_cells_per_pseudobulk_inclusive": 10, "minimum_age_years_inclusive": 20}}
    with pytest.raises(ValueError, match="availability"):
        reference_expression(adata, spec, cell_type="T")

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline.gene_profile_report import generate_gene_profile_report
from scripts.profile_gene import _fit_profile

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = ROOT / "config" / "pipeline.json"
def test_profile_gene_rejects_study_specific_gene_before_fitting(tmp_path) -> None:
    rows = []
    for study in ("synthetic_study_a", "synthetic_study_b"):
        for index in range(4):
            rows.append({
                "study": study,
                "study_site": f"{study}_site",
                "sample": f"{study}_sample_{index}",
                "age": 40 + index * 10,
                "sex": "female" if index % 2 == 0 else "male",
                "n_cells": 15,
                "total_counts": 100 + index,
                "aifi_l2_majority": "CD14 monocyte",
            })
    obs = pd.DataFrame(rows, index=[f"pb_{i}" for i in range(len(rows))])
    var = pd.DataFrame(
        {
            "available_in_synthetic_study_a": [True],
            "available_in_synthetic_study_b": [False],
        },
        index=["SYNTH_STUDY_A_ONLY"],
    )
    input_path = tmp_path / "pseudobulk.h5ad"
    ad.AnnData(X=np.full((len(obs), 1), 20, dtype=np.int64), obs=obs, var=var).write_h5ad(
        input_path,
    )
    with pytest.raises(ValueError, match="not available in every study"):
        _fit_profile(
            input_path=input_path,
            config_path=PIPELINE_PATH,
            cell_type="CD14 monocyte",
            gene="SYNTH_STUDY_A_ONLY",
            output_dir=tmp_path / "profiles",
            cpus=1,
        )


def test_profile_report_renderer_requires_fit_summary(tmp_path) -> None:
    with pytest.raises(FileNotFoundError, match="fit_summary.json"):
        generate_gene_profile_report(
            profile_dir=tmp_path / "missing",
            template_path=ROOT / "reports" / "gene_profile_report.py",
            project_root=ROOT,
        )

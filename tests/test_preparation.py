from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline.preparation import PreparationError, prepare_study, read_prepared_cells

ROOT = Path(__file__).resolve().parents[1]


def _schema():
    import json

    return json.loads((ROOT / "config/harmonized_obs_schema.json").read_text())


def test_onek_adapter_materializes_exactly_one_canonical_row_per_cell(tmp_path):
    source = ad.AnnData(X=np.ones((2, 1)))
    source.obs = pd.DataFrame(
        {
            "cell_type": ["T cell", "B cell"], "donor_id": ["d1", "d2"],
            "age": [40, 50], "sex": ["Female", "Male"], "pool_number": [1, 2],
        }, index=["cell_a", "cell_b"],
    )
    expression = tmp_path / "source.h5ad"
    source.write_h5ad(expression)
    output = tmp_path / "onek.cells.csv.gz"
    report = tmp_path / "onek.prepare.json"

    summary = prepare_study(
        ROOT, "onek1k", {"preparation": {"adapter": "onek1k"}}, _schema(),
        expression, output, report,
    )
    prepared = read_prepared_cells(output, source.obs_names, _schema())

    assert summary["n_cells"] == 2
    assert prepared.index.tolist() == ["cell_a", "cell_b"]
    assert prepared["sample"].tolist() == ["d1", "d2"]
    assert prepared["batch_single_cell"].tolist() == ["onek1k_1", "onek1k_2"]
    assert prepared["study"].tolist() == ["onek1k", "onek1k"]


def test_prepared_metadata_rejects_missing_or_extra_cell_ids(tmp_path):
    table = pd.DataFrame(
        {"sample": ["s1"], "subject": ["d1"]}, index=pd.Index(["cell_a"], name="cell_id")
    )
    path = tmp_path / "prepared.csv.gz"
    table.to_csv(path, compression="infer")

    with pytest.raises(PreparationError, match="missing=1"):
        read_prepared_cells(path, pd.Index(["cell_a", "cell_b"]), _schema())

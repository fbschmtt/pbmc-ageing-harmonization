from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline.preparation import PreparationError, prepare_study, read_prepared_cells
from pbmc_pipeline.studies import nehar_belaid26

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


def test_nehar_belaid_adapter_excludes_reused_public_cohorts_and_uses_supplement(monkeypatch):
    source = pd.DataFrame(
        {
            "Study": ["Nehar-Belaid_et_al", "Deng_et_al_2025", "Nehar-Belaid_et_al"],
            "sample_id": ["HI1", "public_1", "HO1"],
            "LS_L3": ["CD4_naive", "CD8_memory", "CD14_mono"],
        },
        index=["own_1", "public_1", "own_2"],
    )
    supplement = pd.DataFrame(
        {
            "Names": ["HI1", "HO1"], "IDs": ["JB1", "JB2"],
            "Age (in months)": [6, 960], "BMI": ["UD", 25.0],
            "CMV status": ["UD", "Positive"], "Ethnicity": ["W", "AA"],
            "Sex": ["F", "M"], "Site of blood collection": ["NCH", "HSNRI"],
            "runs_10x": ["B1", "B2"], "Chemistry_10x": ["V2", "V3"],
        }
    )
    monkeypatch.setattr(nehar_belaid26.pd, "read_excel", lambda *args, **kwargs: supplement)

    prepared = nehar_belaid26.prepare_cells(source, ROOT)

    assert prepared.index.tolist() == ["own_1", "own_2"]
    assert prepared["sample"].tolist() == ["HI1", "HO1"]
    assert prepared["subject"].tolist() == ["JB1", "JB2"]
    assert prepared["age"].tolist() == [0.5, 80.0]
    assert np.isnan(prepared.loc["own_1", "bmi"])
    assert prepared["cmv"].tolist() == ["not_provided", "yes"]
    assert prepared["country"].tolist() == ["usa", "canada"]
    assert prepared["technology"].tolist() == ["10X3'v2", "10X3'v3"]


def test_prepared_metadata_allows_declared_cell_subset(tmp_path):
    schema = _schema()
    table = pd.DataFrame(
        {
            column: [0.0 if expected_type == "float" else "value"]
            for column, expected_type in schema["required"].items()
            if column not in {"aifi_l1_majority", "aifi_l2_majority", "aifi_l3_majority"}
        },
        index=pd.Index(["cell_a"], name="cell_id"),
    )
    table["sample"] = "s1"
    table["subject"] = "d1"
    path = tmp_path / "prepared.csv.gz"
    table.to_csv(path, compression="infer")

    prepared = read_prepared_cells(
        path, pd.Index(["cell_a", "cell_b"]), schema, allow_cell_subset=True
    )

    assert prepared.index.tolist() == ["cell_a"]

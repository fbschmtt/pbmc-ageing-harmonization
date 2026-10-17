from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline.preparation import PreparationError, prepare_study, read_prepared_cells
from pbmc_pipeline.studies import aida25, cima, fachrul26, nehar_belaid26, perez22, terekhova26

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
            "assay": ["10x 3' v2", "10x 3' v2"],
        }, index=["cell_a", "cell_b"],
    )
    expression = tmp_path / "source.h5ad"
    source.write_h5ad(expression)
    output = tmp_path / "onek.cells.csv.gz"
    source_obs_output = tmp_path / "onek.source_obs.csv.gz"
    report = tmp_path / "onek.prepare.json"

    summary = prepare_study(
        ROOT, "onek1k", {
            "publication_doi": "10.1126/science.abf3041",
            "preparation": {"adapter": "onek1k"},
        }, _schema(),
        expression, output, report, source_obs_output_path=source_obs_output,
    )
    prepared = read_prepared_cells(output, source.obs_names, _schema())

    assert summary["n_cells"] == 2
    assert prepared.index.tolist() == ["cell_a", "cell_b"]
    assert prepared["sample"].tolist() == ["d1", "d2"]
    assert prepared["batch_single_cell"].tolist() == ["onek1k_1", "onek1k_2"]
    assert prepared["technology"].tolist() == ["10X3'v2", "10X3'v2"]
    assert prepared["study"].tolist() == ["onek1k", "onek1k"]
    assert prepared["study_doi"].tolist() == ["10.1126/science.abf3041"] * 2
    source_obs = pd.read_csv(source_obs_output, index_col="cell_id")
    assert source_obs.index.tolist() == ["cell_a", "cell_b"]
    assert source_obs["cell_type"].tolist() == ["T cell", "B cell"]
    assert summary["source_obs_output"] == str(source_obs_output)


def test_fachrul_adapter_uses_embedded_metadata():
    source = pd.DataFrame(
        {
            "cell_type": ["CD14-positive monocyte"], "sample_id": ["sample_1"],
            "donor_id": ["donor_1"], "development_stage": ["54-year-old stage"],
            "sex": ["male"], "self_reported_ethnicity": ["Indonesian"],
            "assay": ["10x 5' v2"], "alignment_software": ["Cell Ranger count v7.2.0"],
            "reference_genome": ["GRCh38"], "sample_preservation_method": ["fresh"],
            "intronic_reads_counted": ["yes"], "Village": ["Pedawa"],
            "library_id": ["library_1"],
        },
        index=["cell_a"],
    )

    prepared = fachrul26.prepare_cells(source, ROOT)

    assert prepared.loc["cell_a", "sample"] == "sample_1"
    assert prepared.loc["cell_a", "age"] == 54.0
    assert prepared.loc["cell_a", "technology"] == "10X5'v2"
    assert prepared.loc["cell_a", "frozen"] == "yes"
    assert prepared.loc["cell_a", "study_site"] == "fachrul26_pedawa"
    assert prepared.loc["cell_a", "batch_single_cell"] == "fachrul26_library_1"


def test_cima_adapter_uses_embedded_metadata_and_highest_resolution_label():
    source = pd.DataFrame(
        {
            "cell_type_l4": ["CD8_Tem_CCR7neg"], "sample": ["CIMA_H001"],
            "age": [42], "sex": ["Female"], "BMI": [22.5],
            "smoking_ststus": ["past"], "library": ["CIMA_RNA_H001_1"],
        }, index=["cell_a"],
    )

    prepared = cima.prepare_cells(source, ROOT)

    assert prepared.loc["cell_a", "study_celltype"] == "CD8_Tem_CCR7neg"
    assert prepared.loc["cell_a", "sample"] == "CIMA_H001"
    assert prepared.loc["cell_a", "technology"] == "DNBelab C Series"
    assert prepared.loc["cell_a", "smoking_status"] == "no"
    assert prepared.loc["cell_a", "batch_single_cell"] == "cima_CIMA_RNA_H001_1"
    assert prepared.loc["cell_a", "aligner"] == "STAR 2.7.2b"
    assert prepared.loc["cell_a", "genome"] == "GRCh38"
    assert prepared.loc["cell_a", "frozen"] == "not_provided"
    assert prepared.loc["cell_a", "include_intronic"] == "no"


def test_terekhova26_adapter_uses_embedded_author_labels():
    source = pd.DataFrame(
        {
            "author_celltype": ["CD8+ Tem GZMB+"], "donor_id": ["donor_1"],
            "age": [66], "sex": ["Male"], "BMI": [27.7], "CMV": ["Positive"],
            "ethnicity": ["European"], "Pool": ["Pool_16"],
        }, index=["cell_a"],
    )

    prepared = terekhova26.prepare_cells(source, ROOT)

    assert prepared.loc["cell_a", "study_celltype"] == "CD8+ Tem GZMB+"
    assert prepared.loc["cell_a", "sample"] == "donor_1"
    assert prepared.loc["cell_a", "cmv"] == "yes"
    assert prepared.loc["cell_a", "batch_single_cell"] == "terekhova26_Pool_16"
    assert prepared.loc["cell_a", "aligner"] == "Cell Ranger multi 9.0.1"
    assert prepared.loc["cell_a", "genome"] == "GRCh38-2020-A"
    assert prepared.loc["cell_a", "demultiplexing"] == "genetic+HTO"
    assert prepared.loc["cell_a", "frozen"] == "yes"
    assert prepared.loc["cell_a", "include_intronic"] == "yes"


def test_aida_adapter_uses_level4_labels_and_prefixes_site(monkeypatch):
    source = pd.DataFrame(
        {
            "donor_id": ["donor_1"], "Country": ["IN"],
            "Annotation_Level4": ["CD4_T_CM_IL7Rhi"], "development_stage": ["45-year-old stage"],
            "sex": ["Female"], "BMI": [24.0], "self_reported_ethnicity": ["Asian"],
            "assay": ["10x 5' v2"], "alignment_software": ["Cell Ranger"],
            "reference_genome": ["GRCh38"], "intronic_reads_counted": ["yes"],
        },
        index=["cell_a"],
    )
    metadata = pd.DataFrame(
        {
            "DCP_ID": ["donor_1"], "Country": ["India"],
            "scRNA-seq Experimental Batch": ["batch_1"], "Smoking Status": ["0"],
        }
    )
    monkeypatch.setattr(aida25.pd, "read_excel", lambda *args, **kwargs: metadata)

    prepared = aida25.prepare_cells(source, ROOT)

    assert prepared.loc["cell_a", "study_celltype"] == "CD4_T_CM_IL7Rhi"
    assert prepared.loc["cell_a", "study_site"] == "aida_IN"


def test_perez_adapter_uses_embedded_metadata():
    source = pd.DataFrame(
        {
            "cell_type": ["classical monocyte"], "sample_uuid": ["sample_1"],
            "donor_id": ["donor_1"], "development_stage": ["33-year-old stage"],
            "sex": ["female"], "self_reported_ethnicity": ["European American"],
            "assay": ["10x 3' v2"], "library_uuid": ["library_1"], "disease": ["normal"],
        },
        index=["cell_a"],
    )

    prepared = perez22.prepare_cells(source, ROOT)

    assert prepared.loc["cell_a", "sample"] == "sample_1"
    assert prepared.loc["cell_a", "age"] == 33.0
    assert prepared.loc["cell_a", "technology"] == "10X3'v2"
    assert prepared.loc["cell_a", "genome"] == "GRCh37"
    assert prepared.loc["cell_a", "batch_single_cell"] == "perez22_library_1"
    assert prepared.loc["cell_a", "disease_status"] == "healthy"


def test_perez_adapter_excludes_sle_cells_before_canonical_metadata():
    source = pd.DataFrame(
        {
            "disease": ["normal", "systemic lupus erythematosus"],
            "cell_type": ["B cell", "T cell"], "sample_uuid": ["healthy", "sle"],
            "donor_id": ["healthy", "sle"], "development_stage": ["30-year-old stage"] * 2,
            "sex": ["female"] * 2, "self_reported_ethnicity": ["European"] * 2,
            "assay": ["10x 3' v2"] * 2, "library_uuid": ["one", "two"],
        },
        index=["healthy_cell", "sle_cell"],
    )

    prepared = perez22.prepare_cells(source, ROOT)

    assert prepared.index.tolist() == ["healthy_cell"]
    assert prepared["disease_status"].tolist() == ["healthy"]


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
            "LS_L4": ["CD4_naive_CCR7+", "CD8_memory_GZMK+", "CD14_mono_S100A8+"],
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
    assert prepared["study_site"].tolist() == ["nehar_belaid26_nch", "nehar_belaid26_hsnri"]
    assert prepared["technology"].tolist() == ["10X3'v2", "10X3'v3"]
    assert prepared["study_celltype"].tolist() == ["CD4_naive_CCR7+", "CD14_mono_S100A8+"]


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

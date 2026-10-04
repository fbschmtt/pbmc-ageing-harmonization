"""Prepare the in-house Nehar-Belaid et al. (2026) cohort only."""
from __future__ import annotations

import anndata as ad
import pandas as pd

from ..metadata import MetadataError
from .common import canonical, safe_left_join

STUDY_COLUMN = "Study"
STUDY_VALUE = "Nehar-Belaid_et_al"
METADATA_PATH = "input_data/nehar_belaid26/41467_2026_73729_MOESM3_ESM.xls"
LABEL_PATH = "input_data/nehar_belaid26/GSE233321_all_PBMCs.h5ad"
LABEL_COLUMNS = ("LS_L1", "LS_L2", "LS_L3", "LS_L4")


def prepare_cells(source, root):
    """Join Supplementary Data 1a to the in-house source cells by sample ID."""
    if STUDY_COLUMN in source:
        cells = source.loc[source[STUDY_COLUMN].astype("string").eq(STUDY_VALUE)].copy()
        if cells.empty:
            raise MetadataError(f"Nehar-Belaid source contains no {STUDY_VALUE!r} cells")
    else:
        cells = source.copy()

    table = pd.read_excel(root / METADATA_PATH, sheet_name="1a", header=1)
    table = table.dropna(subset=["Names"]).copy()
    table["sample_id"] = table["Names"].astype("string")
    table["subject_id"] = table["IDs"].astype("string")
    table["age_years"] = pd.to_numeric(table["Age (in months)"], errors="raise") / 12
    table["bmi_numeric"] = pd.to_numeric(table["BMI"], errors="coerce")
    table["cmv_normalized"] = table["CMV status"].astype("string").map(
        {"Positive": "yes", "Negative": "no", "UD": "not_provided"}
    ).fillna("not_provided")
    table["study_site_normalized"] = table["Site of blood collection"].astype("string").str.lower()
    table["country_normalized"] = table["Site of blood collection"].astype("string").map(
        {"NCH": "usa", "UCHC": "usa", "HSNRI": "canada"}
    ).fillna("not_provided")
    table["technology_normalized"] = table["Chemistry_10x"].astype("string").str.upper().map(
        {"V2": "10X3'v2", "V3": "10X3'v3"}
    )
    if table["technology_normalized"].isna().any():
        unexpected = sorted(table.loc[table["technology_normalized"].isna(), "Chemistry_10x"].dropna().unique())
        raise MetadataError(f"Nehar-Belaid Supplementary Data 1a has unknown 10x chemistry: {unexpected}")

    if "subject_id" in cells:
        joined = safe_left_join(cells, table, ["subject_id"])
        if joined["Names"].isna().any():
            missing = joined.loc[joined["Names"].isna(), "subject_id"].nunique()
            raise MetadataError(f"Nehar-Belaid Supplementary Data 1a misses {missing} archive subjects")
        joined = _attach_h5ad_labels(joined, root)
    elif "sample_id" in cells:
        joined = safe_left_join(cells, table, ["sample_id"])
        if joined["IDs"].isna().any():
            missing = joined.loc[joined["IDs"].isna(), "sample_id"].nunique()
            raise MetadataError(f"Nehar-Belaid Supplementary Data 1a misses {missing} retained samples")
    else:
        raise MetadataError("Nehar-Belaid source is missing 'subject_id' or 'sample_id'")

    return canonical(joined, "nehar_belaid26", {
        "study_celltype": {"source": "LS_L3"},
        "sample": {"source": "sample_id"},
        "subject": {"source": "IDs"},
        "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "age_years", "dtype": "float"},
        "sex": {"source": "Sex", "map": {"F": "female", "M": "male"}, "lower": True},
        "bmi": {"source": "bmi_numeric", "dtype": "float"},
        "cmv": {"source": "cmv_normalized"},
        "ethnicity": {"source": "Ethnicity", "lower": True},
        "ethnicity_fine": {"source": "Ethnicity", "lower": True},
        "technology": {"source": "technology_normalized"},
        "aligner": {"constant": "Cell Ranger 3.0.2"},
        "genome": {"constant": "hg19"},
        # inferred: no demultiplexing method is reported for this release.
        "demultiplexing": {"constant": "none"},
        # inferred: thawed PBMCs imply frozen input material.
        "frozen": {"constant": "yes"},
        # inferred: intronic reads are not reported for this release.
        "include_intronic": {"constant": "no"},
        "smoking_status": {"constant": "not_provided"},
        "disease_status": {"constant": "healthy"},
        "study_site": {"source": "study_site_normalized", "prefix": "nehar_belaid26_"},
        "country": {"source": "country_normalized"},
        "batch_single_cell": {"source": "runs_10x", "prefix": "nehar_belaid26_"},
    })


def _attach_h5ad_labels(cells: pd.DataFrame, root) -> pd.DataFrame:
    """Transfer published labels only for exact sample-plus-barcode matches."""
    if "barcode" not in cells:
        return cells
    labels = ad.read_h5ad(root / LABEL_PATH, backed="r")
    try:
        available = [column for column in LABEL_COLUMNS if column in labels.obs]
        table = labels.obs.loc[:, ["sample_id", *available]].copy()
    finally:
        labels.file.close()
    # H5AD barcode suffixes were reassigned while its source cohorts were
    # concatenated (for example, raw ``-1`` becomes H5AD ``-68``). The 16-base
    # sequence remains stable, and the sample ID keeps the transfer identity-safe.
    table["barcode_sequence"] = table.index.astype(str).str.extract(r"([ACGT]+)-\d+", expand=False)
    table["_key"] = table["sample_id"].astype(str) + "\0" + table["barcode_sequence"]
    table = table.dropna(subset=["barcode_sequence"]).drop_duplicates("_key", keep=False).set_index("_key")
    barcode_sequence = cells["barcode"].astype(str).str.extract(r"([ACGT]+)-\d+", expand=False)
    keys = cells["sample_id"].astype(str) + "\0" + barcode_sequence
    result = cells.copy()
    for column in LABEL_COLUMNS:
        values = table[column] if column in table else pd.Series(dtype="string")
        result[column] = keys.map(values).astype("string").fillna("not_provided").to_numpy()
    return result

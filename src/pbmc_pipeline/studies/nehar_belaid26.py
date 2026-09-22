"""Prepare the in-house Nehar-Belaid et al. (2026) cohort only."""
from __future__ import annotations

import pandas as pd

from ..metadata import MetadataError
from .common import canonical, safe_left_join

STUDY_COLUMN = "Study"
STUDY_VALUE = "Nehar-Belaid_et_al"
METADATA_PATH = "input_data/nehar_belaid26/41467_2026_73729_MOESM3_ESM.xls"


def prepare_cells(source, root):
    """Exclude reused public cohorts and join Supplementary Data 1a by sample ID."""
    if STUDY_COLUMN not in source:
        raise MetadataError(f"Nehar-Belaid source is missing required {STUDY_COLUMN!r} column")
    cells = source.loc[source[STUDY_COLUMN].astype("string").eq(STUDY_VALUE)].copy()
    if cells.empty:
        raise MetadataError(f"Nehar-Belaid source contains no {STUDY_VALUE!r} cells")

    table = pd.read_excel(root / METADATA_PATH, sheet_name="1a", header=1)
    table = table.dropna(subset=["Names"]).copy()
    table["sample_id"] = table["Names"].astype("string")
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

    if "sample_id" not in cells:
        raise MetadataError("Nehar-Belaid source is missing required 'sample_id' column")
    joined = safe_left_join(cells, table, ["sample_id"])
    if joined["IDs"].isna().any():
        missing = joined.loc[joined["IDs"].isna(), "sample_id"].nunique()
        raise MetadataError(f"Nehar-Belaid Supplementary Data 1a misses {missing} retained samples")

    # Supplementary Data 1a records V2 and V3 in ``Chemistry_10x``. Preserve
    # that most-specific reported technology as 10X3'v2/10X3'v3; ``runs_10x``
    # remains available as the per-cell sequencing batch.
    #
    # Cell Ranger 3.0.2 and hg19 are reported in the paper's processing
    # methods.  The remaining technical values below are documented in config
    # provenance as inferences: thawed PBMCs imply frozen input, while neither
    # intronic reads nor sample demultiplexing are reported for this release.
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
        "demultiplexing": {"constant": "none"},
        "frozen": {"constant": "yes"},
        "include_intronic": {"constant": "no"},
        "smoking_status": {"constant": "not_provided"},
        "study_site": {"source": "study_site_normalized"},
        "country": {"source": "country_normalized"},
        "batch_single_cell": {"source": "runs_10x", "prefix": "nehar_belaid26_"},
    })

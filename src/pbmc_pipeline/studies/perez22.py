from __future__ import annotations

from ..metadata import MetadataError
from .common import canonical

DISEASE_STATUS_COLUMN = "_disease_status"
DISEASE_STATUS_MAP = {
    "normal": "healthy",
    "systemic lupus erythematosus": "systemic_lupus_erythematosus",
}


def prepare_cells(source, root):
    """Retain embedded-metadata healthy Perez22 cells and map their disease status."""
    del root
    if "disease" not in source:
        raise MetadataError("Perez22 source is missing required 'disease' column")
    cells = source.copy()
    cells[DISEASE_STATUS_COLUMN] = cells["disease"].astype("string").str.lower().map(DISEASE_STATUS_MAP)
    unknown = sorted(cells.loc[cells[DISEASE_STATUS_COLUMN].isna(), "disease"].dropna().astype(str).unique())
    if unknown:
        raise MetadataError(f"Perez22 source has unrecognized disease values: {unknown}")
    # Selection is performed before canonical metadata is written, just as for
    # reused cohorts in Nehar-Belaid: excluded SLE cells must not leak into any
    # downstream artifact through a metadata-only filter.
    cells = cells.loc[cells[DISEASE_STATUS_COLUMN].eq("healthy")].copy()
    if cells.empty:
        raise MetadataError("Perez22 source contains no healthy cells")
    return canonical(cells, "perez22", {
        "study_celltype": {"source": "cell_type"},
        "sample": {"source": "sample_uuid"},
        "subject": {"source": "donor_id"},
        "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "development_stage", "replace": {"-year-old stage": ""}, "dtype": "float"},
        "sex": {"source": "sex", "lower": True},
        "bmi": {"constant": None},
        "cmv": {"constant": "not_provided"},
        "ethnicity": {"source": "self_reported_ethnicity", "lower": True},
        "ethnicity_fine": {"source": "self_reported_ethnicity", "lower": True},
        "technology": {"source": "assay", "map": {"10x 3' v2": "10X3'v2"}},
        "aligner": {"constant": "not_provided"},
        # inferred: GENCODE 19 implies the GRCh37-era reference assembly.
        "genome": {"constant": "GRCh37"},
        "demultiplexing": {"constant": "not_provided"},
        "frozen": {"constant": "not_provided"},
        "include_intronic": {"constant": "not_provided"},
        "smoking_status": {"constant": "not_provided"},
        "disease_status": {"source": DISEASE_STATUS_COLUMN},
        "study_site": {"constant": "perez22"},
        "country": {"constant": "not_provided"},
        "batch_single_cell": {"source": "library_uuid", "prefix": "perez22_"},
    })

from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    del root
    return canonical(source, "onek1k", {
        "study_celltype": {"source": "cell_type"},
        "sample": {"source": "donor_id"},
        "subject": {"source": "donor_id"},
        "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "age", "dtype": "float"},
        "sex": {"source": "sex", "lower": True},
        "bmi": {"constant": None},
        "cmv": {"constant": "not_provided"},
        "ethnicity": {"constant": "caucasian"},
        "ethnicity_fine": {"constant": "northern_european"},
        "technology": {"source": "assay", "map": {"10x 3' v2": "10X3'v2"}},
        "aligner": {"constant": "STAR"},
        # inferred: retained from the study's legacy processing configuration.
        "genome": {"constant": "GRCh37"},
        "demultiplexing": {"constant": "genetic"},
        # inferred: source protocol interpretation; confirm frozen-state metadata.
        "frozen": {"constant": "no"},
        # inferred: source processing metadata does not expose this per cell.
        "include_intronic": {"constant": "no"},
        "smoking_status": {"constant": "not_provided"},
        "disease_status": {"constant": "healthy"},
        "study_site": {"constant": "onek1k"},
        "country": {"constant": "aus"},
        "batch_single_cell": {"source": "pool_number", "prefix": "onek1k_"},
    })

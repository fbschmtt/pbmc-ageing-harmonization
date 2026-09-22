from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    """Use the embedded CellxGene observation metadata for Perez22."""
    del root
    return canonical(source, "perez22", {
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
        "study_site": {"constant": "perez22"},
        "country": {"constant": "not_provided"},
        "batch_single_cell": {"source": "library_uuid", "prefix": "perez22_"},
    })

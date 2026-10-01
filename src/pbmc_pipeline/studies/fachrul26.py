from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    """Use the embedded CellxGene observation metadata for Fachrul26."""
    del root
    return canonical(source, "fachrul26", {
        "study_celltype": {"source": "cell_type"},
        "sample": {"source": "sample_id"},
        "subject": {"source": "donor_id"},
        "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "development_stage", "replace": {"-year-old stage": ""}, "dtype": "float"},
        "sex": {"source": "sex", "lower": True},
        "bmi": {"constant": None},
        "cmv": {"constant": "not_provided"},
        "ethnicity": {"source": "self_reported_ethnicity", "lower": True},
        "ethnicity_fine": {"source": "self_reported_ethnicity", "lower": True},
        "technology": {"source": "assay", "map": {"10x 5' v2": "10X5'v2"}},
        "aligner": {"source": "alignment_software"},
        "genome": {"source": "reference_genome"},
        "demultiplexing": {"constant": "not_provided"},
        # Paper methods: cryopreserved PBMCs were thawed. CELLxGENE obs conflicts
        # and says sample_preservation_method="fresh"; preserve both facts in provenance.
        "frozen": {"constant": "yes"},
        "include_intronic": {"source": "intronic_reads_counted", "lower": True},
        "smoking_status": {"constant": "not_provided"},
        "disease_status": {"constant": "healthy"},
        "study_site": {"source": "Village", "lower": True},
        # inferred: source description identifies this as an Indonesian cohort.
        "country": {"constant": "indonesia"},
        "batch_single_cell": {"source": "library_id", "prefix": "fachrul26_"},
    })

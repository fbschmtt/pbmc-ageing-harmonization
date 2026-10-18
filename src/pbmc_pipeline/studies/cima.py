from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    """Use the complete cell and donor metadata embedded in the CIMA H5AD."""
    del root
    return canonical(source, "cima", {
        "study_celltype": {"source": "cell_type_l4"},
        "sample": {"source": "sample"},
        "subject": {"source": "sample"},
        "age": {"source": "age", "dtype": "float"},
        "sex": {"source": "sex", "lower": True},
        "bmi": {"source": "BMI", "dtype": "float"},
        "cmv": {"constant": "not_provided"},
        # inferred: CIMA is a Chinese cohort; no per-donor ancestry field is in obs.
        "ethnicity": {"constant": "chinese"},
        "ethnicity_fine": {"constant": "chinese"},
        "technology": {"constant": "DNBelab C Series"},
        "aligner": {"constant": "STAR 2.7.2b"},
        "genome": {"constant": "GRCh38"},
        # inferred: no demultiplexing method is described in the paper.
        "demultiplexing": {"constant": "none"},
        "frozen": {"constant": "not_provided"},
        "include_intronic": {"constant": "no"},
        "smoking_status": {
            "source": "smoking_ststus",
            "map": {"current": "yes", "never": "no", "past": "no"},
            "default": "not_provided",
        },
        "disease_status": {"constant": "healthy"},
        "study_site": {"constant": "cima"},
        "country": {"constant": "china"},
        "batch_single_cell": {"source": "RNA_experiment", "prefix": "cima_"},
    })

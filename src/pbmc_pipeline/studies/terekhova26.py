from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    """Use the author-reannotated KCL observation metadata without a sidecar join."""
    del root
    return canonical(source, "terekhova26", {
        "study_celltype": {"source": "author_celltype"},
        "sample": {"source": "donor_id"},
        "subject": {"source": "donor_id"},
        "age": {"source": "age", "dtype": "float"},
        "sex": {"source": "sex", "lower": True},
        "bmi": {"source": "BMI", "dtype": "float"},
        "cmv": {
            "source": "CMV",
            "map": {"Positive": "yes", "Negative": "no", "No_info": "not_provided"},
            "default": "not_provided",
        },
        "ethnicity": {"source": "ethnicity", "lower": True},
        "ethnicity_fine": {"source": "ethnicity", "lower": True},
        # inferred: the EGA catalogue reports 10x Genomics scRNA-seq but not its chemistry.
        "technology": {"constant": "10X"},
        "aligner": {"constant": "Cell Ranger multi 9.0.1"},
        "genome": {"constant": "GRCh38-2020-A"},
        # Paper methods: genetic demultiplexing, with HTO antibodies mapping genotypes to samples.
        "demultiplexing": {"constant": "genetic+HTO"},
        "frozen": {"constant": "yes"},
        # inferred: Cell Ranger v9 uses the v6+ default arguments that include intronic reads.
        "include_intronic": {"constant": "yes"},
        "smoking_status": {"constant": "not_provided"},
        "disease_status": {"constant": "healthy"},
        "study_site": {"constant": "terekhova26_kcl"},
        "country": {"constant": "united_kingdom"},
        "batch_single_cell": {"source": "Pool", "prefix": "terekhova26_"},
    })

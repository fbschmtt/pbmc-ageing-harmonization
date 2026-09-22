from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    del root
    return canonical(source, "aifi", {
        "study_celltype": {"source": "AIFI_L3"},
        "sample": {"source": "specimen.specimenGuid"},
        "subject": {"source": "specimen.specimenGuid", "split": {"separator": "-", "index": 0}},
        "sampling_timepoint": {"source": "sample.visitName"},
        "age": {"source": "sample.subjectAgeAtDraw", "dtype": "float"},
        "sex": {"source": "subject.biologicalSex", "lower": True},
        "bmi": {"source": "subject.bmi", "dtype": "float"},
        "cmv": {"source": "subject.cmv", "lower": True},
        "ethnicity": {"source": "subject.race", "lower": True},
        "ethnicity_fine": {"concat": ["subject.race", "subject.ethnicity"], "separator": "_", "lower": True},
        "technology": {"constant": "10X3'"}, "aligner": {"constant": "not_provided"},
        "genome": {"constant": "not_provided"}, "demultiplexing": {"constant": "HTO"},
        "frozen": {"constant": "yes"}, "include_intronic": {"constant": "not_provided"},
        "smoking_status": {"constant": "not_provided"}, "study_site": {"constant": "aifi"},
        "country": {"constant": "usa"}, "batch_single_cell": {"source": "well_id"},
        "aifi_l1_majority": {"source": "AIFI_L1"}, "aifi_l2_majority": {"source": "AIFI_L2"},
        "aifi_l3_majority": {"source": "AIFI_L3"},
    })

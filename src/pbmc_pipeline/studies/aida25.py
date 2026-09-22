from __future__ import annotations

import pandas as pd

from .common import canonical, safe_left_join


def prepare_cells(source, root):
    table = pd.read_excel(root / "input_data/aida25/mmc1.xlsx", skiprows=1)
    table = table.rename(columns={"DCP_ID": "donor_id", "scRNA-seq Experimental Batch": "batch_single_cell"})
    table["Country"] = table["Country"].replace({"India": "IN", "Japan": "JP", "Singapore": "SG", "South Korea": "KR", "Thailand": "TH"})
    duplicate = table.duplicated(["Country", "donor_id"], keep=False)
    table.loc[duplicate, "batch_single_cell"] = "unclear"
    table = table.drop_duplicates(["Country", "donor_id"], keep="first")
    joined = safe_left_join(source, table, ["donor_id", "Country"])
    return canonical(joined, "aida25", {
        "study_celltype": {"source": "Annotation_Level2"},
        "sample": {"concat": [{"constant": "aida"}, "Country", "donor_id"], "separator": "_"},
        "subject": {"source": "donor_id"}, "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "development_stage", "replace": {"-year-old stage": ""}, "dtype": "float"},
        "sex": {"source": "sex", "lower": True}, "bmi": {"source": "BMI", "dtype": "float"},
        "cmv": {"constant": "not_provided"},
        "ethnicity": {"source": "self_reported_ethnicity", "map": {"European": "caucasian"}, "default": "asian"},
        "ethnicity_fine": {"source": "self_reported_ethnicity", "lower": True},
        "technology": {"constant": "10X5'"}, "aligner": {"source": "alignment_software"},
        "genome": {"source": "reference_genome"}, "demultiplexing": {"constant": "genetic"},
        "frozen": {"constant": "yes"}, "include_intronic": {"source": "intronic_reads_counted"},
        "smoking_status": {"source": "Smoking Status", "map": {"0": "no", "1": "yes"}, "default": "not_provided"},
        "study_site": {"concat": [{"constant": "aida"}, "Country"], "separator": "_"},
        "country": {"source": "Country", "lower": True},
        "batch_single_cell": {"source": "batch_single_cell", "prefix": "aida_"},
    })

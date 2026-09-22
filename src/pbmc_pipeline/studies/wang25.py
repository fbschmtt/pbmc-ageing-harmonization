from __future__ import annotations

import pandas as pd

from .common import canonical, safe_left_join


def prepare_cells(source, root):
    cells = source.copy()
    cells["Sample ID"] = cells["sampleName"].astype("string").str.replace("sample", "sc_sample_", regex=False)
    table = pd.read_excel(
        root / "input_data/metadata/wang25/41590_2024_2059_MOESM3_ESM.xlsx",
        sheet_name=1, skiprows=2,
    )
    joined = safe_left_join(cells, table, ["Sample ID"])
    return canonical(joined, "wang25", {
        "study_celltype": {"source": "secondary_type"}, "sample": {"source": "Sample ID"},
        "subject": {"source": "Sample ID"}, "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "Age (y)", "dtype": "float"}, "sex": {"source": "Gender", "lower": True},
        "bmi": {"source": "BMI (kg/m2)", "dtype": "float"},
        "cmv": {"source": "CMV IgM antibodies", "lower": True}, "ethnicity": {"constant": "asian"},
        "ethnicity_fine": {"constant": "mainland_china"}, "technology": {"constant": "10X5'"},
        "aligner": {"constant": "CRv3.0.2"}, "genome": {"constant": "GRCh38"},
        "demultiplexing": {"constant": "none"}, "frozen": {"constant": "no"},
        "include_intronic": {"constant": "no"}, "smoking_status": {"constant": "no"},
        "study_site": {"constant": "wang25"}, "country": {"constant": "cn"},
        "batch_single_cell": {"constant": "not_provided"},
    })

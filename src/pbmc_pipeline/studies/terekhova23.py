from __future__ import annotations

import pandas as pd

from ..metadata import MetadataError
from .common import canonical, safe_left_join


def prepare_cells(source, root):
    lookup = pd.read_csv(root / "input_data/metadata/terekhova23/cell_to_tube.csv.gz", usecols=["cell_id", "Tube_id"])
    if lookup["cell_id"].duplicated().any():
        raise MetadataError("Terekhova cell-to-tube lookup contains duplicate cell IDs")
    lookup = lookup.set_index("cell_id")
    missing = source.index.difference(lookup.index)
    if len(missing):
        raise MetadataError(f"Terekhova cell-to-tube lookup misses {len(missing)} source cells")
    cells = source.copy()
    cells["Tube_id"] = lookup.loc[cells.index, "Tube_id"].to_numpy()

    table = pd.read_excel(root / "input_data/metadata/terekhova23/mmc2.xlsx")
    visits = _reshape_visits(table)
    joined = safe_left_join(cells, visits, ["Tube_id"])
    return canonical(joined, "terekhova23", {
        "study_celltype": {"source": "Cluster_names"}, "sample": {"source": "Tube_id"},
        "subject": {"source": "Donor_id"}, "sampling_timepoint": {"source": "Visit", "prefix": "visit_"},
        "age": {"source": "Age", "dtype": "float"}, "sex": {"source": "Gender", "lower": True},
        "bmi": {"source": "BMI", "dtype": "float"}, "cmv": {"constant": "not_provided"},
        "ethnicity": {"source": "Ethnicity", "lower": True}, "ethnicity_fine": {"source": "Ethnicity", "lower": True},
        "technology": {"constant": "10X5'v2"}, "aligner": {"constant": "CRv7.0.0"},
        "genome": {"constant": "GRCh38"}, "demultiplexing": {"constant": "genetic"},
        "frozen": {"constant": "yes"}, "include_intronic": {"constant": "yes"},
        "smoking_status": {"constant": "no"}, "study_site": {"constant": "terekhova23"},
        "country": {"constant": "usa"},
        "batch_single_cell": {"source_index": True, "split": {"separator": "_", "index": 0}, "prefix": "terekhova23_"},
    })


def _reshape_visits(table):
    rename = {}
    for visit in (1, 2, 3):
        rename.update({
            f"Visit_{visit}_Tube_id": f"Tube_id_{visit}",
            f"Visit_{visit}_Age": f"Age_{visit}",
            f"Visit_{visit}_BMI": f"BMI_{visit}",
        })
    table = table.rename(columns=rename)
    return pd.wide_to_long(
        table, stubnames=["Tube_id", "Age", "BMI"],
        i=["Donor_id", "Gender", "Ethnicity", "Age_group"], j="Visit", sep="_", suffix=r"\d+",
    ).reset_index().dropna(subset=["Tube_id"])

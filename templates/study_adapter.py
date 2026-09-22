"""Copy this file to ``src/pbmc_pipeline/studies/<study_id>.py`` and adapt it."""
from __future__ import annotations

from .common import canonical


def prepare_cells(source, root):
    """Return one canonical metadata row for each retained source cell.

    Use ``safe_left_join`` from ``.common`` for a sample-level supplement. If
    intentionally excluding cells (for example, reused public cohorts), filter
    ``source`` here and set ``allow_cell_subset: true`` in the study config.

    Use the most specific reported value for ``technology`` (for example,
    ``10X3'v3`` rather than ``10X3'`` when the reagent-kit revision is known).
    Record the source field and any inferred technical metadata in the study
    config's provenance; retain a run/library identifier in
    ``batch_single_cell`` when the source provides one.
    """
    del root
    return canonical(source, "replace_me", {
        "study_celltype": {"source": "source_celltype_column"},
        "sample": {"source": "source_sample_column"},
        "subject": {"source": "source_subject_column"},
        "sampling_timepoint": {"constant": "not_provided"},
        "age": {"source": "source_age_column", "dtype": "float"},
        "sex": {"source": "source_sex_column", "lower": True},
        "bmi": {"constant": float("nan")},
        "cmv": {"constant": "not_provided"},
        "ethnicity": {"constant": "not_provided"},
        "ethnicity_fine": {"constant": "not_provided"},
        "technology": {"constant": "not_provided"},
        "aligner": {"constant": "not_provided"},
        "genome": {"constant": "not_provided"},
        "demultiplexing": {"constant": "not_provided"},
        "frozen": {"constant": "not_provided"},
        "include_intronic": {"constant": "not_provided"},
        "smoking_status": {"constant": "not_provided"},
        "study_site": {"constant": "not_provided"},
        "country": {"constant": "not_provided"},
        "batch_single_cell": {"constant": "not_provided"},
    })

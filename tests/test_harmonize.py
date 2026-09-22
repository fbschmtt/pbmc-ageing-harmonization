import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from pbmc_pipeline.harmonize import harmonize_study

ROOT = Path(__file__).resolve().parents[1]


def _schema():
    return json.loads((ROOT / "config/harmonized_obs_schema.json").read_text())


def _prepared_metadata(schema, index):
    metadata = {
        column: [0.0 if expected_type == "float" else "not_provided"] * len(index)
        for column, expected_type in schema["required"].items()
    }
    metadata.update({
        "study": ["test"] * len(index), "study_celltype": ["unknown"] * len(index),
        "sample": ["sample_1", "sample_2"], "subject": ["subject_1", "subject_2"],
        "batch_single_cell": ["batch_1", "batch_2"],
    })
    return pd.DataFrame(metadata, index=pd.Index(index, name="cell_id"))


def test_harmonized_output_keeps_only_selected_counts_in_x(tmp_path):
    source = ad.AnnData(
        X=np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell_a", "cell_b"]),
        var=pd.DataFrame(index=["gene_a", "gene_b"]),
    )
    source.raw = source.copy()
    counts = np.array([[1, 2], [3, 4]], dtype=np.int32)
    source.layers["counts"] = counts
    source.layers["normalized"] = source.X.copy()
    expression = tmp_path / "source.h5ad"
    source.write_h5ad(expression)

    schema = _schema()
    prepared_path = tmp_path / "prepared.csv.gz"
    _prepared_metadata(schema, source.obs_names).to_csv(prepared_path, compression="infer")
    output = tmp_path / "harmonized.h5ad"
    harmonize_study(
        ROOT,
        "test",
        {
            "counts_source": "layer:counts", "features": {"duplicate_policy": "error"},
            "preparation": {"allow_cell_subset": False}, "annotation": {"method": "existing"},
        },
        {"processing": {"output_compression": "gzip"}},
        schema,
        input_path=expression,
        prepared_obs_path=prepared_path,
        output_path=output,
        report_path=tmp_path / "report.json",
    )

    result = ad.read_h5ad(output)
    assert np.array_equal(result.X, counts)
    assert result.raw is None
    assert list(result.layers) == []

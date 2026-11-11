"""Check that production Nextflow publication preserves models, diagnostics and reports."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def check_publication(outdir: Path) -> None:
    directory = outdir / "analysis/per_type/test-type/differential_expression"
    expected = (
        "run_metadata.json", "merged.csv", "diagnostics.json",
        "fraction_model_coefficients.tsv", "residual_wgcna/metadata.json",
        "report.html", "executed.ipynb", "report_render.log",
    )
    missing = [name for name in expected if not (directory / name).is_file()]
    if missing:
        raise ValueError(f"DE publication lost required artifacts: {missing}")
    if json.loads((directory / "run_metadata.json").read_text())["cell_type"] != "Test type":
        raise ValueError("Published model metadata differs from the input fixture")
    if "TEST_GENE,0.25" not in (directory / "merged.csv").read_text():
        raise ValueError("Diagnostic/report publication replaced the model table")
    print("DE publication check passed: models, diagnostics, nested artifacts and reports coexist.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    check_publication(args.outdir)


if __name__ == "__main__":
    main()

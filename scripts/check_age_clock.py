"""Check that the synthetic age-clock workflow published its expected files."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pbmc_pipeline.utils import slugify

METHODS = ("elastic_net", "equal_study_elastic_net", "rank_elastic_net")
METHOD_OUTPUTS = (
    "model.joblib",
    "external_predictions.csv",
    "cv_grid_search.csv",
    "cv_by_study.csv",
    "cv_oof_predictions.csv",
    "metadata.json",
    "age_clock.json",
    "executed.ipynb",
    "report.html",
)


def check_age_clock(outdir: Path, cell_type: str) -> None:
    slug = slugify(cell_type)
    clock_root = outdir / "analysis" / "per_type" / slug / "ageing_clock"
    required = [
        clock_root / "age_clock.json",
        clock_root / "executed.ipynb",
        clock_root / "report.html",
        outdir / "index.html",
    ]
    required.extend(
        clock_root / method / filename
        for method in METHODS
        for filename in METHOD_OUTPUTS
    )
    missing = [
        str(path.relative_to(outdir))
        for path in required
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise SystemExit("Missing or empty age-clock outputs:\n- " + "\n- ".join(missing))
    manifest = json.loads((clock_root / "age_clock.json").read_text())
    n_jobs = manifest.get("cv_parallel_jobs")
    if not isinstance(n_jobs, int) or n_jobs < 1:
        raise SystemExit(f"{clock_root / 'age_clock.json'} has invalid cv_parallel_jobs: {n_jobs!r}")
    for method in METHODS:
        method_root = clock_root / method
        metadata = json.loads((method_root / "metadata.json").read_text())
        diagnostics = metadata.get("fit_diagnostics", {})
        if "iteration_limit_reached" not in diagnostics:
            raise SystemExit(f"{method_root / 'metadata.json'} lacks convergence diagnostics")
        cv_by_study = pd.read_csv(method_root / "cv_by_study.csv")
        required_fold_columns = {"n_iterations", "dual_gap", "iteration_limit_reached"}
        if missing_columns := sorted(required_fold_columns - set(cv_by_study.columns)):
            raise SystemExit(
                f"{method_root / 'cv_by_study.csv'} lacks convergence columns: {missing_columns}"
            )
        cv_grid = pd.read_csv(method_root / "cv_grid_search.csv")
        required_grid_columns = {"n_iteration_limit_folds"}
        for study in cv_by_study["held_out_study"].astype(str):
            required_grid_columns.update({
                f"n_iterations_{study}", f"dual_gap_{study}",
                f"iteration_limit_reached_{study}",
            })
        if missing_columns := sorted(required_grid_columns - set(cv_grid.columns)):
            raise SystemExit(
                f"{method_root / 'cv_grid_search.csv'} lacks convergence columns: {missing_columns}"
            )
        report_html = (method_root / "report.html").read_text(errors="replace")
        if "Convergence diagnostics" not in report_html:
            raise SystemExit(f"{method_root / 'report.html'} lacks convergence diagnostics")
    print(f"Age-clock output files are present for {cell_type} ({len(required)} files)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("output/test"))
    parser.add_argument("--cell-type", default="cd14-monocyte")
    args = parser.parse_args()
    check_age_clock(args.outdir, args.cell_type)


if __name__ == "__main__":
    main()

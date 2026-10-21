"""Check that the synthetic age-clock workflow published its expected files."""
from __future__ import annotations

import argparse
from pathlib import Path

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
    print(f"Age-clock output files are present for {cell_type} ({len(required)} files)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("output/test"))
    parser.add_argument("--cell-type", default="cd14-monocyte")
    args = parser.parse_args()
    check_age_clock(args.outdir, args.cell_type)


if __name__ == "__main__":
    main()

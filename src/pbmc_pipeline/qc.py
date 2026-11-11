from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .logging_utils import configure_logging
from .notebook import render_notebook

LOGGER = logging.getLogger(__name__)


def generate_qc_report(
    *, root: Path, input_path: Path, run_report_path: Path, output_dir: Path,
    study: str, template_path: Path | None = None,
) -> Path:
    """Execute the QC template and return the generated HTML path."""
    for path in (input_path, run_report_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    return render_notebook(
        template_path=template_path or root / "reports/qc_report.py",
        output_dir=output_dir, project_root=root,
        environment={"QC_INPUT_H5AD": str(input_path.resolve()),
                     "QC_RUN_REPORT": str(run_report_path.resolve()), "QC_STUDY": study},
        title=f"{study} QC | PBMC ageing",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate an executed HTML QC report")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, required=True, help="Harmonized H5AD")
    parser.add_argument("--run-report", type=Path, required=True, help="Harmonization JSON report")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--template", type=Path, help="Override the QC report template")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO",
        help="Progress-log severity written to standard error (default: INFO)",
    )
    args = parser.parse_args()
    configure_logging(args.log_level)
    root = args.project_root.resolve()
    html = generate_qc_report(
        root=root,
        input_path=args.input,
        run_report_path=args.run_report,
        output_dir=args.output_dir,
        study=args.study,
        template_path=args.template,
    )
    print(html)


if __name__ == "__main__":
    main()

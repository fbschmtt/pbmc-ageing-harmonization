from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

from .logging_utils import configure_logging

LOGGER = logging.getLogger(__name__)


def generate_qc_report(
    *,
    root: Path,
    input_path: Path,
    run_report_path: Path,
    output_dir: Path,
    study: str,
    template_path: Path | None = None,
) -> Path:
    """Execute the QC notebook and return the generated HTML path."""
    input_path = input_path.resolve()
    run_report_path = run_report_path.resolve()
    output_dir = output_dir.resolve()
    template = (template_path or root / "reports" / "qc_report.ipynb").resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"QC input does not exist: {input_path}")
    if not run_report_path.exists():
        raise FileNotFoundError(f"Run report does not exist: {run_report_path}")
    if not template.exists():
        raise FileNotFoundError(f"QC notebook template does not exist: {template}")
    total_started = time.perf_counter()
    LOGGER.info("study=%s step=qc_start input=%s", study, input_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    executed = output_dir / "executed.ipynb"
    html = output_dir / "report.html"
    env = os.environ.copy()
    env.update(
        {
            "QC_INPUT_H5AD": str(input_path),
            "QC_RUN_REPORT": str(run_report_path),
            "QC_STUDY": study,
            "MPLCONFIGDIR": str((root / ".cache" / "matplotlib").resolve()),
        }
    )
    execute_command = [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--execute",
        "--to",
        "notebook",
        "--ExecutePreprocessor.timeout=-1",
        f"--output={executed.name}",
        f"--output-dir={output_dir}",
        str(template),
    ]
    started = time.perf_counter()
    subprocess.run(execute_command, cwd=root, env=env, check=True)
    LOGGER.info(
        "study=%s step=execute_qc_notebook duration_seconds=%.2f path=%s",
        study, time.perf_counter() - started, executed,
    )
    export_command = [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "html",
        f"--output={html.name}",
        f"--output-dir={output_dir}",
        str(executed),
    ]
    started = time.perf_counter()
    subprocess.run(export_command, cwd=root, env=env, check=True)
    LOGGER.info(
        "study=%s step=export_qc_html duration_seconds=%.2f path=%s",
        study, time.perf_counter() - started, html,
    )
    LOGGER.info(
        "study=%s step=qc_complete duration_seconds=%.2f",
        study, time.perf_counter() - total_started,
    )
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate an executed HTML QC report")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, required=True, help="Harmonized H5AD")
    parser.add_argument("--run-report", type=Path, required=True, help="Harmonization JSON report")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--template", type=Path, help="Override the QC notebook template")
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

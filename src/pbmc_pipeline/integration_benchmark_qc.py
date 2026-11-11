"""Render the global integration-benchmark notebook and HTML report."""
from __future__ import annotations

import argparse
from pathlib import Path

from .logging_utils import configure_logging
from .notebook import render_notebook


def generate_integration_benchmark_report(
    *, root: Path, input_path: Path, run_report: Path, output_dir: Path,
    template_path: Path | None = None,
) -> Path:
    return render_notebook(
        template_path=template_path or root / "reports/integration_benchmark_report.py",
        output_dir=output_dir, project_root=root,
        environment={"INTEGRATION_BENCHMARK_H5AD": str(input_path.resolve()),
                     "INTEGRATION_BENCHMARK_REPORT": str(run_report.resolve())},
        title="Integration benchmark | PBMC ageing",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Render the global integration benchmark report")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--run-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()
    configure_logging(args.log_level)
    print(generate_integration_benchmark_report(
        root=args.project_root.resolve(), input_path=args.input, run_report=args.run_report,
        output_dir=args.output_dir, template_path=args.template,
    ))

"""Render the global integration-benchmark notebook and HTML report."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .logging_utils import configure_logging
from .report_toc import populate_html_toc, set_notebook_title


def generate_integration_benchmark_report(
    *, root: Path, input_path: Path, run_report: Path, output_dir: Path,
    template_path: Path | None = None,
) -> Path:
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    template = (template_path or root / "reports" / "integration_benchmark_report.py").resolve()
    env = os.environ.copy()
    env.update({
        "INTEGRATION_BENCHMARK_H5AD": str(input_path.resolve()),
        "INTEGRATION_BENCHMARK_REPORT": str(run_report.resolve()),
        "MPLCONFIGDIR": str((root / ".cache" / "matplotlib").resolve()),
    })
    source_path = str(root / "src")
    env["PYTHONPATH"] = os.pathsep.join(
        path for path in (source_path, env.get("PYTHONPATH", "")) if path
    )
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    with tempfile.TemporaryDirectory(prefix="pbmc-integration-benchmark-notebook-") as temporary_dir:
        temporary_path = Path(temporary_dir)
        jupyter_data = temporary_path / "jupyter-data"
        jupyter_data.mkdir()
        env["JUPYTER_DATA_DIR"] = str(jupyter_data)
        notebook = temporary_path / f"{template.stem}.ipynb"
        subprocess.run([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(notebook), str(template),
        ], cwd=root, env=env, check=True)
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={executed.name}",
            f"--output-dir={output_dir}", str(notebook),
        ], cwd=root, env=env, check=True)
    set_notebook_title(executed, "Integration benchmark | PBMC ageing")
    subprocess.run([
        sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
        "--HTMLExporter.exclude_input=True",
        "--HTMLExporter.exclude_input_prompt=True",
        "--HTMLExporter.exclude_output_prompt=True",
        f"--output={html.name}", f"--output-dir={output_dir}", str(executed),
    ], cwd=root, env=env, check=True)
    populate_html_toc(html)
    return html


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

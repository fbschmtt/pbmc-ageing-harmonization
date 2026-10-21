"""Render the cell-type age-clock summary and diagnostic report."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .report_toc import populate_html_toc, set_notebook_title
from .utils import slugify


def generate_age_clock_report(
    *, report_dir: Path, template_path: Path, project_root: Path,
) -> Path:
    """Execute the report template and publish a standalone HTML report."""
    for path in (report_dir, template_path, project_root):
        if not path.exists():
            raise FileNotFoundError(path)
    if not (report_dir / "age_clock.json").is_file():
        raise FileNotFoundError(report_dir / "age_clock.json")
    report_dir = report_dir.resolve()
    template_path = template_path.resolve()
    project_root = project_root.resolve()
    executed = report_dir / "executed.ipynb"
    html = report_dir / "report.html"
    env = os.environ.copy()
    env["AGE_CLOCK_REPORT_DIR"] = str(report_dir)
    env["MPLCONFIGDIR"] = str((report_dir / ".matplotlib").resolve())
    env["JUPYTER_DATA_DIR"] = str((report_dir / ".jupyter-data").resolve())
    source_path = str(project_root / "src")
    env["PYTHONPATH"] = os.pathsep.join(
        path for path in (source_path, env.get("PYTHONPATH", "")) if path
    )
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["JUPYTER_DATA_DIR"]).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pbmc-age-clock-report-") as temp_dir:
        env["IPYTHONDIR"] = str(Path(temp_dir) / "ipython")
        env["JUPYTER_CONFIG_DIR"] = str(Path(temp_dir) / "jupyter-config")
        env["JUPYTER_RUNTIME_DIR"] = str(Path(temp_dir) / "jupyter-runtime")
        for key in ("IPYTHONDIR", "JUPYTER_CONFIG_DIR", "JUPYTER_RUNTIME_DIR"):
            Path(env[key]).mkdir()
        notebook = Path(temp_dir) / f"{template_path.stem}.ipynb"
        subprocess.run([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(notebook), str(template_path),
        ], cwd=project_root, env=env, check=True)
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={executed.name}",
            "--KernelManager.transport_encryption=required",
            f"--output-dir={report_dir}", str(notebook),
        ], cwd=project_root, env=env, check=True)
        set_notebook_title(executed, "Cell-type transcriptomic age clocks | PBMC ageing")
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
            "--HTMLExporter.exclude_input=True",
            "--HTMLExporter.exclude_input_prompt=True",
            "--HTMLExporter.exclude_output_prompt=True",
            f"--output={html.name}", f"--output-dir={report_dir}", str(executed),
        ], cwd=project_root, env=env, check=True)
    populate_html_toc(html)
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--report-dir", type=Path)
    group.add_argument("--report-root", type=Path)
    parser.add_argument("--cell-type", help="Optional cell-type label or slug when scanning a report root")
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--comparison-template", type=Path)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.report_dir:
        report_dirs = [args.report_dir.resolve()]
    else:
        report_root = args.report_root.resolve()
        manifests = sorted(report_root.glob("*/ageing_clock/*/age_clock.json"))
        if not manifests:
            manifests = sorted(report_root.glob("*/age_clock.json"))
        if args.cell_type:
            requested = args.cell_type.strip()
            matches = []
            for path in manifests:
                manifest = json.loads(path.read_text())
                if any(
                    record["cell_type"] == requested
                    or slugify(record["cell_type"]) == requested
                    for record in manifest.get("cell_types", [])
                ):
                    matches.append(path)
            manifests = matches
        report_dirs = [path.parent for path in manifests]
    if not report_dirs:
        raise FileNotFoundError(f"No per-cell-type age-clock manifests under {args.report_root}")
    for report_dir in report_dirs:
        print(generate_age_clock_report(
            report_dir=report_dir,
            template_path=args.template.resolve(),
            project_root=args.project_root.resolve(),
        ))
    if args.report_root and args.comparison_template:
        for comparison_dir in sorted({directory.parent for directory in report_dirs}):
            print(generate_age_clock_report(
                report_dir=comparison_dir,
                template_path=args.comparison_template.resolve(),
                project_root=args.project_root.resolve(),
            ))


if __name__ == "__main__":
    main()

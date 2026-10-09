"""Render one combined QC document for pseudobulk and optional single-cell merges."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .logging_utils import configure_logging
from .report_toc import populate_html_toc


def generate_merge_qc_report(
    *, root: Path, inputs: list[Path], reports: list[Path], output_dir: Path,
    template_path: Path | None = None, studies_config_path: Path | None = None,
) -> Path:
    reports_by_stem = {report_path.stem: report_path for report_path in reports}
    if len(reports_by_stem) != len(reports):
        raise ValueError("Merge report names must be unique")
    by_kind: dict[str, tuple[Path, Path]] = {}
    for input_path in inputs:
        report_path = reports_by_stem.get(input_path.stem)
        if report_path is None:
            raise ValueError(f"{input_path}: matching merge report is absent")
        document = json.loads(report_path.read_text())
        kind = document.get("kind")
        if kind not in {"pseudobulk_merge", "single_cell_merge"}:
            raise ValueError(f"{report_path}: unsupported merge report kind {kind!r}")
        if kind in by_kind:
            raise ValueError(f"More than one {kind} report supplied")
        by_kind[kind] = (input_path.resolve(), report_path.resolve())
    if "pseudobulk_merge" not in by_kind:
        raise ValueError("Combined merge QC requires a pseudobulk_merge report")

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    template = (template_path or root / "reports" / "merge_qc_report.py").resolve()
    env = os.environ.copy()
    pseudo_input, pseudo_report = by_kind["pseudobulk_merge"]
    env.update({
        "QC_PSEUDOBULK_H5AD": str(pseudo_input),
        "QC_PSEUDOBULK_REPORT": str(pseudo_report),
        "QC_STUDIES_CONFIG": str((studies_config_path or root / "config" / "studies.json").resolve()),
    })
    if "single_cell_merge" in by_kind:
        single_input, single_report = by_kind["single_cell_merge"]
        env.update({"QC_SINGLE_CELL_H5AD": str(single_input), "QC_SINGLE_CELL_REPORT": str(single_report)})
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    with tempfile.TemporaryDirectory(prefix="pbmc-merge-qc-notebook-") as temporary_dir:
        matplotlib_cache = Path(temporary_dir) / "matplotlib"
        matplotlib_cache.mkdir()
        env["MPLCONFIGDIR"] = str(matplotlib_cache)
        materialized_template = Path(temporary_dir) / f"{template.stem}.ipynb"
        subprocess.run([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(materialized_template), str(template),
        ], cwd=root, env=env, check=True)
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={executed.name}",
            f"--output-dir={output_dir}", str(materialized_template),
        ], cwd=root, env=env, check=True)
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
    parser = argparse.ArgumentParser(description="Render one combined cross-study merge QC report")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--run-report", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--studies-config", type=Path)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()
    configure_logging(args.log_level)
    print(generate_merge_qc_report(
        root=args.project_root.resolve(), inputs=args.input, reports=args.run_report,
        output_dir=args.output_dir, template_path=args.template,
        studies_config_path=args.studies_config,
    ))

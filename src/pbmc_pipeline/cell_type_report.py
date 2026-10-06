"""Execute the per-cell-type analysis notebook and export self-contained HTML."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def generate_cell_type_report(
    *, input_path: Path, analysis_path: Path, output_dir: Path, template_path: Path,
    config_path: Path, differential_expression_dir: Path | None = None,
) -> Path:
    for path in (input_path, analysis_path, template_path, config_path):
        if not path.exists():
            raise FileNotFoundError(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "age_trajectory_clusters.csv",
        "age_trajectory_cluster_means.csv",
        "residual_sample_clusters.csv",
        "residual_sample_pc_scores.csv",
    ):
        (output_dir / name).unlink(missing_ok=True)
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    env = os.environ.copy()
    env.update({
        "CELL_TYPE_H5AD": str(input_path.resolve()),
        "CELL_TYPE_ANALYSIS_H5AD": str(analysis_path.resolve()),
        "CELL_TYPE_OUTPUT_DIR": str(output_dir.resolve()),
        "CELL_TYPE_CONFIG": str(config_path.resolve()),
    })
    env.pop("CELL_TYPE_DIFFERENTIAL_EXPRESSION_DIR", None)
    if differential_expression_dir is not None:
        if not differential_expression_dir.is_dir():
            raise NotADirectoryError(differential_expression_dir)
        env["CELL_TYPE_DIFFERENTIAL_EXPRESSION_DIR"] = str(
            differential_expression_dir.resolve()
        )
    with tempfile.TemporaryDirectory(prefix="pbmc-cell-type-notebook-") as temporary_dir:
        materialized_template = Path(temporary_dir) / f"{template_path.stem}.ipynb"
        subprocess.run([sys.executable, "-m", "jupytext", "--to", "notebook", "--output", str(materialized_template), str(template_path.resolve())], check=True, env=env)
        subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook", "--output", executed.name, "--output-dir", str(output_dir), str(materialized_template)], check=True, env=env)
    if executed.is_file():
        import nbformat

        report_metadata_path = output_dir / "report.json"
        report_metadata = json.loads(report_metadata_path.read_text()) if report_metadata_path.is_file() else {}
        notebook = nbformat.read(executed, as_version=4)
        notebook.metadata["title"] = (
            f"{report_metadata.get('cell_type', input_path.stem)} | PBMC ageing"
        )
        nbformat.write(notebook, executed)
    subprocess.run([
        sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
        "--HTMLExporter.exclude_input=True",
        "--HTMLExporter.exclude_input_prompt=True",
        "--HTMLExporter.exclude_output_prompt=True",
        "--output", html.name, "--output-dir", str(output_dir), str(executed),
    ], check=True, env=env)
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a cell-type analysis HTML report")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--differential-expression-dir", type=Path)
    args = parser.parse_args()
    print(generate_cell_type_report(
        input_path=args.input,
        analysis_path=args.analysis,
        output_dir=args.output_dir,
        template_path=args.template,
        config_path=args.config,
        differential_expression_dir=args.differential_expression_dir,
    ))


if __name__ == "__main__":
    main()

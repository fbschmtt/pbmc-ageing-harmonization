"""Render the established pseudobulk-DE portion of the cell-type report."""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from .report_toc import populate_html_toc, set_notebook_title


def _run_report_command(command, *, cwd: Path, env: dict[str, str], log_path: Path) -> None:
    """Record notebook rendering commands and surface a useful failure."""
    result = subprocess.run(
        command, cwd=cwd, env=env, capture_output=True, text=True, check=False,
    )
    with log_path.open("a") as log:
        log.write(f"$ {shlex.join(command)}\n")
        for stream in (result.stdout, result.stderr):
            if stream:
                log.write(stream)
                if not stream.endswith("\n"):
                    log.write("\n")
    if result.returncode:
        if result.stdout:
            sys.stdout.write(result.stdout)
        if result.stderr:
            sys.stderr.write(result.stderr)
        raise subprocess.CalledProcessError(
            result.returncode, command, output=result.stdout, stderr=result.stderr,
        )


def _de_notebook_source(template_path: Path) -> str:
    """Keep the split report's DE cells identical to the original report.

    ``reports/cell_type_report.py`` remains the canonical source for these
    helpers and the pseudobulk-DE section. The bootstrap supplies the state
    that the DE section previously received from the combined notebook.
    """
    source = template_path.read_text(encoding="utf-8")
    setup_end = source.index("primitive = sc.read_h5ad")
    helpers_start = source.index("def display_collapsible_table")
    helpers_end = source.index("def natural_sort_key", helpers_start)
    de_start = source.index("if de_dir:", helpers_end)
    imports = source[:setup_end]
    helpers = source[helpers_start:helpers_end]
    de_section = source[de_start:]
    # The combined report imports Scanpy for its single-cell sections. The DE
    # section never uses it, and importing it in this standalone notebook can
    # trigger an irrelevant Numba cache failure in read-only containers.
    imports = imports.replace("import scanpy as sc\n", "")
    imports = imports.replace(
        "from pbmc_pipeline.cell_type import add_unintegrated_umap\n", ""
    )
    bootstrap = '''# %%
de_dir = Path(os.environ["DE_REPORT_DIR"])
output_dir = Path(os.environ["DE_REPORT_OUTPUT_DIR"])
output_dir.mkdir(parents=True, exist_ok=True)
pipeline = read_json(Path(os.environ["DE_REPORT_CONFIG"]))
trajectory_settings = AgeTrajectorySettings.from_mapping(
    pipeline["differential_expression"]["age_trajectory"]
)
de_metadata_path = de_dir / "run_metadata.json"
de_run_metadata = read_json(de_metadata_path) if de_metadata_path.is_file() else {}
result_record_path = de_dir / "cell_type_result.json"
result_record = read_json(result_record_path) if result_record_path.is_file() else {}
cell_type_name = str(
    de_run_metadata.get("cell_type")
    or result_record.get("cell_type")
    or de_dir.parent.name.replace("-", " ").title()
)
completed_models = [
    model for model in de_run_metadata.get("models", [])
    if model.get("status") == "complete"
]
summary_metrics = [
    ("Pseudobulk samples", str(result_record.get("n_pseudobulks", "—"))),
    ("Completed models", str(len(completed_models))),
    ("Analysis mode", str(de_run_metadata.get("analysis_mode", "unknown")).replace("_", " ").title()),
]
display(HTML(render_report_header(
    title=cell_type_name,
    eyebrow="PBMC ageing · pseudobulk differential expression",
    subtitle=(
        "Sample-level pseudobulk differential-expression models, age-bin trajectories, "
        "and residual structure for this cell type."
    ),
    metrics=summary_metrics,
)))
if de_run_metadata.get("analysis_mode") == "test_only":
    display(Markdown(
        "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
        + de_run_metadata["interpretation_warning"]
    ))

'''
    return imports + bootstrap + helpers + "\n# %%\n" + de_section


def generate_differential_expression_report(
    *, result_dir: Path, output_dir: Path, template_path: Path, config_path: Path,
    project_root: Path,
) -> Path:
    """Execute the original DE report section as a standalone notebook report."""
    for path in (result_dir, template_path, config_path, project_root):
        if not path.exists():
            raise FileNotFoundError(path)
    result_dir = result_dir.resolve()
    output_dir = output_dir.resolve()
    template_path = template_path.resolve()
    config_path = config_path.resolve()
    project_root = project_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "age_trajectory_clusters.csv",
        "age_trajectory_cluster_means.csv",
        "residual_sample_clusters.csv",
        "residual_sample_pc_scores.csv",
    ):
        (output_dir / name).unlink(missing_ok=True)
    render_log = output_dir / "report_render.log"
    render_log.unlink(missing_ok=True)
    output_notebook = output_dir / "executed.ipynb"
    output_html = output_dir / "report.html"
    environment = os.environ.copy()
    environment.update({
        "DE_REPORT_DIR": str(result_dir),
        "DE_REPORT_OUTPUT_DIR": str(output_dir),
        "DE_REPORT_CONFIG": str(config_path),
    })
    source_path = str(project_root / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        value for value in (source_path, environment.get("PYTHONPATH", "")) if value
    )

    with tempfile.TemporaryDirectory(prefix="pbmc-de-notebook-") as temporary_dir:
        temporary_path = Path(temporary_dir)
        matplotlib_cache = temporary_path / "matplotlib"
        jupyter_data = temporary_path / "jupyter-data"
        numba_cache = temporary_path / "numba"
        matplotlib_cache.mkdir()
        jupyter_data.mkdir()
        numba_cache.mkdir()
        environment["MPLCONFIGDIR"] = str(matplotlib_cache)
        environment["JUPYTER_DATA_DIR"] = str(jupyter_data)
        environment["NUMBA_CACHE_DIR"] = str(numba_cache)
        materialized_source = temporary_path / "differential_expression_report.py"
        materialized_source.write_text(_de_notebook_source(template_path), encoding="utf-8")
        materialized_notebook = temporary_path / "differential_expression_report.ipynb"
        _run_report_command([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(materialized_notebook), str(materialized_source),
        ], cwd=project_root, env=environment, log_path=render_log)
        _run_report_command([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={output_notebook.name}",
            f"--output-dir={output_dir}", str(materialized_notebook),
        ], cwd=project_root, env=environment, log_path=render_log)

    import json

    metadata_path = result_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    cell_type = str(metadata.get("cell_type") or result_dir.parent.name)
    set_notebook_title(output_notebook, f"{cell_type} pseudobulk differential expression | PBMC ageing")
    _run_report_command([
        sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
        "--HTMLExporter.exclude_input=True",
        "--HTMLExporter.exclude_input_prompt=True",
        "--HTMLExporter.exclude_output_prompt=True",
        f"--output={output_html.name}", f"--output-dir={output_dir}",
        str(output_notebook),
    ], cwd=project_root, env=environment, log_path=render_log)
    populate_html_toc(output_html)
    return output_html


def render_per_type_reports(
    input_dir: Path, output_root: Path, *, template_path: Path, config_path: Path,
    project_root: Path,
) -> list[Path]:
    """Render a standalone DE report for every completed cell-type fit."""
    rendered = []
    for result_dir in sorted(input_dir.glob("*/differential_expression")):
        if not (result_dir / "run_metadata.json").is_file():
            continue
        output_dir = output_root / result_dir.parent.name / "differential_expression"
        rendered.append(generate_differential_expression_report(
            result_dir=result_dir,
            output_dir=output_dir,
            template_path=template_path,
            config_path=config_path,
            project_root=project_root,
        ))
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    reports = render_per_type_reports(
        args.input_dir,
        args.output_root,
        template_path=args.template,
        config_path=args.config,
        project_root=args.project_root,
    )
    print(f"Rendered {len(reports)} per-type DE reports")


if __name__ == "__main__":
    main()

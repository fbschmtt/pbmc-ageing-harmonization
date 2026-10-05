"""Prepare cross-cell-type trajectory clusters and render their notebook report."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .config import AgeTrajectorySettings, read_json
from .trajectory_analysis import analyze_age_trajectories


def generate_trajectory_report(
    *,
    differential_expression_dir: Path,
    output_dir: Path,
    config_path: Path,
    template_path: Path,
    project_root: Path,
) -> Path:
    """Cluster shared-bin trajectories, then execute their HTML notebook."""
    for path in (differential_expression_dir, config_path, template_path):
        if not path.exists():
            raise FileNotFoundError(path)
    pipeline = read_json(config_path)
    settings = AgeTrajectorySettings.from_mapping(
        pipeline["differential_expression"]["age_trajectory"]
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str((project_root / ".cache" / "matplotlib").resolve())
    env["JUPYTER_DATA_DIR"] = str((project_root / ".cache" / "jupyter-data").resolve())
    env["NUMBA_CACHE_DIR"] = str((project_root / ".cache" / "numba").resolve())
    env["TRAJECTORY_ANALYSIS_DIR"] = str(output_dir.resolve())
    os.environ.update({
        "MPLCONFIGDIR": env["MPLCONFIGDIR"],
        "JUPYTER_DATA_DIR": env["JUPYTER_DATA_DIR"],
        "NUMBA_CACHE_DIR": env["NUMBA_CACHE_DIR"],
        "TRAJECTORY_ANALYSIS_DIR": env["TRAJECTORY_ANALYSIS_DIR"],
    })
    for key in ("MPLCONFIGDIR", "JUPYTER_DATA_DIR", "NUMBA_CACHE_DIR"):
        Path(env[key]).mkdir(parents=True, exist_ok=True)
    analyze_age_trajectories(
        differential_expression_dir,
        output_dir,
        settings,
    )
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    with tempfile.TemporaryDirectory(prefix="pbmc-trajectory-notebook-") as temp_dir:
        notebook = Path(temp_dir) / f"{template_path.stem}.ipynb"
        subprocess.run([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(notebook), str(template_path.resolve()),
        ], cwd=project_root, env=env, check=True)
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={executed.name}",
            f"--output-dir={output_dir}", str(notebook),
        ], cwd=project_root, env=env, check=True)
    subprocess.run([
        sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
        "--HTMLExporter.exclude_input=True",
        "--HTMLExporter.exclude_input_prompt=True",
        "--HTMLExporter.exclude_output_prompt=True",
        f"--output={html.name}", f"--output-dir={output_dir}", str(executed),
    ], cwd=project_root, env=env, check=True)
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--differential-expression-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    print(generate_trajectory_report(
        differential_expression_dir=args.differential_expression_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        config_path=args.config.resolve(),
        template_path=args.template.resolve(),
        project_root=args.project_root.resolve(),
    ))


if __name__ == "__main__":
    main()

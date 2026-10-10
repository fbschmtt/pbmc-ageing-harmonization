"""Execute and render the focused gene-profile notebook report."""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .report_toc import populate_html_toc, set_notebook_title


def generate_gene_profile_report(
    *, profile_dir: Path, template_path: Path, project_root: Path,
) -> Path:
    """Execute the gene-profile notebook and build its heading-based HTML TOC."""
    for path in (profile_dir / "fit_summary.json", template_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    profile_dir = profile_dir.resolve()
    project_root = project_root.resolve()
    output_notebook = profile_dir / "executed.ipynb"
    output_html = profile_dir / "report.html"
    environment = os.environ.copy()
    environment["GENE_PROFILE_DIR"] = str(profile_dir)
    source_path = str(project_root / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        path for path in (source_path, environment.get("PYTHONPATH", "")) if path
    )

    with tempfile.TemporaryDirectory(prefix="pbmc-gene-profile-notebook-") as temp_dir:
        temp_path = Path(temp_dir)
        matplotlib_cache = temp_path / "matplotlib"
        jupyter_data = temp_path / "jupyter-data"
        matplotlib_cache.mkdir()
        jupyter_data.mkdir()
        environment["MPLCONFIGDIR"] = str(matplotlib_cache)
        environment["JUPYTER_DATA_DIR"] = str(jupyter_data)
        materialized_template = temp_path / f"{template_path.stem}.ipynb"
        subprocess.run([
            sys.executable, "-m", "jupytext", "--to", "notebook",
            "--output", str(materialized_template), str(template_path.resolve()),
        ], cwd=project_root, env=environment, check=True)
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook",
            "--ExecutePreprocessor.timeout=-1", f"--output={output_notebook.name}",
            f"--output-dir={profile_dir}", str(materialized_template),
        ], cwd=project_root, env=environment, check=True)
        import json

        summary = json.loads((profile_dir / "fit_summary.json").read_text())
        set_notebook_title(
            output_notebook,
            f"{summary['gene']} · {summary['cell_type']} | PBMC ageing",
        )
        subprocess.run([
            sys.executable, "-m", "jupyter", "nbconvert", "--to", "html",
            "--HTMLExporter.exclude_input=True",
            "--HTMLExporter.exclude_input_prompt=True",
            "--HTMLExporter.exclude_output_prompt=True",
            f"--output={output_html.name}", f"--output-dir={profile_dir}",
            str(output_notebook),
        ], cwd=project_root, env=environment, check=True)
    populate_html_toc(output_html, minimum_heading_level=1)
    return output_html

"""Execute the per-cell-type analysis notebook and export self-contained HTML."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def generate_cell_type_report(*, input_path: Path, output_dir: Path, template_path: Path, config_path: Path) -> Path:
    for path in (input_path, template_path, config_path):
        if not path.exists():
            raise FileNotFoundError(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    executed, html = output_dir / "executed.ipynb", output_dir / "report.html"
    env = os.environ.copy()
    env.update({
        "CELL_TYPE_H5AD": str(input_path.resolve()),
        "CELL_TYPE_OUTPUT_DIR": str(output_dir.resolve()),
        "CELL_TYPE_CONFIG": str(config_path.resolve()),
    })
    with tempfile.TemporaryDirectory(prefix="pbmc-cell-type-notebook-") as temporary_dir:
        materialized_template = Path(temporary_dir) / f"{template_path.stem}.ipynb"
        subprocess.run([sys.executable, "-m", "jupytext", "--to", "notebook", "--output", str(materialized_template), str(template_path.resolve())], check=True, env=env)
        subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--execute", "--to", "notebook", "--output", executed.name, "--output-dir", str(output_dir), str(materialized_template)], check=True, env=env)
    subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "html", "--output", html.name, "--output-dir", str(output_dir), str(executed)], check=True, env=env)
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a cell-type analysis HTML report")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    print(generate_cell_type_report(input_path=args.input, output_dir=args.output_dir, template_path=args.template, config_path=args.config))


if __name__ == "__main__":
    main()

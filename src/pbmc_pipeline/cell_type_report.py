"""Execute the per-cell-type analysis notebook and export self-contained HTML."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from pbmc_pipeline.config import load_analysis_configuration
from pbmc_pipeline.notebook import render_notebook


def generate_cell_type_report(
    *, input_path: Path, analysis_path: Path, output_dir: Path, template_path: Path,
    config_path: Path,
) -> Path:
    for path in (input_path, analysis_path, template_path, config_path):
        if not path.exists():
            raise FileNotFoundError(path)
    load_analysis_configuration(config_path, sections=("cell_type_analysis", "processing"))
    def title():
        metadata_path = output_dir / "report.json"
        metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
        return f"{metadata.get('cell_type', input_path.stem)} | PBMC ageing"
    return render_notebook(
        template_path=template_path, output_dir=output_dir, project_root=Path.cwd(),
        environment={"CELL_TYPE_H5AD": str(input_path.resolve()),
                     "CELL_TYPE_ANALYSIS_H5AD": str(analysis_path.resolve()),
                     "CELL_TYPE_OUTPUT_DIR": str(output_dir.resolve()),
                     "CELL_TYPE_CONFIG": str(config_path.resolve())}, title=title,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a cell-type analysis HTML report")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    print(generate_cell_type_report(
        input_path=args.input,
        analysis_path=args.analysis,
        output_dir=args.output_dir,
        template_path=args.template,
        config_path=args.config,
    ))


if __name__ == "__main__":
    main()

"""Execute and export the per-cell-type pseudobulk DE report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_analysis_configuration
from .de_diagnostics import copy_diagnostics, validate_diagnostic_sources
from .notebook import render_notebook


def generate_differential_expression_report(
    *, result_dir: Path, output_dir: Path, template_path: Path, config_path: Path,
    project_root: Path, analysis_dir: Path,
) -> Path:
    """Render prepared DE diagnostics without fitting or publishing analysis tables."""
    for path in (result_dir, template_path, config_path, project_root):
        if not path.exists():
            raise FileNotFoundError(path)
    validate_diagnostic_sources(result_dir, analysis_dir)
    copy_diagnostics(analysis_dir, output_dir, load_analysis_configuration(config_path))
    metadata_path = result_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    cell_type = str(metadata.get("cell_type") or result_dir.parent.name)
    return render_notebook(
        template_path=template_path, output_dir=output_dir, project_root=project_root,
        environment={
            "DE_REPORT_DIR": str(result_dir.resolve()),
            "DE_REPORT_OUTPUT_DIR": str(output_dir.resolve()),
            "DE_REPORT_CONFIG": str(config_path.resolve()),
            "DE_REPORT_ANALYSIS_DIR": str(output_dir.resolve()),
        },
        title=f"{cell_type} pseudobulk differential expression | PBMC ageing",
        log_path=output_dir / "report_render.log", stream=True,
        progress_variable="DE_REPORT_PROGRESS_LOG",
    )


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
            analysis_dir=result_dir,
        ))
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description="Render per-cell-type pseudobulk DE reports")
    input_mode = parser.add_mutually_exclusive_group(required=True)
    input_mode.add_argument(
        "--input-dir", type=Path,
        help="render all completed per-type fits found under this directory",
    )
    input_mode.add_argument(
        "--result-dir", type=Path,
        help="render one completed per-type fit directory",
    )
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--analysis-dir", type=Path, help="prepared diagnostic artifacts for one fit")
    args = parser.parse_args()
    if args.result_dir is not None:
        if args.output_dir is None:
            parser.error("--result-dir requires --output-dir")
        if args.analysis_dir is None:
            parser.error("--result-dir requires --analysis-dir from pbmc-de-diagnostics")
        report = generate_differential_expression_report(
            result_dir=args.result_dir,
            output_dir=args.output_dir,
            template_path=args.template,
            config_path=args.config,
            project_root=args.project_root,
            analysis_dir=args.analysis_dir,
        )
        print(report)
        return
    if args.output_root is None:
        parser.error("--input-dir requires --output-root")
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

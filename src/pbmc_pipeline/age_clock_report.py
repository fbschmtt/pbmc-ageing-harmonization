"""Render the cell-type age-clock summary and diagnostic report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .notebook import render_notebook
from .utils import slugify


def generate_age_clock_report(
    *, report_dir: Path, template_path: Path, project_root: Path,
) -> Path:
    """Render one clock method or the cell-type method comparison."""
    manifest = json.loads((report_dir / "age_clock.json").read_text())
    record = manifest.get("cell_types", [{}])[0]
    return render_notebook(
        template_path=template_path, output_dir=report_dir, project_root=project_root,
        environment={"AGE_CLOCK_REPORT_DIR": str(report_dir.resolve())},
        title=f"{record.get('cell_type', 'cell type')} · "
              f"{record.get('method_name', 'Transcriptomic')} age clock | PBMC ageing",
    )


def generate_cross_type_age_clock_report(
    *, report_dir: Path, template_path: Path, project_root: Path,
) -> Path:
    """Render the cross-cell-type clock summary when an aggregate is available."""
    summary = json.loads((report_dir / "summary.json").read_text())
    if not summary.get("methods"):
        raise ValueError(f"{report_dir}: cross-cell-type clock summary has no evaluated methods")
    return render_notebook(
        template_path=template_path, output_dir=report_dir, project_root=project_root,
        environment={"AGE_CLOCK_CROSS_TYPE_DIR": str(report_dir.resolve())},
        title="Cross-cell-type transcriptomic clocks | PBMC ageing",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--report-dir", type=Path)
    group.add_argument("--report-root", type=Path)
    parser.add_argument("--cell-type", help="Optional cell-type label or slug when scanning a report root")
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--comparison-template", type=Path)
    parser.add_argument("--cross-type-template", type=Path)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.report_dir:
        report_dirs = [args.report_dir.resolve()]
    else:
        report_root = args.report_root.resolve()
        type_manifests = sorted(report_root.glob("*/ageing_clock/age_clock.json"))
        manifests = []
        if type_manifests:
            # The type-level manifest is authoritative; stale method folders can
            # remain after an input change or a focused refit.
            for type_manifest in type_manifests:
                type_document = json.loads(type_manifest.read_text())
                if args.cell_type:
                    requested = args.cell_type.strip()
                    if type_document.get("cell_type") != requested and slugify(
                        str(type_document.get("cell_type", ""))
                    ) != requested:
                        continue
                manifests.extend(
                    type_manifest.parent / str(method["method"]) / "age_clock.json"
                    for method in type_document.get("methods", [])
                    if isinstance(method, dict)
                    and method.get("method")
                    and (type_manifest.parent / str(method["method"]) / "age_clock.json").is_file()
                )
        else:
            manifests = sorted(report_root.glob("*/ageing_clock/*/age_clock.json"))
            if not manifests:
                manifests = sorted(report_root.glob("*/age_clock.json"))
        if args.cell_type and not type_manifests:
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
    if args.report_root and args.cross_type_template:
        cross_type_dir = args.report_root.resolve() / "age_clock_cross_type"
        if (cross_type_dir / "summary.json").is_file():
            print(generate_cross_type_age_clock_report(
                report_dir=cross_type_dir,
                template_path=args.cross_type_template.resolve(),
                project_root=args.project_root.resolve(),
            ))


if __name__ == "__main__":
    main()

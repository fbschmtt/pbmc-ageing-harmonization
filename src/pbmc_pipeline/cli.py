from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

from .config import load_configuration
from .harmonize import harmonize_study
from .logging_utils import configure_logging
from .qc import generate_qc_report

LOGGER = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description="Harmonize one or all configured PBMC studies")
    parser.add_argument("--study", required=True, help="Study ID or 'all'")
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, help="Override the configured input (one study only)")
    parser.add_argument("--output", type=Path, help="Override the H5AD output (one study only)")
    parser.add_argument("--report-output", type=Path, help="Override the JSON report (one study only)")
    parser.add_argument("--test", action="store_true", help="Use downsampled test inputs")
    parser.add_argument("--validate-only", action="store_true", help="Skip CellTypist and H5AD output")
    parser.add_argument("--qc", action="store_true", help="Execute the QC notebook after writing H5AD")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO",
        help="Progress-log severity written to standard error (default: INFO)",
    )
    args = parser.parse_args()
    if args.qc and args.validate_only:
        parser.error("--qc cannot be combined with --validate-only")
    if args.study == "all" and any((args.input, args.output, args.report_output)):
        parser.error("path overrides require exactly one --study")
    if args.validate_only and args.output:
        parser.error("--output cannot be combined with --validate-only")

    configure_logging(args.log_level)
    root = args.project_root.resolve()
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".cache" / "matplotlib"))
    os.environ.setdefault("CELLTYPIST_FOLDER", str(root / ".cache" / "celltypist"))
    pipeline, studies_document, schema = load_configuration(root, args.config)
    studies = studies_document["studies"]
    selected = list(studies) if args.study == "all" else [args.study]
    unknown = set(selected) - set(studies)
    if unknown:
        parser.error(f"unknown study: {', '.join(sorted(unknown))}")
    for study_id in selected:
        try:
            report = harmonize_study(
                root, study_id, studies[study_id], pipeline, schema,
                test=args.test, validate_only=args.validate_only,
                input_path=args.input,
                output_path=args.output,
                report_path=args.report_output,
            )
        except Exception:
            LOGGER.exception("study=%s step=failed", study_id)
            raise
        print(json.dumps({"study": study_id, "status": report["status"], "cells": report["n_cells"]}))
        if args.qc:
            suffix = ".test" if args.test else ""
            output_path = args.output or root / pipeline["output_dir"] / f"{study_id}{suffix}.h5ad"
            report_path = (args.report_output or
                           root / pipeline["report_dir"] / f"{study_id}{suffix}.json")
            html = generate_qc_report(
                root=root,
                input_path=output_path,
                run_report_path=report_path,
                output_dir=root / "output" / "qc" / f"{study_id}{suffix}",
                study=study_id,
            )
            print(json.dumps({"study": study_id, "qc_html": str(html)}))

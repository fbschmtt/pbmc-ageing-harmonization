from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .config import load_configuration
from .harmonize import harmonize_study
from .qc import generate_qc_report


def main() -> None:
    parser = argparse.ArgumentParser(description="Harmonize one or all configured PBMC studies")
    parser.add_argument("--study", required=True, help="Study ID or 'all'")
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--test", action="store_true", help="Use downsampled test inputs")
    parser.add_argument("--validate-only", action="store_true", help="Skip CellTypist and H5AD output")
    parser.add_argument("--qc", action="store_true", help="Execute the QC notebook after writing H5AD")
    args = parser.parse_args()
    if args.qc and args.validate_only:
        parser.error("--qc cannot be combined with --validate-only")

    root = Path.cwd()
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".cache" / "matplotlib"))
    os.environ.setdefault("CELLTYPIST_FOLDER", str(root / ".cache" / "celltypist"))
    pipeline, studies_document, schema = load_configuration(root, args.config)
    studies = studies_document["studies"]
    selected = list(studies) if args.study == "all" else [args.study]
    unknown = set(selected) - set(studies)
    if unknown:
        parser.error(f"unknown study: {', '.join(sorted(unknown))}")
    for study_id in selected:
        report = harmonize_study(
            root, study_id, studies[study_id], pipeline, schema,
            test=args.test, validate_only=args.validate_only,
        )
        print(json.dumps({"study": study_id, "status": report["status"], "cells": report["n_cells"]}))
        if args.qc:
            suffix = ".test" if args.test else ""
            html = generate_qc_report(
                root=root,
                input_path=root / pipeline["output_dir"] / f"{study_id}{suffix}.h5ad",
                run_report_path=root / pipeline["report_dir"] / f"{study_id}{suffix}.json",
                output_dir=root / "output" / "qc" / f"{study_id}{suffix}",
                study=study_id,
            )
            print(json.dumps({"study": study_id, "qc_html": str(html)}))

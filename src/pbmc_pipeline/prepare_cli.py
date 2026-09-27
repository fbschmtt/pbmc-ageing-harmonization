from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_configuration
from .preparation import prepare_study


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare one canonical metadata row per expression cell")
    parser.add_argument("--study", required=True)
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--input", type=Path, help="Expression H5AD to index")
    parser.add_argument("--test", action="store_true", help="Use the configured test H5AD")
    parser.add_argument("--output", type=Path, required=True, help="Prepared .csv.gz table")
    parser.add_argument(
        "--source-obs-output", type=Path, required=True,
        help="Retained source obs sidecar (.csv.gz) for adapter auditability",
    )
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    pipeline, document, schema = load_configuration(root, args.config)
    if args.study not in document["studies"]:
        parser.error(f"unknown study: {args.study}")
    if args.input is None:
        if not args.test:
            parser.error("--input is required unless --test is supplied")
        input_path = root / pipeline["test_input_root"] / document["studies"][args.study]["test_input"]
    else:
        input_path = args.input.resolve()
    report = prepare_study(
        root, args.study, document["studies"][args.study], schema,
        input_path, args.output, args.report_output,
        source_obs_output_path=args.source_obs_output,
    )
    print(json.dumps({"study": args.study, "cells": report["n_cells"], "status": "prepared"}))

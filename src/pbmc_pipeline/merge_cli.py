from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import read_json
from .logging_utils import configure_logging
from .merge import merge_pseudobulks, merge_single_cells, pseudobulk_study


def main() -> None:
    parser = argparse.ArgumentParser(description="Create and merge PBMC study pseudobulks or single cells")
    parser.add_argument("--mode", choices=("pseudobulk", "pseudobulk-merge", "single-cell-merge"), required=True)
    parser.add_argument("--input", type=Path, required=True, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=Path("config/pipeline.json"))
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()
    configure_logging(args.log_level)
    root = args.project_root.resolve()
    pipeline = read_json(root / args.config)
    if args.mode == "pseudobulk":
        if len(args.input) != 1:
            parser.error("pseudobulk mode requires exactly one --input")
        report = pseudobulk_study(args.input[0], args.output, args.report_output, pipeline)
    elif args.mode == "pseudobulk-merge":
        report = merge_pseudobulks(args.input, args.output, args.report_output, pipeline)
    else:
        report = merge_single_cells(args.input, args.output, args.report_output, pipeline)
    print(json.dumps({"kind": report["kind"], "status": report["status"], "output": str(args.output)}))

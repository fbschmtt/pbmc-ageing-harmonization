"""Inspect a saved workflow without rerunning tasks or changing its outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.workflow_diagnostics import format_failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=Path('output'))
    parser.add_argument('--summary', type=Path, help='Archived workflow summary instead of the latest one')
    args = parser.parse_args()
    path = args.summary if args.summary is not None else args.outdir / 'workflow_summary.json'
    summary = json.loads(path.read_text())
    print(f"Workflow {summary['status']}: {path}")
    print(format_failures(summary))


if __name__ == '__main__':
    main()

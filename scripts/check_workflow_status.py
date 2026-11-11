"""Require a complete current workflow before accepting existing scientific artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def require_complete(outdir: Path) -> dict:
    path = outdir / 'workflow_summary.json'
    summary = json.loads(path.read_text())
    if summary['status'] != 'complete' or not summary['tasks']:
        raise SystemExit(f'Integration workflow is incomplete; see {path}')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    require_complete(parser.parse_args().outdir)


if __name__ == '__main__':
    main()

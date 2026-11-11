"""Validate the completed composed integration suite without rerunning workflows."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from check_age_clock import check_age_clock
from check_cell_type_reports import check_reports as check_cell_types
from check_gene_profile import check_profile
from check_synthetic_de import check_reports, check_results
from check_workflow_status import require_complete


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--cell-type', default='cd14-monocyte')
    parser.add_argument('--cell-types', required=True)
    parser.add_argument('--methods', required=True)
    args = parser.parse_args()
    require_complete(args.outdir)
    subprocess.run([sys.executable, 'scripts/check_integration_benchmark.py', '--outdir', str(args.outdir)], check=True)
    check_results(args.outdir)
    check_age_clock(args.outdir, args.cell_type, tuple(args.methods.split(',')))
    check_cell_types(args.outdir, cell_types=args.cell_types)
    check_reports(args.outdir)
    check_profile(args.outdir / 'analysis/per_gene/cd14-monocyte/synth-age-marker-1', 'SYNTH_AGE_MARKER_1')
    print('Composed workflow integration checks passed.')


if __name__ == '__main__':
    main()

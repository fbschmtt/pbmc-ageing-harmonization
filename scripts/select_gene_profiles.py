"""Select profile genes from one explicit current-run DE result directory."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from profile_top_age_genes import _top_genes_by_cell_type


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result-dir', type=Path, required=True)
    parser.add_argument('--top-n', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error('--top-n must be positive')
    # Reuse the established ranking with an explicit path, avoiding directory scans.
    selected = _top_genes_by_cell_type(args.result_dir, args.top_n, result_paths=[args.result_dir / 'merged.csv'])
    with args.output.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerows(selected)


if __name__ == '__main__':
    main()

"""Select profile genes from one explicit current-run DE result directory."""
from __future__ import annotations

import argparse
import csv
import json
import sys
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
    result_path = args.result_dir / 'merged.csv'
    metadata_path = args.result_dir / 'run_metadata.json'
    skipped_model = None
    if not result_path.exists() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text())
        merged_models = [model for model in metadata['models'] if model['model'] == 'merged']
        if len(merged_models) == 1 and merged_models[0]['status'] == 'skipped':
            skipped_model = merged_models[0]
    if skipped_model is not None:
        print(f"Skipping gene selection for {metadata['cell_type']}: {skipped_model['reason']}", file=sys.stderr)
        selected = []
    else:
        # Reuse the established ranking with an explicit path, avoiding directory scans.
        selected = _top_genes_by_cell_type(args.result_dir, args.top_n, result_paths=[result_path])
    with args.output.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerows(selected)


if __name__ == '__main__':
    main()

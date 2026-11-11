"""List clock types without requiring DE model or trajectory settings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from pbmc_pipeline.config import load_analysis_configuration
from pbmc_pipeline.differential_expression import _cell_type_index


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input', 'config', 'output', 'expected'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--selection', default='')
    args = parser.parse_args()
    pipeline = load_analysis_configuration(args.config, sections=('age_clock',))
    adata = ad.read_h5ad(args.input, backed='r')
    try:
        types = _cell_type_index(adata.obs[pipeline['differential_expression']['split_by']])
    finally:
        adata.file.close()
    args.selection = args.selection.strip()
    if args.selection:
        matches = [(slug, label) for slug, label in types if args.selection in (slug, label)]
        if not matches:
            raise ValueError(f'Unknown clock cell type {args.selection!r}. Available: {types}')
        types = matches
    if not types:
        raise ValueError('Pseudobulk input contains no cell types')
    args.output.write_text(''.join(f'{slug}\t{label}\n' for slug, label in types))
    args.expected.write_text(json.dumps([label for _, label in types]) + '\n')


if __name__ == '__main__':
    main()

"""Assert that failed analysis tasks leave successful siblings and partial summaries."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def check_failures(outdir: Path, phase: str = 'fits') -> None:
    required = [
        'analysis/per_gene/test-type/test-gene/gene_age_result.csv',
        'analysis/per_gene/test-type/test-gene/report.html',
        'analysis/per_type/test-type/ageing_clock/age_clock.json',
        'analysis/per_type/test-type/ageing_clock/report.html',
        'analysis/per_type/test-type/differential_expression/merged.csv',
        'analysis/per_type/test-type/differential_expression/report.html',
        'analysis/per_type/cell_type_manifest.json', 'index.html',
    ]
    if phase != 'all-cells':
        required.append('analysis/per_type/test-type/analysis/report.html')
    missing = [name for name in required if not (outdir / name).is_file()]
    assert not missing, missing
    summaries = list((outdir / 'provenance/workflows').glob('*/summary.json'))
    assert len(summaries) == 1, summaries
    summary = json.loads(summaries[0].read_text())
    assert summary['exitcode'] == 0 and summary['status'] == 'partial', summary
    failed_names = {row['name'].split(':')[-1].split(' (')[0] for row in summary['failed_tasks']}
    expected = {
        'fits': {'FIT_GENE_PROFILE', 'FIT_AGE_CLOCK', 'ANALYSE_CELL_TYPE', 'RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION'},
        'reports': {'RENDER_GENE_PROFILE', 'RENDER_AGE_CLOCK', 'RENDER_CELL_TYPE_REPORT', 'RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS'},
        'all-cells': {'ANALYSE_CELL_TYPE'},
    }[phase]
    assert failed_names == expected, failed_names
    if phase == 'reports':
        assert (outdir / 'analysis/per_type/failed-type/analysis/analysis.h5ad').is_file()
        assert (outdir / 'analysis/per_gene/failed-type/test-gene/gene_age_result.csv').is_file()
        assert (outdir / 'analysis/per_type/failed-type/ageing_clock/age_clock.json').is_file()
        assert (outdir / 'analysis/per_type/failed-type/differential_expression/merged.csv').is_file()
    manifest = json.loads((outdir / 'analysis/per_type/cell_type_manifest.json').read_text())
    assert manifest['status'] == 'partial', manifest
    if phase == 'all-cells':
        assert manifest['cell_types'] == []
        assert len(manifest['failed_cell_types']) == 2
    else:
        assert [r['slug'] for r in manifest['cell_types']] == ['test-type'], manifest
        assert manifest['failed_cell_types'] == [
            {'slug': 'failed-type', 'cell_type': 'Failed type', 'failed_stage': 'report' if phase == 'reports' else 'analysis'},
        ], manifest
    print(f'Failure isolation passed ({phase}): successful siblings, preserved data and partial summaries.')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--phase', choices=('fits', 'reports', 'all-cells'), default='fits')
    args = parser.parse_args()
    check_failures(args.outdir, args.phase)


if __name__ == '__main__':
    main()

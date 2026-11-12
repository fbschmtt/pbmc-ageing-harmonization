"""Validate current composed test tasks and every selected gene profile without fitting."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ANALYSIS_PROCESSES = frozenset({
    'VALIDATE_CONFIGURATION',
    'GENERATE_SYNTHETIC_PSEUDOBULK',
    'LIST_PSEUDOBULK_CELL_TYPES',
    'RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION',
    'WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST',
    'ANALYSE_CELL_TYPE_DE_DIAGNOSTICS',
    'RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS',
    'RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT',
    'COLLECT_PROFILE_DE_CONTEXT',
    'SELECT_PROFILE_GENES',
    'PREPARE_GENE_PROFILE_MODEL',
    'FIT_GENE_PROFILE',
    'RENDER_GENE_PROFILE',
    'SPLIT_CELL_TYPES',
    'ANALYSE_CELL_TYPE',
    'RENDER_CELL_TYPE_REPORT',
    'WRITE_CELL_TYPE_MANIFEST',
    'LIST_CLOCK_CELL_TYPES',
    'FIT_AGE_CLOCK',
    'RENDER_AGE_CLOCK',
    'AGGREGATE_AGE_CLOCKS',
    'RUN_INTEGRATION_METHOD',
    'AGGREGATE_INTEGRATION_METHODS',
    'RENDER_INTEGRATION_BENCHMARK_REPORT',
})
CORE_PROCESSES = frozenset({
    'PREPARE_CELLS',
    'HARMONIZE',
    'PSEUDOBULK',
    'MERGE_PSEUDOBULKS',
    'MERGE_SINGLE_CELLS',
    'WRITE_RUN_MANIFEST',
    'RENDER_STUDY_QC',
    'RENDER_MERGE_QC',
})


def task_identity(task: dict) -> tuple[str, str]:
    name, separator, tag = task['name'].partition(' (')
    return name.split(':')[-1], tag.removesuffix(')') if separator else ''


def task_index(tasks: list[dict], process: str) -> dict[str, dict]:
    records = [task for task in tasks if task_identity(task)[0] == process]
    tags = [task_identity(task)[1] for task in records]
    if len(tags) != len(set(tags)):
        raise ValueError(f'Duplicate {process} tasks: {tags}')
    return dict(zip(tags, records, strict=True))


def require_branches(summary: dict, mode: str) -> list[dict]:
    if mode not in {'full', 'analysis'}:
        raise ValueError(f'Unsupported test mode: {mode}')
    tasks = summary.get('tasks', [])
    if summary.get('status') != 'complete' or not tasks or summary.get('exitcode') != 0:
        raise ValueError('Expected a complete current workflow with successful tasks')
    for task in tasks:
        if task.get('status') not in {'COMPLETED', 'CACHED'} or str(task.get('exit')) != '0':
            raise ValueError(f'Unsuccessful current task: {task}')
    found = {task_identity(task)[0] for task in tasks}
    required = ANALYSIS_PROCESSES | (CORE_PROCESSES if mode == 'full' else set())
    if missing := required - found:
        raise ValueError(f'Missing current workflow processes: {sorted(missing)}')
    if mode == 'analysis' and found & CORE_PROCESSES:
        raise ValueError('Analysis mode unexpectedly executed core processing')
    return tasks


def _same_file(source: Path, published: Path) -> None:
    if not source.is_file() or not published.is_file():
        raise ValueError(f'Missing current/published artifact: {source} -> {published}')
    if hashlib.sha256(source.read_bytes()).digest() != hashlib.sha256(published.read_bytes()).digest():
        raise ValueError(f'Published artifact differs from current task: {published}')


def check_gene_profiles(tasks: list[dict], outdir: Path, *, planted_gene: bool = True) -> int:
    requests = {}
    for slug, task in task_index(tasks, 'SELECT_PROFILE_GENES').items():
        path = Path(task['workdir']) / 'genes.tsv'
        with path.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if not row:
                    continue
                if len(row) != 2 or not all(row):
                    raise ValueError(f'Invalid selected gene row in {path}: {row}')
                cell_type, gene = row
                gene_slug = re.sub(r'[^a-z0-9]+', '-', gene.lower()).strip('-') or 'unnamed'
                key = f'{slug}/{gene_slug}'
                if key in requests:
                    raise ValueError(f'Duplicate or colliding selected gene: {key}')
                requests[key] = (cell_type, gene)
    if not requests:
        raise ValueError('No ranked gene profiles were requested')
    if planted_gene:
        key = 'cd14-monocyte/synth-age-marker-1'
        expected = ('CD14 monocyte', 'SYNTH_AGE_MARKER_1')
        if key in requests and requests[key] != expected:
            raise ValueError(f'Planted gene slug collision: {key}')
        requests[key] = expected
    stages = {name: task_index(tasks, name) for name in ('FIT_GENE_PROFILE', 'RENDER_GENE_PROFILE')}
    for name, records in stages.items():
        if set(records) != set(requests):
            raise ValueError(f'{name} requests differ: missing={sorted(set(requests)-set(records))}, '
                             f'unexpected={sorted(set(records)-set(requests))}')
    prepared = task_index(tasks, 'PREPARE_GENE_PROFILE_MODEL')
    if set(prepared) != {key.split('/')[0] for key in requests}:
        raise ValueError('Model preparations do not match requested profile cell types')
    for key, (cell_type, gene) in requests.items():
        for name, files in (
            ('FIT_GENE_PROFILE', ('fit_summary.json', 'gene_age_result.csv')),
            ('RENDER_GENE_PROFILE', ('report.html', 'executed.ipynb', 'report_render.log')),
        ):
            directory = Path(stages[name][key]['workdir']) / 'profiles' / key
            for filename in files:
                _same_file(directory / filename, outdir / 'analysis/per_gene' / key / filename)
            if name == 'FIT_GENE_PROFILE':
                document = json.loads((directory / 'fit_summary.json').read_text())
                if (document.get('cell_type'), document.get('gene')) != (cell_type, gene):
                    raise ValueError(f'Profile identity differs from selection: {key}')
    return len(requests)


def check_inventory(outdir: Path, *, mode: str, cell_types: str, clock_cell_type: str) -> None:
    summary = json.loads((outdir / 'workflow_summary.json').read_text())
    tasks = require_branches(summary, mode)
    for task in tasks:
        if (Path(task['workdir']) / '.exitcode').read_text().strip() != '0':
            raise ValueError(f'Task exit code disagrees with trace: {task["name"]}')
    selected = {value.strip() for value in cell_types.split(',') if value.strip()}
    for process, expected in (('ANALYSE_CELL_TYPE', selected), ('RENDER_CELL_TYPE_REPORT', selected),
                              ('FIT_AGE_CLOCK', {clock_cell_type}), ('RENDER_AGE_CLOCK', {clock_cell_type})):
        actual = set(task_index(tasks, process))
        if not expected or actual != expected:
            raise ValueError(f'{process} selection differs: expected={sorted(expected)}, actual={sorted(actual)}')
    count = check_gene_profiles(tasks, outdir)
    # Benchmark task count is configuration-dependent; never pin the whole graph's task total.
    counts = Counter(task_identity(task)[0] for task in tasks)
    print(f'Current workflow inventory passed ({mode}): {len(tasks)} successful tasks, '
          f'{count} selected/planted gene profiles, {counts["RUN_INTEGRATION_METHOD"]} benchmark methods.')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--mode', choices=('full', 'analysis'), default='full')
    parser.add_argument('--cell-types', required=True)
    parser.add_argument('--clock-cell-type', default='cd14-monocyte')
    args = parser.parse_args()
    check_inventory(args.outdir, mode=args.mode, cell_types=args.cell_types, clock_cell_type=args.clock_cell_type)


if __name__ == '__main__':
    main()

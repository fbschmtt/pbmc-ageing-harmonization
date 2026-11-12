"""Reject disconnected branches and missing/stale requests without scientific execution."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.check_workflow_inventory import (
    ANALYSIS_PROCESSES,
    CORE_PROCESSES,
    check_gene_profiles,
    require_branches,
)


def workflow_summary(mode='full'):
    processes = ANALYSIS_PROCESSES | (CORE_PROCESSES if mode == 'full' else set())
    return {'status': 'complete', 'exitcode': 0, 'tasks': [
        {'name': f'WORKFLOW:{process}', 'status': 'COMPLETED', 'exit': '0'} for process in sorted(processes)
    ]}


@pytest.mark.parametrize('mode', ['full', 'analysis'])
def test_expected_inventory_accepts_successful_or_cached_tasks(mode):
    summary = workflow_summary(mode)
    summary['tasks'][0]['status'] = 'CACHED'
    assert require_branches(summary, mode) == summary['tasks']


@pytest.mark.parametrize('process', ['HARMONIZE', 'SPLIT_CELL_TYPES', 'FIT_GENE_PROFILE',
                                     'RUN_INTEGRATION_METHOD', 'GENERATE_SYNTHETIC_PSEUDOBULK'])
def test_complete_run_cannot_hide_a_disconnected_branch(process):
    summary = workflow_summary()
    summary['tasks'] = [task for task in summary['tasks'] if task['name'] != f'WORKFLOW:{process}']
    with pytest.raises(ValueError, match=f'Missing current workflow processes:.*{process}'):
        require_branches(summary, 'full')


def test_analysis_mode_must_not_rebuild_core():
    with pytest.raises(ValueError, match='unexpectedly executed core'):
        require_branches(workflow_summary(), 'analysis')


@pytest.mark.parametrize(('status', 'exitcode'), [('FAILED', '1'), ('ABORTED', '-'), ('COMPLETED', '2')])
def test_complete_label_does_not_override_task_failure(status, exitcode):
    summary = workflow_summary()
    summary['tasks'][0].update(status=status, exit=exitcode)
    with pytest.raises(ValueError, match='Unsuccessful current task'):
        require_branches(summary, 'full')


@pytest.fixture
def gene_run(tmp_path):
    tasks = []
    outdir = tmp_path / 'published'

    def task(process, tag):
        workdir = tmp_path / 'work' / process / tag
        workdir.mkdir(parents=True)
        record = {'name': f'TOP_GENE_PROFILES:GENE_PROFILES:{process} ({tag})',
                  'status': 'COMPLETED', 'exit': '0', 'workdir': str(workdir)}
        tasks.append(record)
        return workdir

    selected = [('cd14-monocyte', 'CD14 monocyte', 'SYNTH_AGE_MARKER_2', 'synth-age-marker-2'),
                ('naive-cd4-t-cell', 'Naive CD4 T cell', 'SYNTH_AGE_MARKER_1', 'synth-age-marker-1')]
    for slug, cell_type, gene, _ in selected:
        workdir = task('SELECT_PROFILE_GENES', slug)
        (workdir / 'genes.tsv').write_text(f'{cell_type}\t{gene}\n')
    requested = [*selected, ('cd14-monocyte', 'CD14 monocyte', 'SYNTH_AGE_MARKER_1', 'synth-age-marker-1')]
    for slug in sorted({row[0] for row in requested}):
        task('PREPARE_GENE_PROFILE_MODEL', slug)
    for slug, cell_type, gene, gene_slug in requested:
        key = f'{slug}/{gene_slug}'
        published = outdir / 'analysis/per_gene' / key
        published.mkdir(parents=True)
        for process, files in [
            ('FIT_GENE_PROFILE', {'fit_summary.json': json.dumps({'cell_type': cell_type, 'gene': gene}),
                                 'gene_age_result.csv': f'gene,estimate\n{gene},0.1\n'}),
            ('RENDER_GENE_PROFILE', {'report.html': f'<html>{key}</html>', 'executed.ipynb': '{}',
                                    'report_render.log': f'Rendered {key}'}),
        ]:
            directory = task(process, key) / 'profiles' / key
            directory.mkdir(parents=True)
            for name, content in files.items():
                (directory / name).write_text(content)
                (published / name).write_text(content)
    return tasks, outdir


def test_all_selected_and_planted_genes_have_matching_fit_and_report_outputs(gene_run):
    tasks, outdir = gene_run
    assert check_gene_profiles(tasks, outdir) == 3


@pytest.mark.parametrize('process', ['FIT_GENE_PROFILE', 'RENDER_GENE_PROFILE'])
def test_stale_published_gene_cannot_replace_a_missing_current_task(gene_run, process):
    tasks, outdir = gene_run
    missing = f'{process} (naive-cd4-t-cell/synth-age-marker-1)'
    tasks[:] = [task for task in tasks if not task['name'].endswith(missing)]
    with pytest.raises(ValueError, match=f'{process} requests differ: missing='):
        check_gene_profiles(tasks, outdir)


def test_changed_selection_cannot_silently_drop_a_gene_request(gene_run):
    tasks, outdir = gene_run
    selector = next(task for task in tasks if 'SELECT_PROFILE_GENES (cd14-monocyte)' in task['name'])
    with (Path(selector['workdir']) / 'genes.tsv').open('a') as handle:
        handle.write('CD14 monocyte\tUNEXPECTED_GENE\n')
    with pytest.raises(ValueError, match='FIT_GENE_PROFILE requests differ'):
        check_gene_profiles(tasks, outdir)


def test_duplicate_model_preparation_is_rejected(gene_run):
    tasks, outdir = gene_run
    tasks.append(next(task for task in tasks if 'PREPARE_GENE_PROFILE_MODEL' in task['name']))
    with pytest.raises(ValueError, match='Duplicate PREPARE_GENE_PROFILE_MODEL'):
        check_gene_profiles(tasks, outdir)


def test_slug_collision_is_rejected_before_accepting_profile_outputs(gene_run):
    tasks, outdir = gene_run
    selector = next(task for task in tasks if 'SELECT_PROFILE_GENES (cd14-monocyte)' in task['name'])
    with (Path(selector['workdir']) / 'genes.tsv').open('a') as handle:
        handle.write('CD14 monocyte\tSYNTH-AGE-MARKER-2\n')
    with pytest.raises(ValueError, match='colliding selected gene'):
        check_gene_profiles(tasks, outdir)


def test_published_report_must_match_its_current_producer(gene_run):
    tasks, outdir = gene_run
    (outdir / 'analysis/per_gene/naive-cd4-t-cell/synth-age-marker-1/report.html').write_text('Old report')
    with pytest.raises(ValueError, match='differs from current task'):
        check_gene_profiles(tasks, outdir)


def test_old_published_report_cannot_mask_missing_current_output(gene_run):
    tasks, outdir = gene_run
    renderer = next(task for task in tasks if 'RENDER_GENE_PROFILE (naive-cd4-t-cell/' in task['name'])
    (Path(renderer['workdir']) / 'profiles/naive-cd4-t-cell/synth-age-marker-1/report.html').unlink()
    with pytest.raises(ValueError, match='Missing current/published artifact'):
        check_gene_profiles(tasks, outdir)


def test_profile_identity_must_match_selected_gene_even_when_files_agree(gene_run):
    tasks, outdir = gene_run
    fit = next(task for task in tasks if 'FIT_GENE_PROFILE (naive-cd4-t-cell/' in task['name'])
    key = 'naive-cd4-t-cell/synth-age-marker-1'
    for directory in (Path(fit['workdir']) / 'profiles' / key, outdir / 'analysis/per_gene' / key):
        (directory / 'fit_summary.json').write_text('{"gene":"WRONG_GENE","cell_type":"Naive CD4 T cell"}')
    with pytest.raises(ValueError, match='Profile identity differs from selection'):
        check_gene_profiles(tasks, outdir)


def test_planted_gene_cannot_hide_an_entirely_disconnected_ranking_step(gene_run):
    tasks, outdir = gene_run
    tasks[:] = [task for task in tasks if 'SELECT_PROFILE_GENES' not in task['name']]
    with pytest.raises(ValueError, match='No ranked gene profiles'):
        check_gene_profiles(tasks, outdir)

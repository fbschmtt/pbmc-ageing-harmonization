"""Public orchestration contracts: paths, failure reporting, and task ownership."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from scripts.run_workflow import run_workflow, summarize_trace

ROOT = Path(__file__).resolve().parents[1]


def test_trace_summary_keeps_failed_tasks_visible_without_changing_successful_exit(tmp_path):
    trace = tmp_path / 'trace.tsv'
    trace.write_text('name\tstatus\texit\nfit (good)\tCOMPLETED\t0\nfit (bad)\tFAILED\t17\n')
    summary = summarize_trace(trace, returncode=0, command=['nextflow', 'run'])
    assert summary['status'] == 'partial'
    assert summary['exitcode'] == 0
    assert summary['failed_tasks'] == [{'name': 'fit (bad)', 'status': 'FAILED', 'exit': '17'}]
    assert len(summary['tasks']) == 2


def test_workflow_finalizer_runs_once_even_after_essential_failure(tmp_path, monkeypatch):
    calls = []

    def launch(command, **kwargs):
        trace = Path(command[command.index('-with-trace') + 1])
        trace.write_text('name\tstatus\texit\ncore\tFAILED\t2\n')
        return subprocess.CompletedProcess(command, 2)

    monkeypatch.setattr('scripts.run_workflow.subprocess.run', launch)
    monkeypatch.setattr('scripts.run_workflow.render_report_index', lambda outdir: calls.append(outdir))
    assert run_workflow(['nextflow', 'run', 'main.nf'], tmp_path, 'core') == 2
    assert calls == [tmp_path]
    summary = json.loads((tmp_path / 'workflow_summary.json').read_text())
    assert summary['status'] == 'failed' and summary['exitcode'] == 2
    assert list((tmp_path / 'provenance/workflows').glob('core-*/summary.json'))


def test_index_failure_does_not_hide_the_original_workflow_failure(tmp_path, monkeypatch):
    monkeypatch.setattr('scripts.run_workflow.subprocess.run', lambda command, **kw: subprocess.CompletedProcess(command, 7))

    def broken_index(outdir):
        raise ValueError('bad inventory')

    monkeypatch.setattr('scripts.run_workflow.render_report_index', broken_index)
    assert run_workflow(['nextflow', 'run', 'main.nf'], tmp_path, 'core') == 7
    summary = json.loads((tmp_path / 'workflow_summary.json').read_text())
    assert summary['index_error'] == 'bad inventory'


@pytest.mark.parametrize('target', [
    'run-test', 'run-de-test', 'run-cell-type-analysis-test', 'run-integration-benchmark-test',
    'run-cell-type-test', 'render-cell-type-test', 'split-cell-types-test',
    'check-age-clock-synthetic-test', 'check-synthetic-de-results', 'check-integration-benchmark-test',
])
def test_test_targets_respect_the_requested_output_root(target):
    result = subprocess.run(
        ['make', '-n', target, 'OUTDIR=output/custom-contract-test', 'STUDIES=cima'],
        cwd=ROOT, text=True, capture_output=True, check=True,
    )
    assert 'output/custom-contract-test' in result.stdout
    assert 'output/test' not in result.stdout
    assert '-resume' not in result.stdout


def test_composed_suite_has_one_nextflow_invocation_and_no_direct_scientific_batches():
    result = subprocess.run(
        ['make', '-n', 'run-all-test', 'OUTDIR=output/custom-contract-test', 'STUDIES=cima'],
        cwd=ROOT, text=True, capture_output=True, check=True,
    )
    assert result.stdout.count('nextflow run pipeline.nf') == 1
    assert result.stdout.count('docker build -f docker/python.Dockerfile') == 1
    assert 'nextflow run main.nf' not in result.stdout
    assert 'python -m pbmc_pipeline.synthetic_de' not in result.stdout
    assert 'python scripts/profile_top_age_genes.py' not in result.stdout
    assert 'output/test' not in result.stdout
    assert '-resume' not in result.stdout


def test_partial_workflow_summary_is_not_treated_as_success_in_strict_suite(tmp_path):
    # The public strict validator must reject partial runs even when old artifacts exist.
    (tmp_path / 'workflow_summary.json').write_text(json.dumps({'status': 'partial', 'tasks': [{}]}))
    result = subprocess.run([
        str(ROOT / '.venv/bin/python'), 'scripts/check_workflow_suite.py', '--outdir', str(tmp_path),
        '--cell-types', 'cd14-monocyte', '--methods', 'elastic_net',
    ], cwd=ROOT, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert 'Integration workflow is incomplete' in result.stderr


def test_successful_task_with_partial_analysis_marks_the_run_partial(tmp_path):
    workdir = tmp_path / 'task'
    workdir.mkdir()
    (workdir / 'integration_benchmark.json').write_text('{"status":"partial","methods":[]}')
    trace = tmp_path / 'trace.tsv'
    trace.write_text('name\tstatus\texit\tworkdir\n'
                     f'INTEGRATION_BENCHMARK:AGGREGATE_INTEGRATION_METHODS\tCOMPLETED\t0\t{workdir}\n')
    summary = summarize_trace(trace, returncode=0, command=['nextflow', 'run'])
    assert summary['status'] == 'partial'
    assert summary['failed_tasks'] == []
    assert summary['partial_artifacts'][0]['status'] == 'partial'


@pytest.mark.parametrize('contents', [None, '{invalid', '[]'])
def test_unreadable_application_summary_cannot_make_a_run_complete(tmp_path, contents):
    workdir = tmp_path / 'task'
    workdir.mkdir()
    if contents is not None:
        (workdir / 'cell_type_manifest.json').write_text(contents)
    trace = tmp_path / 'trace.tsv'
    trace.write_text('name\tstatus\texit\tworkdir\n'
                     f'CELL_TYPE_ANALYSIS:WRITE_CELL_TYPE_MANIFEST\tCOMPLETED\t0\t{workdir}\n')
    summary = summarize_trace(trace, returncode=0, command=['nextflow', 'run'])
    assert summary['status'] == 'partial'
    assert summary['partial_artifacts'][0]['status'] == 'error'


def test_ineligible_clock_task_retains_cohort_provenance(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from pbmc_pipeline.age_clock import ClockEligibilityError
    from scripts.clock_task import fit

    captured = {}

    def ineligible(input_path, config_path, external, **kwargs):
        captured.update(external=external, methods=kwargs['method_names'])
        raise ClockEligibilityError('No eligible external samples')

    monkeypatch.setattr('scripts.clock_task.fit_age_clocks', ineligible)
    args = SimpleNamespace(
        input=tmp_path/'input.h5ad', config=tmp_path/'pipeline.json', cell_type='CD14 monocyte',
        slug='cd14-monocyte', input_identity='source/pseudobulk.h5ad', external_studies='a, b',
        alphas='0.1', l1_ratios='0.5', methods='elastic_net, rank_elastic_net',
        calibration_samples=10, calibration_seed=42, n_jobs=1, status=tmp_path/'status.json',
    )
    fit(args)
    record = json.loads(args.status.read_text())
    assert record['status'] == 'ineligible'
    assert record['input'] == 'source/pseudobulk.h5ad'
    assert record['external_studies'] == ['a', 'b']
    assert captured == {'external': ['a', 'b'], 'methods': ['elastic_net', 'rank_elastic_net']}


def test_clock_aggregate_keeps_provenance_when_every_type_is_ineligible(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from scripts.clock_task import aggregate

    status_dir = tmp_path/'statuses'
    status_dir.mkdir()
    (status_dir/'type.json').write_text(json.dumps({
        'cell_type': 'CD14 monocyte', 'status': 'ineligible',
        'input': 'source/pseudobulk.h5ad', 'external_studies': ['external'],
        'message': 'No eligible external samples',
    }))
    expected = tmp_path/'expected.json'
    expected.write_text('["CD14 monocyte"]')
    monkeypatch.setattr('scripts.clock_task.aggregate_cross_type_clocks', lambda *a, **kw: {'methods': []})
    args = SimpleNamespace(input_root=tmp_path/'types', status_dir=status_dir, expected=expected,
                           existing=False, template=None, output=tmp_path/'summary.json')
    aggregate(args)
    summary = json.loads(args.output.read_text())
    assert summary['status'] == 'complete'
    assert summary['input'] == 'source/pseudobulk.h5ad'
    assert summary['external_studies'] == ['external']
    assert summary['timestamp_utc']
    assert len(summary['eligibility_skipped_cell_types']) == 1
    assert summary['failed_cell_types'] == []

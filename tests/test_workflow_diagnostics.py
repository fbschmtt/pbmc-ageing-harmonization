from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.run_workflow import run_workflow, summarize_trace
from scripts.workflow_diagnostics import failure_details, format_failures

ROOT = Path(__file__).resolve().parents[1]


def failed_model_record(status='failed'):
    return {'cell_type': 'Naive CD4 T cell', 'models': [{
        'model': 'combined', 'study': None, 'covariates': ['cmv'],
        'design': '~ age + cmv', 'n_samples': 8,
        'status': status, 'error': 'LinAlgError: Singular matrix',
    }]}


def write_model_task(tmp_path, *, status='failed'):
    directory = tmp_path / 'fit'
    metadata = directory / 'cell_type_result/naive-cd4-t-cell/run_metadata.json'
    metadata.parent.mkdir(parents=True)
    metadata.write_text(json.dumps(failed_model_record(status)))
    return {'name': 'DE:RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION (Naive CD4 T cell)',
            'status': 'COMPLETED', 'exit': '0', 'workdir': str(directory)}


def test_workflow_end_prints_and_saves_task_error_without_changing_exit(tmp_path, monkeypatch, capsys):
    directory = tmp_path / 'task'
    directory.mkdir()
    error = directory / '.command.err'
    error.write_text('Traceback (most recent call last):\nRuntimeError: fitting failed\n')

    def launch(command, **kwargs):
        trace = Path(command[command.index('-with-trace') + 1])
        trace.write_text(f'name\tstatus\texit\tworkdir\nfit (bad)\tFAILED\t17\t{directory}\n')
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr('scripts.run_workflow.subprocess.run', launch)
    monkeypatch.setattr('scripts.run_workflow.render_report_index', lambda _: None)

    assert run_workflow(['nextflow', 'run', 'main.nf'], tmp_path, 'test') == 0

    summary = json.loads((tmp_path / 'workflow_summary.json').read_text())
    assert summary['status'] == 'partial'
    assert summary['task_diagnostics'][0]['log'] == str(error)
    assert summary['task_diagnostics'][0]['message'] == 'RuntimeError: fitting failed'
    output = capsys.readouterr().out
    assert 'fit (bad) [exit 17]' in output and 'RuntimeError: fitting failed' in output
    assert str(error) in output


@pytest.mark.parametrize('collector_status', ['COMPLETED', 'FAILED'])
def test_model_failure_is_visible_for_partial_or_rejected_collection(tmp_path, collector_status):
    producer = write_model_task(tmp_path)
    collector = tmp_path / 'collector'
    manifest = collector / 'analysis/combined/differential_expression/differential_expression.json'
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({'status': 'partial', 'cell_types': [failed_model_record()]}))
    trace = tmp_path / 'trace.tsv'
    trace.write_text('name\tstatus\texit\tworkdir\n'
                     f"{producer['name']}\tCOMPLETED\t0\t{producer['workdir']}\n"
                     f'DE:WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST\t{collector_status}\t'
                     f'{0 if collector_status == "COMPLETED" else 1}\t{collector}\n')

    summary = summarize_trace(trace, returncode=0, command=['nextflow', 'run'])

    assert summary['status'] == 'partial'
    assert len(summary['analysis_failures']) == 1  # Collector repeats the same fit.
    assert summary['analysis_failures'][0]['n_samples'] == 8
    output = format_failures(summary)
    assert 'Naive CD4 T cell / combined / cmv (8 samples)' in output
    assert 'LinAlgError: Singular matrix' in output
    assert len(summary['failed_tasks']) == (collector_status == 'FAILED')


def test_failed_model_marks_successful_standalone_fit_partial(tmp_path):
    task = write_model_task(tmp_path)
    trace = tmp_path / 'trace.tsv'
    trace.write_text(f"name\tstatus\texit\tworkdir\n{task['name']}\tCOMPLETED\t0\t{task['workdir']}\n")

    summary = summarize_trace(trace, returncode=0, command=['nextflow', 'run'])

    assert summary['status'] == 'partial'
    assert summary['exitcode'] == 0 and summary['failed_tasks'] == []
    assert len(summary['analysis_failures']) == 1


def test_skipped_models_are_not_reported_as_failures(tmp_path):
    task = write_model_task(tmp_path, status='skipped')

    assert failure_details([task]) == {'task_diagnostics': [], 'analysis_failures': []}


@pytest.mark.parametrize('contents', ['{invalid', '[]', '{"models":null}'])
def test_bad_diagnostic_metadata_does_not_hide_task_failures(tmp_path, contents):
    task = write_model_task(tmp_path)
    metadata = next(Path(task['workdir']).glob('cell_type_result/*/run_metadata.json'))
    metadata.write_text(contents)
    task.update(status='FAILED', exit='1')

    details = failure_details([task])

    assert len(details['task_diagnostics']) == 1
    assert 'No error text available' in details['task_diagnostics'][0]['message']


def test_saved_excerpt_is_bounded_and_falls_back_to_combined_log(tmp_path):
    (tmp_path / '.command.err').write_text('')
    combined = tmp_path / '.command.log'
    combined.write_text('irrelevant progress\n' * 10000 + 'Tool crashed: out of memory\n')
    details = failure_details([{'name': 'fit', 'status': 'FAILED', 'exit': '137', 'workdir': str(tmp_path)}])

    diagnostic = details['task_diagnostics'][0]
    assert diagnostic['log'] == str(combined)
    assert diagnostic['message'] == 'Tool crashed: out of memory'
    assert len(diagnostic['log_excerpt']) <= 4096


def test_collector_metadata_reports_model_errors_without_producer_tasks(tmp_path):
    manifest = tmp_path / 'analysis/combined/differential_expression/differential_expression.json'
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({'cell_types': [failed_model_record()]}))

    details = failure_details([{'name': 'WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST',
                               'status': 'FAILED', 'exit': '1', 'workdir': str(tmp_path)}])

    assert details['analysis_failures'][0]['message'] == 'LinAlgError: Singular matrix'


def test_read_only_cli_supports_old_summaries_without_changing_them(tmp_path):
    task = write_model_task(tmp_path)
    path = tmp_path / 'summary.json'
    original = json.dumps({'status': 'partial', 'tasks': [task]})
    path.write_text(original)

    result = subprocess.run([sys.executable, '-m', 'scripts.show_workflow_failures', '--summary', str(path)],
                            cwd=ROOT, capture_output=True, text=True, check=True)

    assert 'Failed analysis models (1)' in result.stdout
    assert 'LinAlgError: Singular matrix' in result.stdout
    assert path.read_text() == original


def test_saved_diagnostics_remain_visible_after_work_directory_is_removed():
    summary = {'task_diagnostics': [{'name': 'fit', 'exit': '1', 'message': 'Recorded failure'}],
               'analysis_failures': []}

    assert 'Recorded failure' in format_failures(summary)

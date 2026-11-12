from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.run_workflow import run_workflow
from scripts.workflow_test_output import expected_failure_test


def test_expected_failure_output_is_labeled_and_pass_waits_for_assertions(tmp_path, capsys):
    log = tmp_path / 'console.log'

    with expected_failure_test('main/de-report', 'inject report exit 43', log):
        print('Failed tasks (1): deliberate report failure')
        print('Task failure stderr', file=sys.stderr)
        assert '[PASS]' not in capsys.readouterr().out

    output = capsys.readouterr().out
    assert '[PASS] main/de-report: expected failures and required outcomes verified.' in output
    assert 'Failed tasks (1)' not in output
    assert 'Failed tasks (1)' in log.read_text()
    assert 'Task failure stderr' in log.read_text()


def test_expected_failure_banner_precedes_captured_errors(tmp_path, capsys):
    log = tmp_path / 'console.log'

    with expected_failure_test('main/de-report', 'inject report exit 43', log):
        banner = capsys.readouterr().out
        assert '[EXPECTED FAILURE TEST] main/de-report: inject report exit 43' in banner
        assert str(log) in banner


def test_unexpected_test_result_shows_diagnostics_and_still_raises(tmp_path, capsys):
    log = tmp_path / 'console.log'

    with (
        pytest.raises(AssertionError, match='unexpected task failure'),
        expected_failure_test('main/de-report', 'inject report exit 43', log),
    ):
        print('RuntimeError: an unrelated task crashed')
        raise AssertionError('unexpected task failure')

    output = capsys.readouterr().out
    assert '[FAIL] main/de-report: unexpected result' in output
    assert 'RuntimeError: an unrelated task crashed' in output
    assert str(log) in output
    assert '[PASS]' not in output
    print('stdout restored')
    assert 'stdout restored' in capsys.readouterr().out


def test_captured_workflow_retains_real_failure_status_and_child_output(tmp_path, monkeypatch, capsys):
    run = subprocess.run

    def launch(command, **kwargs):
        trace = Path(command[command.index('-with-trace') + 1])
        trace.write_text('name\tstatus\texit\nreport (failed-type)\tFAILED\t43\n')
        return run([sys.executable, '-c',
                    'import sys; print("child stdout"); print("child stderr", file=sys.stderr); sys.exit(43)'],
                   **kwargs)

    monkeypatch.setattr('scripts.run_workflow.subprocess.run', launch)
    monkeypatch.setattr('scripts.run_workflow.render_report_index', lambda _: None)
    log = tmp_path / 'console.log'

    with expected_failure_test('fatal-test', 'inject task exit 43', log) as output:
        assert run_workflow(['nextflow', 'run', 'main.nf'], tmp_path, 'fatal-test', output=output) == 43
        summary = json.loads((tmp_path / 'workflow_summary.json').read_text())
        assert summary['status'] == 'failed' and summary['exitcode'] == 43
        assert summary['failed_tasks'][0]['exit'] == '43'

    contents = log.read_text()
    assert 'child stdout' in contents and 'child stderr' in contents
    assert 'Failed tasks (1)' in contents
    assert 'Workflow failed:' not in capsys.readouterr().out

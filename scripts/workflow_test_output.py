"""Clearly separate deliberate workflow failures from failed verification checks."""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path

from scripts.workflow_diagnostics import _log_tail


@contextmanager
def expected_failure_test(name: str, description: str, log_path: Path):
    """Show a pass only after every assertion in the test case succeeds."""
    print(f'[EXPECTED FAILURE TEST] {name}: {description}', flush=True)
    print(f'  Failure is deliberately injected; full console output: {log_path}', flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open('w') as output, redirect_stdout(output), redirect_stderr(output):
            yield output
    except Exception:
        print(f'[FAIL] {name}: unexpected result in failure-isolation test.', flush=True)
        print(f'  Full diagnostics: {log_path}', flush=True)
        tail = _log_tail(log_path, maximum_bytes=8192)
        if tail:
            print(tail, flush=True)
        raise
    print(f'[PASS] {name}: expected failures and required outcomes verified.', flush=True)

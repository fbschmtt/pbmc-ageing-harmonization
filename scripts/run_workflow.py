"""Run a workflow, record partial failures, and refresh the workspace index once."""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO
from uuid import uuid4

from scripts.render_report_index import render_report_index
from scripts.workflow_diagnostics import failure_details, format_failures

# These are declared output contracts, resolved only inside current task work directories.
APPLICATION_SUMMARIES = {
    'WRITE_CELL_TYPE_MANIFEST': 'cell_type_manifest.json',
    'WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST': 'analysis/combined/differential_expression/differential_expression.json',
    'AGGREGATE_INTEGRATION_METHODS': 'integration_benchmark.json',
    'AGGREGATE_AGE_CLOCKS': 'age_clock_all_types_summary.json',
}


def summarize_trace(trace: Path, *, returncode: int, command: list[str]) -> dict:
    records = []
    if trace.is_file():
        with trace.open() as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
    failed = [row for row in records if row.get('status') in {'FAILED', 'ABORTED'}]
    partial_artifacts = []
    for row in records:
        process = row.get('name', '').split(':')[-1].split(' (')[0]
        if process not in APPLICATION_SUMMARIES or not row.get('workdir'):
            continue
        if row.get('status') not in {'COMPLETED', 'CACHED'}:
            continue
        path = Path(row['workdir']) / APPLICATION_SUMMARIES[process]
        try:
            document = json.loads(path.read_text())
            if not isinstance(document, dict):
                raise TypeError('Application summary must be an object')
        except (OSError, ValueError, TypeError) as error:
            partial_artifacts.append({'process': row['name'], 'path': str(path), 'status': 'error', 'reason': str(error)})
            continue
        if document.get('status') in {'partial', 'failed', 'error'}:
            partial_artifacts.append({'process': row['name'], 'path': str(path), 'status': document['status']})
    details = failure_details(records)
    return {
        'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'failed' if returncode else ('partial' if failed or partial_artifacts or details['analysis_failures'] else 'complete'),
        'exitcode': returncode, 'command': command, 'tasks': records, 'failed_tasks': failed,
        'partial_artifacts': partial_artifacts,
        **details,
    }


def run_workflow(command: list[str], outdir: Path, name: str, *, output: TextIO | None = None) -> int:
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid4().hex[:8]
    record_dir = outdir.resolve() / 'provenance/workflows' / f'{name}-{run_id}'
    record_dir.mkdir(parents=True)
    trace = record_dir / 'trace.tsv'
    # An explicit trace belongs to this invocation even when other runs share a work directory.
    run_position = command.index('run')
    log = record_dir / 'nextflow.log'
    invocation = [*command[:run_position], '-log', str(log), *command[run_position:], '-with-trace', str(trace)]
    result = subprocess.run(invocation, check=False, stdout=output,
                            stderr=subprocess.STDOUT if output is not None else None)
    summary = summarize_trace(trace, returncode=result.returncode, command=invocation)
    summary.update(name=name, run_id=run_id)
    if log.is_file():
        session = re.search(r'Session UUID: ([a-f0-9-]+)', log.read_text())
        if session:
            summary['session_id'] = session.group(1)
    finalization_error = None
    try:
        render_report_index(outdir)
    except Exception as error:  # noqa: BLE001 - preserve the workflow exit status and record finalization failure
        finalization_error = str(error)
        summary['index_error'] = finalization_error
        summary['status'] = 'failed'
    document = json.dumps(summary, indent=2) + '\n'
    (record_dir / 'summary.json').write_text(document)
    latest = outdir / f'.workflow_summary.{run_id}.tmp'
    latest.write_text(document)
    latest.replace(outdir / 'workflow_summary.json')
    print(f"Workflow {summary['status']}: {record_dir / 'summary.json'}", file=output, flush=True)
    if summary['status'] != 'complete':
        print(format_failures(summary), file=output, flush=True)
    if finalization_error:
        print(f'Index finalization failed: {finalization_error}', file=output if output is not None else sys.stderr)
    return result.returncode or (1 if finalization_error else 0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('Specify the Nextflow command after --')
    raise SystemExit(run_workflow(command, args.outdir, args.name))


if __name__ == '__main__':
    main()

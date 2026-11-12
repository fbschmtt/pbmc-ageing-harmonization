"""Describe task and model failures from explicit current-run artifacts."""
from __future__ import annotations

import json
from pathlib import Path


def _log_tail(path: Path, maximum_bytes: int = 4096) -> str:
    try:
        with path.open('rb') as handle:
            handle.seek(0, 2)
            size = handle.tell()
            handle.seek(max(0, size - maximum_bytes))
            return handle.read(maximum_bytes).decode('utf-8', errors='replace').strip()
    except OSError:
        return ''


def failure_details(tasks: list[dict]) -> dict:
    """Keep diagnostics additive; trace rows and application exit codes stay intact."""
    task_diagnostics = []
    analysis_failures = {}
    for task in tasks:
        directory = Path(task['workdir']) if task.get('workdir') else None
        if task.get('status') in {'FAILED', 'ABORTED'}:
            detail = {key: task[key] for key in ('name', 'hash', 'status', 'exit', 'workdir') if key in task}
            excerpt = ''
            if directory is not None:
                for filename in ('.command.err', '.command.log', '.command.out'):
                    path = directory / filename
                    excerpt = _log_tail(path)
                    if excerpt:
                        detail.update(log=str(path), log_excerpt=excerpt)
                        break
            detail['message'] = (
                excerpt.splitlines()[-1] if excerpt else
                f"No error text available; task status {task.get('status')}, exit {task.get('exit', 'unknown')}"
            )
            task_diagnostics.append(detail)
        if directory is None:
            continue
        process = task.get('name', '').split(':')[-1].split(' (')[0]
        if process == 'RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION':
            paths = sorted(directory.glob('cell_type_result/*/run_metadata.json'))
        elif process == 'WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST':
            # The collector writes this contract before rejecting partial fits.
            paths = [directory / 'analysis/combined/differential_expression/differential_expression.json']
        else:
            continue
        for path in paths:
            try:
                document = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            if not isinstance(document, dict):
                continue
            records = document.get('cell_types', []) if process.startswith('WRITE_') else [document]
            if not isinstance(records, list):
                continue
            for record in records:
                if not isinstance(record, dict) or not isinstance(record.get('models'), list):
                    continue
                for model in record['models']:
                    if not isinstance(model, dict) or model.get('status') != 'failed':
                        continue
                    detail = {
                        'process': task['name'], 'path': str(path),
                        'cell_type': record.get('cell_type'),
                        **{key: model.get(key) for key in ('model', 'study', 'covariates', 'design', 'n_samples')},
                        'message': model.get('error') or model.get('reason') or 'Model failed without a recorded error',
                    }
                    # A collector repeats the per-type records. Report each fit once.
                    key = json.dumps({k: v for k, v in detail.items() if k not in {'process', 'path'}}, sort_keys=True)
                    analysis_failures.setdefault(key, detail)
    return {'task_diagnostics': task_diagnostics, 'analysis_failures': list(analysis_failures.values())}


def format_failures(summary: dict, *, limit: int = 10) -> str:
    """Format saved diagnostics, also supporting summaries from older runs."""
    details = summary if 'task_diagnostics' in summary and 'analysis_failures' in summary else {
        **summary, **failure_details(summary.get('tasks', [])),
    }
    lines = []
    tasks = details['task_diagnostics']
    if tasks:
        lines.append(f'Failed tasks ({len(tasks)}):')
        for task in tasks[:limit]:
            lines.append(f"  {task['name']} [exit {task.get('exit', 'unknown')}]")
            lines.append(f"    {task['message']}")
            if task.get('log'):
                lines.append(f"    Log: {task['log']}")
            if task.get('workdir'):
                lines.append(f"    Work directory: {task['workdir']}")
        if len(tasks) > limit:
            lines.append(f'  {len(tasks) - limit} more failed tasks; see the saved workflow summary.')
    models = details['analysis_failures']
    if models:
        lines.append(f'Failed analysis models ({len(models)}):')
        for model in models[:limit]:
            label = ' / '.join(str(value) for value in (
                model['cell_type'], model['model'], model['study'], ', '.join(model['covariates'] or []),
            ) if value)
            lines.append(f"  {label} ({model['n_samples']} samples)")
            lines.append(f"    {model['message']}")
            lines.append(f"    Metadata: {model['path']}")
        if len(models) > limit:
            lines.append(f'  {len(models) - limit} more failed models; see the saved workflow summary.')
    for artifact in summary.get('partial_artifacts', []):
        lines.append(f"Analysis {artifact['status']}: {artifact['process']}")
        if artifact.get('reason'):
            lines.append(f"  {artifact['reason']}")
        lines.append(f"  Summary: {artifact['path']}")
    return '\n'.join(lines) or 'No recorded task or analysis failures.'

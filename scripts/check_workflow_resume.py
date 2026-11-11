"""Check same-run caching and report-only invalidation using the real profile graph."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import tempfile
from pathlib import Path

from scripts.run_workflow import run_workflow

ROOT = Path(__file__).resolve().parents[1]


def check_resume(nextflow: str) -> None:
    (ROOT / 'output/tests').mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='workflow-resume.', dir=ROOT / 'output/tests'))
    template = directory / 'template.py'
    template.write_text('# Initial report template\n')
    outdir = directory / 'published'
    command = [*shlex.split(nextflow), 'run', 'tests/workflows/gene_resume.nf', '-c', 'nextflow.config',
               '-profile', 'standard', '-stub-run', '-ansi-log', 'false', '-work-dir', str(directory / 'work'),
               '--outdir', str(outdir), '--project_dir', str(ROOT), '--resume_template', str(template)]
    os.environ['NXF_OFFLINE'] = 'true'

    def invoke(name, arguments):
        assert run_workflow(arguments, outdir, name) == 0
        result = json.loads((outdir / 'workflow_summary.json').read_text())
        assert result['status'] == 'complete', result
        assert len(result['tasks']) == 5, result
        return result

    first = invoke('initial', command)
    assert sum('PREPARE_GENE_PROFILE_MODEL' in task['name'] for task in first['tasks']) == 1
    resumed_command = [*command, '-resume', first['session_id']]
    second = invoke('resumed', resumed_command)
    assert {task['status'] for task in second['tasks']} == {'CACHED'}, second
    template.write_text('# Changed report template; scientific inputs remain identical\n')
    third = invoke('report-change', resumed_command)
    for task in third['tasks']:
        expected = 'COMPLETED' if 'RENDER_GENE_PROFILE' in task['name'] else 'CACHED'
        assert task['status'] == expected, task
    print(f'Resume checks passed: one shared fit, full reuse, report-only invalidation. Logs: {directory}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check_resume(parser.parse_args().nextflow)


if __name__ == '__main__':
    main()

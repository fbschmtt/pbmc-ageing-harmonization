"""Exercise the real validator and fatal fan-out barrier without scientific execution."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import tempfile
from contextlib import nullcontext
from pathlib import Path

from scripts.check_workflow_inventory import task_identity
from scripts.check_workflow_modes import ROOT, prepare_project
from scripts.run_workflow import run_workflow
from scripts.workflow_test_output import expected_failure_test


def check_preflight(nextflow: str) -> None:
    (ROOT / 'output/tests').mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='configuration-preflight.', dir=ROOT / 'output/tests'))
    project = directory / 'project'
    prepare_project(project)
    config = project / 'config/pipeline.json'
    original = config.read_bytes()
    os.environ['NXF_OFFLINE'] = 'true'
    for case in ('valid', 'invalid-analysis', 'invalid-clock', 'invalid-core-registry'):
        config.write_bytes(original)
        outdir = directory / case / 'published'
        if case == 'valid':
            entry = 'tests/workflows/configuration_preflight.nf'
            arguments = ['--preflight_config', str(config)]
        else:
            entry = 'main.nf'
            arguments = ['--mode', 'analysis', '--project_dir', str(project),
                         '--merged_input', str(project / 'merged.h5ad'), '--pseudobulk_input', str(project / 'pseudobulk.h5ad')]
            if case == 'invalid-analysis':
                document = json.loads(original)
                document['differential_expression']['alpha'] = 0
                config.write_text(json.dumps(document))
            elif case == 'invalid-core-registry':
                study_file = project / 'config/studies.json'
                studies = json.loads(study_file.read_text())
                studies['studies']['fixture']['conversion'] = {'source': 'input.rds'}
                study_file.write_text(json.dumps(studies))
                (project / 'input.rds').write_text('No conversion may start before validation')
                (project / 'scripts/convert_rds.R').write_text('# Must never execute in this test')
                arguments = ['--mode', 'full', '--studies', 'fixture', '--project_dir', str(project)]
            else:
                arguments += ['--run_clocks', '--clock_n_jobs', '0']
        command = [*shlex.split(nextflow), 'run', entry, '-profile', 'standard', '-ansi-log', 'false',
                   '-work-dir', str(directory / case / 'work'), '--outdir', str(outdir), *arguments]
        context = nullcontext(None) if case == 'valid' else expected_failure_test(
            f'configuration/{case}', 'reject invalid configuration before scientific tasks',
            directory / case / 'console.log',
        )
        with context as output:
            result = run_workflow(command, outdir, case, output=output)
            summary = json.loads((outdir / 'workflow_summary.json').read_text())
            if case == 'valid':
                assert result == 0 and summary['status'] == 'complete', summary
                assert (outdir / 'accepted.json').read_bytes() == original
                assert len(summary['tasks']) == 2
            else:
                assert result != 0 and summary['status'] == 'failed', summary
                assert len(summary['tasks']) == 1, summary['tasks']
                assert task_identity(summary['tasks'][0])[0] == 'VALIDATE_CONFIGURATION'
                workdir = Path(summary['tasks'][0]['workdir'])
                assert (workdir / '.exitcode').read_text().strip() == '2'
                error = (workdir / '.command.err').read_text()
                assert 'Configuration preflight failed' in error
                if case != 'invalid-core-registry':
                    assert ('differential_expression.alpha' if case == 'invalid-analysis' else 'clock.n_jobs') in error
            print(f'Real configuration preflight passed: {case}', flush=True)
    print(f'Configuration preflight logs: {directory}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check_preflight(parser.parse_args().nextflow)


if __name__ == '__main__':
    main()

"""Run the analysis failure-isolation cases with explicit expected-failure output."""
from __future__ import annotations

import argparse
import os
import shlex
import tempfile
from pathlib import Path

from scripts.check_analysis_failures import check_failures
from scripts.run_workflow import run_workflow
from scripts.workflow_test_output import expected_failure_test

ROOT = Path(__file__).resolve().parents[1]


def check_workflow_failures(nextflow: str) -> None:
    parent = ROOT / 'output/tests'
    parent.mkdir(parents=True, exist_ok=True)
    os.environ['NXF_OFFLINE'] = 'true'
    for phase, arguments, description in [
        ('fits', ['--stub_fail_gene', '--stub_fail_clock', '--stub_fail_cell_type', '--stub_fail_de'],
         'inject gene, clock, cell-type and DE fitting failures'),
        ('reports', ['--stub_fail_gene_report', '--stub_fail_clock_report', '--stub_fail_cell_report', '--stub_fail_de_report'],
         'inject gene, clock, cell-type and DE rendering failures'),
        ('all-cells', ['--stub_fail_all_cell_types'], 'inject failures for every cell-type analysis'),
    ]:
        directory = Path(tempfile.mkdtemp(prefix='analysis-failures.', dir=parent))
        outdir = directory / 'published'
        command = [*shlex.split(nextflow), 'run', 'tests/workflows/analysis_failures.nf',
                   '-c', 'nextflow.config', '-profile', 'standard', '-stub-run', '-ansi-log', 'false',
                   '-work-dir', str(directory / 'work'), '--outdir', str(outdir),
                   '--project_dir', str(ROOT), *arguments]
        with expected_failure_test(f'analysis/{phase}', description, directory / 'console.log') as output:
            result = run_workflow(command, outdir, 'analysis-failures', output=output)
            assert result == 0, f'Expected independent failures to leave Nextflow exit zero; got {result}'
            check_failures(outdir, phase)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check_workflow_failures(parser.parse_args().nextflow)


if __name__ == '__main__':
    main()

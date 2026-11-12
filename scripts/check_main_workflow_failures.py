"""Exercise single-task failure isolation through the production main entry point."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import tempfile
from pathlib import Path

from scripts.check_workflow_inventory import check_gene_profiles, task_identity, task_index
from scripts.check_workflow_modes import ROOT, prepare_project
from scripts.check_workflow_status import require_complete
from scripts.run_workflow import run_workflow
from scripts.workflow_test_output import expected_failure_test


def check_main_failures(nextflow: str) -> None:
    (ROOT / 'output/tests').mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='main-failures.', dir=ROOT / 'output/tests'))
    project = directory / 'project'
    prepare_project(project)
    # The manifest stub uses the real stdlib-only writer, never a fabricated success marker.
    (project / 'src').symlink_to(ROOT / 'src', target_is_directory=True)
    settings = json.loads((project / 'config/pipeline.json').read_text())
    methods = [method['name'] for method in settings['integration_benchmark']['methods']]
    if len(methods) < 2:
        raise ValueError('Failure isolation needs at least two configured benchmark methods')
    os.environ['NXF_OFFLINE'] = 'true'
    for name, arguments, failed_process, exitcode in [
        ('de-fit', ['--stub_fail_de'], 'RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION', '47'),
        ('de-report', ['--stub_fail_de_report'], 'RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS', '43'),
        ('benchmark', ['--stub_fail_benchmark_method', '*'], 'RUN_INTEGRATION_METHOD', '53'),
        ('preflight', ['--stub_fail_preflight'], 'VALIDATE_CONFIGURATION', '61'),
    ]:
        case = directory / name
        outdir = case / 'published'
        command = [*shlex.split(nextflow), 'run', 'main.nf',
                   '-profile', 'standard,cell_type_analysis,differential_expression,integration_benchmark,analysis',
                   '-stub-run', '-ansi-log', 'false', '-work-dir', str(case / 'work'),
                   '--mode', 'analysis', '--project_dir', str(project), '--outdir', str(outdir),
                   '--merged_input', str(project / 'merged.h5ad'),
                   '--pseudobulk_input', str(project / 'pseudobulk.h5ad'),
                   '--run_clocks', '--run_benchmark', *arguments]
        description = f'inject {failed_process} exit {exitcode}'
        with expected_failure_test(f'main/{name}', description, case / 'console.log') as output:
            result = run_workflow(command, outdir, name, output=output)
            summary = json.loads((outdir / 'workflow_summary.json').read_text())
            if name == 'preflight':
                assert result != 0 and summary['status'] == 'failed', summary
                assert len(summary['tasks']) == 1, summary['tasks']
                task = summary['tasks'][0]
                assert task_identity(task)[0] == failed_process and task['exit'] == exitcode, task
                assert (Path(task['workdir']) / '.exitcode').read_text().strip() == exitcode
                print('Main preflight failure passed: fatal exit, no scientific tasks dispatched.', flush=True)
                continue
            assert result == 0 and summary['status'] == 'partial', summary
            assert len(summary['failed_tasks']) == 1, summary['failed_tasks']
            failed = summary['failed_tasks'][0]
            assert task_identity(failed)[0] == failed_process and failed['exit'] == exitcode, failed
            for task in summary['tasks']:
                actual = (Path(task['workdir']) / '.exitcode').read_text().strip()
                if task['status'] == 'FAILED':
                    assert actual == exitcode, task
                else:
                    assert task['status'] == 'COMPLETED' and actual == '0', task
            assert not any(task_identity(task)[0] == 'HARMONIZE' for task in summary['tasks'])
            assert (outdir / 'index.html').is_file()
            manifest = json.loads((outdir / 'analysis/per_type/cell_type_manifest.json').read_text())
            assert manifest['status'] == 'complete', manifest
            assert {row['slug'] for row in manifest['cell_types']} == {'test-type', 'failed-type'}
            for slug in ('test-type', 'failed-type'):
                root = outdir / 'analysis/per_type' / slug
                assert (root / 'analysis/report.html').is_file()
                assert (root / 'ageing_clock/age_clock.json').is_file()
                assert (root / 'ageing_clock/report.html').is_file()
                received = root / 'analysis/received_merged_de.csv'
                if name == 'de-fit' and slug == 'failed-type':
                    assert not received.exists()
                    assert not (root / 'differential_expression/merged.csv').exists()
                else:
                    assert received.read_bytes() == (root / 'differential_expression/merged.csv').read_bytes()
                expected_report = not (slug == 'failed-type' and name in {'de-fit', 'de-report'})
                assert (root / 'differential_expression/report.html').exists() == expected_report
            profiles = check_gene_profiles(summary['tasks'], outdir, planted_gene=False)
            assert profiles == (1 if name == 'de-fit' else 2)
            benchmark = json.loads((outdir / 'harmonization/integration_benchmark/integration_benchmark.json').read_text())
            failed_method = task_identity(failed)[1].removeprefix('integration method ')
            expected_methods = set(methods) - ({failed_method} if name == 'benchmark' else set())
            assert {row['name'] for row in benchmark['methods']} == expected_methods, benchmark
            assert benchmark['status'] == ('partial' if name == 'benchmark' else 'complete'), benchmark
            assert (outdir / 'harmonization/integration_benchmark/report.html').is_file()
            if name == 'benchmark':
                # The production one-fork limit and first-task failure exercise queued siblings.
                siblings = list(task_index(summary['tasks'], failed_process).values())
                assert siblings[0]['status'] == 'FAILED', siblings
                assert all(task['status'] == 'COMPLETED' for task in siblings[1:]) and len(siblings) > 1
            try:
                require_complete(outdir)
            except SystemExit:
                pass
            else:
                raise AssertionError('Strict verification accepted a partial main workflow')
            print(f'Main failure isolation passed ({name}): one failure, continuing siblings and consumers, partial run.', flush=True)
    print(f'Main workflow failure logs: {directory}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check_main_failures(parser.parse_args().nextflow)


if __name__ == '__main__':
    main()

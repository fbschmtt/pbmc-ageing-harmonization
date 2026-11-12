"""Check the public main entry point's modes without Docker or scientific data."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def prepare_project(project: Path) -> None:
    settings = json.loads((ROOT / 'config/pipeline.json').read_text())
    # Inert inputs for previews and process stubs; never use them for scientific execution.
    paths = [settings['obs_schema'], settings['input_sources'], *settings['models'].values(),
             'config/studies.schema.json', 'config/input_sources.schema.json',
             'input.h5ad', 'merged.h5ad', 'pseudobulk.h5ad']
    for folder in ('scripts', 'reports'):
        paths.extend(str(path.relative_to(ROOT)) for path in (ROOT / folder).glob('*.py'))
    for name in paths:
        path = project / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('Graph preview placeholder; not scientific data.\n')
    (project / 'input_data/go').mkdir(parents=True)
    settings['merge']['single_cell']['enabled'] = False
    (project / 'config/pipeline.json').write_text(json.dumps(settings))
    (project / settings['studies_config']).write_text(json.dumps({
        'studies': {'fixture': {'input': 'input.h5ad', 'annotation': {'method': 'existing'}}},
    }))


def check_modes(nextflow: str) -> None:
    (ROOT / 'output/tests').mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='workflow-modes.', dir=ROOT / 'output/tests'))
    project = directory / 'project'
    prepare_project(project)
    inputs = ['--merged_input', str(project / 'merged.h5ad'),
              '--pseudobulk_input', str(project / 'pseudobulk.h5ad')]
    cases = [
        ('default', [], {'CORE', 'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES', 'CELL_TYPE_ANALYSIS', 'INTEGRATION_BENCHMARK'}, True),
        ('full', ['--mode', 'full', '--run_clocks', 'false'], {'CORE', 'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES', 'CELL_TYPE_ANALYSIS', 'INTEGRATION_BENCHMARK'}, True),
        ('full-options', ['--mode', 'full', '--run_benchmark', 'false', '--run_clocks'],
         {'CORE', 'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES', 'CELL_TYPE_ANALYSIS', 'AGE_CLOCKS'}, True),
        ('harmonization', ['--mode', 'harmonization', '--run_benchmark', '--run_clocks', '--merge_single_cell', 'false'], {'CORE'}, False),
        ('harmonization-merge', ['--mode', 'harmonization', '--merge_single_cell'], {'CORE'}, True),
        ('analysis', ['--mode', 'analysis', *inputs],
         {'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES', 'CELL_TYPE_ANALYSIS'}, False),
        ('analysis-benchmark', ['--mode', 'analysis', *inputs, '--run_benchmark'],
         {'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES', 'CELL_TYPE_ANALYSIS', 'INTEGRATION_BENCHMARK'}, False),
    ]
    environment = {**os.environ, 'NXF_OFFLINE': 'true'}

    def invoke(name, arguments):
        case = directory / name
        case.mkdir()
        command = [*shlex.split(nextflow), '-log', str(case / 'nextflow.log'), 'run', 'main.nf',
                   '-profile', 'standard', '-preview', '-ansi-log', 'false',
                   '-with-dag', str(case / 'graph.dot'), '-work-dir', str(case / 'work'),
                   '--outdir', str(case / 'published'), '--project_dir', str(project),
                   '--studies', 'fixture', *arguments]
        result = subprocess.run(command, cwd=ROOT, env=environment, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        (case / 'console.log').write_text(result.stdout)
        if list(case.rglob('.command.sh')) or (case / 'published').exists():
            raise RuntimeError(f'{name}: preview dispatched tasks or published outputs')
        return case, result

    for name, arguments, expected, merged in cases:
        case, result = invoke(name, arguments)
        if result.returncode:
            raise RuntimeError(f'{name}: graph preview failed\n{result.stdout}')
        graph = (case / 'graph.dot').read_text()
        workflows = {'CORE', 'DIFFERENTIAL_EXPRESSION', 'TOP_GENE_PROFILES',
                     'CELL_TYPE_ANALYSIS', 'AGE_CLOCKS', 'INTEGRATION_BENCHMARK'}
        actual = {name for name in workflows if f'label="{name}:' in graph}
        if actual != expected or ('CORE:MERGE_SINGLE_CELLS' in graph) != merged:
            raise RuntimeError(f'{name}: unexpected branches {actual}; see {case / "graph.dot"}')
        print(f'Passed workflow mode: {name}', flush=True)

    for name, arguments, message in [
        ('unknown', ['--mode', 'typo'], 'Unknown --mode'),
        ('zero-profiles', ['--top_age_gene_count', '0'], 'must be a positive integer'),
        ('invalid-profiles', ['--top_age_gene_count', 'typo'], 'must be a positive integer'),
        ('invalid-benchmark', ['--run_benchmark', 'maybe'], '--run_benchmark must be true or false'),
        ('missing-inputs', ['--mode', 'analysis'], 'requires both --merged_input and --pseudobulk_input'),
        ('missing-pseudobulk', ['--mode', 'analysis', *inputs[:2]], 'requires both --merged_input and --pseudobulk_input'),
        ('missing-file', ['--mode', 'analysis', *inputs, '--merged_input', str(project / 'absent.h5ad')], 'No such file or directory'),
    ]:
        _, result = invoke(name, arguments)
        if result.returncode == 0 or message not in result.stdout:
            raise RuntimeError(f'{name}: expected preflight rejection {message!r}\n{result.stdout}')
        print(f'Passed mode preflight rejection: {name}', flush=True)
    print(f'Workflow mode logs and graphs: {directory}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check_modes(parser.parse_args().nextflow)


if __name__ == '__main__':
    main()

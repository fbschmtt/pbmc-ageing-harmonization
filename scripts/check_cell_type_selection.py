"""Exercise production cell-type preflight without opening scientific inputs."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def check_selections(nextflow: str) -> None:
    (ROOT / "output/tests").mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="cell-type-selection.", dir=ROOT / "output/tests"))
    for slug in ("cd14-monocyte", "naive-cd4-t-cell"):
        split = directory / "inputs" / slug / "split.h5ad"
        split.parent.mkdir(parents=True)
        split.write_text("Preflight must not open this file.\n")
    duplicate = directory / "inputs/cd14-monocyte/duplicate.h5ad"
    duplicate.write_text("Duplicate slug.\n")
    (directory / "analyses/cd14-monocyte/analysis").mkdir(parents=True)
    splits = str(directory / "inputs/*/split.h5ad")
    analyses = str(directory / "analyses/*/analysis")
    cases = [
        ("empty", ["--cell_types", " , "], "must contain at least one"),
        ("unknown", ["--cell_type_input", splits, "--cell_types", "typo"], "absent from the input: typo. Available:"),
        ("mixed", ["--cell_type_input", splits, "--cell_types", "cd14-monocyte,typo"], "absent from the input: typo. Available:"),
        ("duplicate", ["--cell_type_input", str(duplicate.parent / "*.h5ad")], "Duplicate cell-type inputs"),
        ("missing-analysis", ["--cell_type_input", splits, "--cell_type_analysis_input", analyses], "Missing analyses: [naive-cd4-t-cell]"),
        ("missing-split", ["--cell_type_input", str(directory / "inputs/naive-cd4-t-cell/split.h5ad"), "--cell_type_analysis_input", analyses], "missing splits: [cd14-monocyte]"),
        ("split-only", ["--split_only", "--cell_types", "cd14-monocyte"], "--cell_types and --split_only cannot be used together"),
    ]
    environment = {**os.environ, "NXF_OFFLINE": "true"}
    for name, arguments, expected in cases:
        case = directory / name
        command = [*shlex.split(nextflow), "run", "cell_type_analysis.nf", "-c", "nextflow.config",
                   "-profile", "standard", "-ansi-log", "false", "-work-dir", str(case / "work"),
                   "--outdir", str(case / "published"), *arguments]
        result = subprocess.run(command, cwd=ROOT, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        (directory / f"{name}.log").write_text(result.stdout)
        if result.returncode == 0 or expected not in result.stdout:
            raise RuntimeError(f"{name}: expected preflight rejection {expected!r}\n{result.stdout}")
        if list(case.rglob(".command.sh")) or (case / "published").exists():
            raise RuntimeError(f"{name}: preflight dispatched work or published outputs")
        print(f"Passed preflight rejection: {name}", flush=True)
    # A saved analysis is a directory; producer outputs are a list of files. Both must stage identically.
    analysis = directory / 'analyses/cd14-monocyte/analysis'
    (analysis / 'analysis.h5ad').write_text('Saved analysis must remain unchanged')
    (analysis / 'report.json').write_text(json.dumps({'status': 'complete', 'cell_type': 'CD14 monocyte'}))
    outdir = directory / 'rerender/published'
    command = [*shlex.split(nextflow), 'run', 'cell_type_analysis.nf', '-c', 'nextflow.config',
               '-profile', 'standard', '-stub-run', '-ansi-log', 'false',
               '-work-dir', str(directory / 'rerender/work'), '--outdir', str(outdir),
               '--cell_type_input', str(directory / 'inputs/cd14-monocyte/split.h5ad'),
               '--cell_type_analysis_input', str(analysis)]
    result = subprocess.run(command, cwd=ROOT, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    (directory / 'rerender.log').write_text(result.stdout)
    published = outdir / 'analysis/per_type/cd14-monocyte/analysis'
    if result.returncode or not (published / 'analysis.h5ad').is_file() or not (published / 'report.html').is_file():
        raise RuntimeError(f'Saved analysis did not stage/render correctly\n{result.stdout}')
    if (analysis / 'report.html').exists() or (analysis / 'analysis.h5ad').read_text() != 'Saved analysis must remain unchanged':
        raise RuntimeError('Rerendering wrote back into the input analysis')
    print('Passed existing-analysis directory staging and immutable-input check', flush=True)
    print(f"Selection logs: {directory}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nextflow", default="nextflow")
    check_selections(parser.parse_args().nextflow)


if __name__ == "__main__":
    main()

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def select_genes(result_dir: Path, output: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(ROOT / 'scripts/select_gene_profiles.py'),
         '--result-dir', str(result_dir), '--top-n', '2', '--output', str(output)],
        capture_output=True, text=True, check=False,
    )


def write_metadata(directory: Path, status: str) -> None:
    (directory / 'run_metadata.json').write_text(json.dumps({
        'cell_type': 'CD14 monocyte',
        'models': [
            {'model': 'per_study', 'status': 'skipped', 'reason': 'too few samples'},
            {'model': 'merged', 'status': status,
             'reason': 'design matrix is rank deficient for the available samples'},
        ],
    }))


def test_selector_emits_empty_selection_for_explicitly_skipped_merged_model(tmp_path) -> None:
    write_metadata(tmp_path, 'skipped')
    output = tmp_path / 'genes.tsv'

    result = select_genes(tmp_path, output)

    assert result.returncode == 0, result.stderr
    assert output.read_text() == ''
    assert 'Skipping gene selection for CD14 monocyte' in result.stderr
    assert 'design matrix is rank deficient' in result.stderr


@pytest.mark.parametrize('status', ['complete', 'failed', None])
def test_selector_rejects_missing_results_without_explicit_skip(tmp_path, status) -> None:
    if status is not None:
        write_metadata(tmp_path, status)
    output = tmp_path / 'genes.tsv'

    result = select_genes(tmp_path, output)

    assert result.returncode != 0
    assert 'FileNotFoundError' in result.stderr
    assert not output.exists()


def test_selector_preserves_age_ranking_and_gene_tie_breaks(tmp_path) -> None:
    (tmp_path / 'merged.csv').write_text(
        'cell_type,gene,covariate,padj\n'
        'CD14 monocyte,SEX_GENE,sex,0.0001\n'
        'CD14 monocyte,GENE_Z,age,0.01\n'
        'CD14 monocyte,GENE_A,age,0.01\n'
        'CD14 monocyte,GENE_B,age,0.02\n'
        'CD14 monocyte,GENE_A,age,0.03\n'
        'CD14 monocyte,UNTESTED,age,\n'
    )
    output = tmp_path / 'genes.tsv'

    result = select_genes(tmp_path, output)

    assert result.returncode == 0, result.stderr
    assert output.read_text() == 'CD14 monocyte\tGENE_A\nCD14 monocyte\tGENE_Z\n'


@pytest.mark.parametrize('contents, error', [
    ('cell_type,gene,covariate,padj\nCD14 monocyte,GENE_A,age,\n', 'No finite merged age-DE'),
    ('cell_type,gene\nCD14 monocyte,GENE_A\n', 'missing required columns'),
])
def test_selector_keeps_invalid_result_tables_as_errors(tmp_path, contents, error) -> None:
    (tmp_path / 'merged.csv').write_text(contents)
    output = tmp_path / 'genes.tsv'

    result = select_genes(tmp_path, output)

    assert result.returncode != 0
    assert error in result.stderr
    assert not output.exists()

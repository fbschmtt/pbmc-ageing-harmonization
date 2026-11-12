import hashlib
import json
from pathlib import Path

import pytest

from pbmc_pipeline.config import ConfigurationError
from pbmc_pipeline.preflight import validate_clock_options, validate_preflight

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture
def configuration(tmp_path):
    path = tmp_path / 'pipeline.json'
    path.write_bytes((ROOT / 'config/pipeline.json').read_bytes())
    return path


def test_preflight_preserves_configuration_bytes_and_needs_no_scientific_inputs(configuration):
    original = configuration.read_bytes()
    record = validate_preflight(configuration, sections=('processing', 'cell_type_analysis', 'differential_expression', 'integration_benchmark'))
    assert record['status'] == 'valid'
    assert record['configuration_sha256'] == hashlib.sha256(original).hexdigest()
    assert configuration.read_bytes() == original


@pytest.mark.parametrize(('section', 'key', 'value'), [
    ('differential_expression', 'alpha', 0), ('cell_type_analysis', 'min_cells', 0),
    ('processing', 'target_sum', -1), ('integration_benchmark', 'max_cells', -1),
])
def test_preflight_rejects_invalid_settings_before_data_access(configuration, section, key, value):
    document = json.loads(configuration.read_text())
    document[section][key] = value
    configuration.write_text(json.dumps(document))
    with pytest.raises(ConfigurationError):
        validate_preflight(configuration, sections=(section,))


def test_clock_preflight_ignores_unconsumed_cell_and_de_fitting_settings(configuration):
    document = json.loads(configuration.read_text())
    document.pop('cell_type_analysis')
    document['differential_expression'].pop('age_trajectory')
    configuration.write_text(json.dumps(document))
    assert validate_preflight(configuration, sections=('age_clock',))['status'] == 'valid'


@pytest.fixture
def clock_options():
    return {'n_jobs': 2, 'calibration_samples': 0, 'calibration_seed': 42, 'alphas': '0.1',
            'l1_ratios': '0.5', 'external_studies': 'one,two', 'methods': 'elastic_net'}


def test_baseline_clock_allows_zero_calibration_donors(clock_options):
    validate_clock_options(clock_options)


@pytest.mark.parametrize(('key', 'value'), [('n_jobs', 0), ('alphas', 'nan'), ('l1_ratios', '1.2'),
                                           ('external_studies', 'one,one'), ('methods', 'typo')])
def test_clock_preflight_rejects_invalid_cli_options(clock_options, key, value):
    clock_options[key] = value
    with pytest.raises(ConfigurationError):
        validate_clock_options(clock_options)


def test_core_preflight_reuses_registry_and_schema_validation_without_opening_data():
    assert validate_preflight(ROOT / 'config/pipeline.json', project_root=ROOT, sections=('processing',))['status'] == 'valid'


def test_all_clock_methods_require_calibration_if_age_aware_is_included(clock_options):
    clock_options['methods'] = ''
    with pytest.raises(ConfigurationError, match='calibration'):
        validate_clock_options(clock_options)

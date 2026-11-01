from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from pbmc_pipeline.config import AgeTrajectorySettings, read_json
from pbmc_pipeline.fraction_age import fit_fraction_age_bins

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = AgeTrajectorySettings.from_mapping(
    read_json(ROOT / 'config/pipeline.json')['differential_expression']['age_trajectory']
)


def _input() -> pd.DataFrame:
    rows = []
    for study, offset in [('a', 0), ('b', 0.1)]:
        for age in (25, 35, 45):
            for i in range(8):
                sex = i % 2
                site = (i // 2) % 2
                fraction = 0.2 + offset + (0.1 if age == 35 else 0.3 if age == 45 else 0) + 0.02 * sex + 0.01 * site
                rows.append({'study': study, 'sample': f'{study}-{age}-{i}', 'age': age,
                             'study_site': f'{study}-{site}', 'sex': 'male' if sex else 'female',
                             'n_cells': round(fraction * 1000), 'n_cells_in_sample': 1000,
                             'n_cells_in_sample_l1_parent': 1000})
    return pd.DataFrame(rows)


def test_per_study_age_effects_align_reference_and_preserve_zero_counts() -> None:
    abundance = _input()
    abundance.loc[0, 'n_cells'] = 0
    result, status = fit_fraction_age_bins(abundance, settings=SETTINGS,
                                          min_denominator_cells=10, minimum_samples_per_bin=5)
    parent = result.loc[result['denominator'] == 'l1_parent']
    assert len(parent) == 6
    assert set(parent.loc[parent['age_bin'] == '20-30', 'estimate']) == {0.0}
    study_b = parent.loc[parent['study'] == 'b'].set_index('age_bin')
    assert np.isclose(study_b.loc['30-40', 'estimate'], 0.1)
    assert np.isclose(study_b.loc['40-50', 'estimate'], 0.3)
    assert np.isclose(study_b.loc['20-30', 'reference_mean_fraction'], 0.315)
    assert pd.isna(study_b.loc['30-40', 'reference_mean_fraction'])
    assert set(parent['n_samples_bin']) == {8}
    assert all(item['status'] == 'complete' for item in status)
    assert all('sex' in item['covariates'] and 'study_site' in item['covariates'] for item in status)


def test_missing_reference_and_unsupported_bins_are_not_fabricated() -> None:
    abundance = _input()
    abundance = abundance.loc[~((abundance['study'] == 'a') & (abundance['age'] == 25))]
    abundance = abundance.loc[~((abundance['study'] == 'b') & (abundance['age'] == 45) & (abundance['sex'] == 'male'))]
    result, status = fit_fraction_age_bins(abundance, settings=replace(SETTINGS),
                                          min_denominator_cells=10, minimum_samples_per_bin=5)
    assert set(result['study']) == {'b'}
    assert '40-50' not in set(result['age_bin'])
    assert any(item['study'] == 'a' and item['status'] == 'skipped' for item in status)


def test_missing_site_is_explicit_and_rank_deficiency_skips() -> None:
    abundance = _input()
    _, missing = fit_fraction_age_bins(abundance.drop(columns='study_site'), settings=SETTINGS,
                                       min_denominator_cells=10, minimum_samples_per_bin=5)
    assert 'study_site' in missing[0]['reason']
    abundance['sex'] = np.where(abundance['age'] == 25, 'female', 'male')
    result, status = fit_fraction_age_bins(abundance, settings=SETTINGS,
                                          min_denominator_cells=10, minimum_samples_per_bin=5)
    assert result.empty
    assert all('rank-deficient' in item['reason'] for item in status)


def test_adjusted_reference_averages_model_covariates_across_study_samples() -> None:
    rng = np.random.default_rng(19)
    rows = []
    for age in (25, 35, 45):
        for index in range(30):
            bmi = int(rng.integers(20, 36))
            male = int(rng.random() < (0.2 if age == 25 else 0.8))
            cmv = int(rng.random() < (0.2 if age == 25 else 0.8))
            site = int(rng.integers(0, 2))
            fraction = 0.1 + bmi * 0.002 + male * 0.06 + cmv * 0.04 + site * 0.01
            fraction += {25: 0.0, 35: 0.05, 45: 0.12}[age]
            rows.append({
                'study': 'a', 'sample': f'{age}-{index}', 'age': age,
                'sex': 'male' if male else 'female', 'cmv': 'yes' if cmv else 'no',
                'bmi': bmi, 'study_site': str(site), 'n_cells': round(fraction * 1000),
                'n_cells_in_sample': 1000, 'n_cells_in_sample_l1_parent': 1000,
            })
    abundance = pd.DataFrame(rows)
    results, _ = fit_fraction_age_bins(
        abundance, settings=SETTINGS, min_denominator_cells=10, minimum_samples_per_bin=5,
    )
    reference = results.loc[
        (results['denominator'] == 'l1_parent') & (results['age_bin'] == '20-30')
    ].iloc[0]
    expected = (
        0.1 + 0.002 * abundance['bmi'].mean()
        + 0.06 * (abundance['sex'] == 'male').mean()
        + 0.04 * (abundance['cmv'] == 'yes').mean()
        + 0.01 * (abundance['study_site'] == '1').mean()
    )
    assert np.isclose(reference['reference_adjusted_fraction'], expected)
    assert abs(reference['reference_adjusted_fraction'] - reference['reference_mean_fraction']) > 0.015

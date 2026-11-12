"""Validate consumed configuration and clock options before workflow fan-out."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

from .config import ConfigurationError, load_analysis_configuration, load_configuration


def validate_clock_options(options: dict) -> None:
    if not isinstance(options, dict):
        raise ConfigurationError('Clock options must be an object')
    for key, minimum in [('n_jobs', 1), ('calibration_samples', 0), ('calibration_seed', 0)]:
        value = options.get(key)
        try:
            parsed = int(value)
        except (TypeError, ValueError) as error:
            raise ConfigurationError(f'clock.{key} must be an integer >= {minimum}') from error
        if isinstance(value, bool) or str(parsed) != str(value) or parsed < minimum:
            raise ConfigurationError(f'clock.{key} must be an integer >= {minimum}')
    grids = {}
    for key in ('alphas', 'l1_ratios'):
        try:
            values = [float(value.strip()) for value in str(options.get(key, '')).split(',')]
        except ValueError as error:
            raise ConfigurationError(f'clock.{key} must contain finite numbers') from error
        if not values or any(not math.isfinite(value) for value in values):
            raise ConfigurationError(f'clock.{key} must contain finite numbers')
        invalid = any(value < 0 for value in values) if key == 'alphas' else any(not 0 <= value <= 1 for value in values)
        if invalid:
            raise ConfigurationError(f'clock.{key} values are out of range')
        grids[key] = values
    external = [value.strip() for value in str(options.get('external_studies', '')).split(',')]
    if not all(external) or len(external) != len(set(external)):
        raise ConfigurationError('clock.external_studies must list unique, nonempty study IDs')
    from .age_clock import _clock_methods

    available_methods = _clock_methods(grids['alphas'], grids['l1_ratios'])
    available = {method.name for method in available_methods}
    methods = options.get('methods', '')
    selected = [value.strip() for value in str(methods).split(',')] if methods else list(available)
    if not all(selected) or set(selected) - available:
        raise ConfigurationError(f'Unknown clock methods: {sorted(set(selected)-available)}; available: {sorted(available)}')
    if int(options['calibration_samples']) == 0 and any(method.age_aware_study_offset for method in available_methods if method.name in selected):
        raise ConfigurationError('Age-aware clocks require at least one calibration sample')


def validate_preflight(config: Path, *, sections: tuple[str, ...], project_root: Path | None = None,
                       clock_options: dict | None = None) -> dict:
    if project_root is not None:
        load_configuration(project_root, config.resolve())
    load_analysis_configuration(config, sections=sections)
    if clock_options is not None:
        validate_clock_options(clock_options)
    return {'schema_version': 1, 'status': 'valid', 'sections': list(sections),
            'configuration_sha256': hashlib.sha256(config.read_bytes()).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--sections', required=True)
    parser.add_argument('--project-root', type=Path)
    parser.add_argument('--clock-options', default='{}')
    parser.add_argument('--output-dir', type=Path, default=Path('validated'))
    args = parser.parse_args()
    try:
        record = validate_preflight(args.config, sections=tuple(args.sections.split(',')),
                                    project_root=args.project_root, clock_options=json.loads(args.clock_options) or None)
    except (ConfigurationError, ValueError, OSError, KeyError) as error:
        parser.exit(2, f'Configuration preflight failed: {error}\n')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.config, args.output_dir / 'pipeline.json')
    (args.output_dir / 'preflight.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Configuration preflight passed: '+', '.join(record['sections']))


if __name__ == '__main__':
    main()

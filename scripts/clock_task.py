"""Fit one clock or aggregate only the clock artifacts from this workflow run."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from pbmc_pipeline.age_clock import ClockEligibilityError, fit_age_clocks, publish_age_clock_fit
from pbmc_pipeline.age_clock_cross_type import aggregate_cross_type_clocks
from pbmc_pipeline.age_clock_report import generate_cross_type_age_clock_report


def fit(args) -> None:
    status = {
        "cell_type": args.cell_type, "slug": args.slug, "status": "complete",
        "input": args.input_identity,
        "external_studies": [v.strip() for v in args.external_studies.split(',') if v.strip()],
    }
    try:
        fitted = fit_age_clocks(
            args.input, args.config, [v.strip() for v in args.external_studies.split(',') if v.strip()], cell_type_filter=args.cell_type,
            alphas=[float(v) for v in args.alphas.split(',')],
            l1_ratios=[float(v) for v in args.l1_ratios.split(',')],
            method_names=[v.strip() for v in args.methods.split(',') if v.strip()] if args.methods else None,
            calibration_samples=args.calibration_samples, calibration_seed=args.calibration_seed,
            n_jobs=args.n_jobs,
        )
    except ClockEligibilityError as error:
        status.update(status="ineligible", reason=str(error), exception=type(error).__name__, message=str(error))
    else:
        # Staging paths differ per task; the cohort identity belongs to the shared source.
        fitted.manifest['input'] = args.input_identity
        publish_age_clock_fit(fitted, args.output_dir)
    args.status.write_text(json.dumps(status, indent=2) + '\n')


def aggregate(args) -> None:
    args.input_root.mkdir(parents=True, exist_ok=True)
    if args.existing:
        statuses = [
            {'cell_type': json.loads(path.read_text())['cell_type'], 'status': 'complete'}
            for path in args.input_root.glob('*/ageing_clock/age_clock.json')
        ]
    else:
        statuses = [json.loads(p.read_text()) for p in args.status_dir.glob('*.json')]
    completed = [r['cell_type'] for r in statuses if r['status'] == 'complete']
    ineligible = [r for r in statuses if r['status'] == 'ineligible']
    observed = {r['cell_type'] for r in statuses}
    expected = completed if args.existing else json.loads(args.expected.read_text())
    failed = [{'cell_type': name, 'exception': 'WorkflowTaskFailure', 'message': 'Clock task failed; see the workflow trace.'}
              for name in expected if name not in observed]
    summary = aggregate_cross_type_clocks(args.input_root, cell_types=completed)
    if summary['methods'] and args.template:
        generate_cross_type_age_clock_report(
            report_dir=args.input_root / 'age_clock_cross_type', template_path=args.template, project_root=Path.cwd(),
        )
    manifests = sorted(args.input_root.glob('*/ageing_clock/age_clock.json'))
    cohort = json.loads(manifests[0].read_text()) if manifests else (statuses[0] if statuses else {})
    args.output.write_text(json.dumps({
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'input': cohort.get('input'), 'external_studies': cohort.get('external_studies', []),
        'kind': 'all_cell_type_age_clock_run', 'successful_cell_types': completed,
        'eligibility_skipped_cell_types': ineligible, 'failed_cell_types': failed,
        'status': 'partial' if failed else 'complete', 'cross_cell_type_ensemble': summary,
    }, indent=2) + '\n')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    fitting = modes.add_parser('fit')
    for name in ('input', 'config', 'output-dir', 'status'):
        fitting.add_argument('--' + name, type=Path, required=True)
    for name in ('cell-type', 'slug', 'input-identity', 'external-studies'):
        fitting.add_argument('--' + name, required=True)
    fitting.add_argument('--methods', default='')
    fitting.add_argument('--alphas', default='0.001,0.01,0.1,1,10')
    fitting.add_argument('--l1-ratios', default='0.1,0.5,0.9,1')
    fitting.add_argument('--calibration-samples', type=int, default=10)
    fitting.add_argument('--calibration-seed', type=int, default=42)
    fitting.add_argument('--n-jobs', type=int, default=20)
    collecting = modes.add_parser('aggregate')
    for name in ('input-root', 'status-dir', 'expected', 'output'):
        collecting.add_argument('--' + name, type=Path, required=True)
    collecting.add_argument('--template', type=Path)
    collecting.add_argument('--existing', action='store_true')
    args = parser.parse_args()
    (fit if args.mode == 'fit' else aggregate)(args)


if __name__ == '__main__':
    main()

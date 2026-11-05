"""Combine cell-type clock predictions for shared donors and compare residuals."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

LOGGER = logging.getLogger(__name__)
MIN_CELL_TYPES = 2


def _read_predictions(path: Path, cell_type: str) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    required = {"study", "age", "predicted_age"}
    if not required.issubset(frame.columns):
        raise ValueError(f"{path} is missing prediction columns: {sorted(required - set(frame.columns))}")
    frame = frame.copy()
    frame["study"] = frame["study"].astype(str)
    frame["age"] = pd.to_numeric(frame["age"], errors="coerce")
    frame["predicted_age"] = pd.to_numeric(frame["predicted_age"], errors="coerce")
    frame = frame.dropna(subset=["age", "predicted_age"])
    sample_column = next((name for name in ("sample_id", "sample") if name in frame), None)
    if sample_column:
        sample = frame[sample_column].astype("string")
        frame["donor_key"] = [
            f"{study}::{value}" if pd.notna(value) and value.strip() else str(identifier).rsplit("::", 1)[0]
            for study, value, identifier in zip(frame.study, sample, frame.index, strict=True)
        ]
    else:
        frame["donor_key"] = [str(identifier).rsplit("::", 1)[0] for identifier in frame.index]
    frame["cell_type"] = cell_type
    if frame["donor_key"].duplicated().any():
        raise ValueError(f"{path} has multiple predictions for a donor within {cell_type}")
    return frame.loc[:, ["donor_key", "study", "age", "predicted_age", "cell_type"]]


def _metrics(frame: pd.DataFrame) -> dict:
    observed = frame["age"].to_numpy(float)
    predicted = frame["predicted_age"].to_numpy(float)
    error = predicted - observed
    by_study = frame.assign(abs_error=np.abs(error)).groupby("study", observed=True).abs_error.mean()
    return {
        "n_donors": len(frame),
        "n_studies": int(frame.study.nunique()),
        "mean_cell_types_per_donor": float(frame.n_cell_types.mean()),
        "min_cell_types_per_donor": int(frame.n_cell_types.min()),
        "max_cell_types_per_donor": int(frame.n_cell_types.max()),
        "mae_years": float(mean_absolute_error(observed, predicted)),
        "rmse_years": float(np.sqrt(mean_squared_error(observed, predicted))),
        "r_squared": float(r2_score(observed, predicted)) if len(frame) > 1 else None,
        "mean_study_mae_years": float(by_study.mean()),
        "bias_years": float(error.mean()),
    }


def _combine_predictions(frames: list[pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    per_donor: dict[str, list[dict]] = {}
    for frame in frames:
        for row in frame.to_dict("records"):
            per_donor.setdefault(row["donor_key"], []).append(row)
    ensemble_rows = []
    residual_rows = []
    rejected_metadata = 0
    for donor_key, rows in per_donor.items():
        if len(rows) < MIN_CELL_TYPES:
            continue
        if len({row["study"] for row in rows}) != 1 or not np.allclose(
            [row["age"] for row in rows], rows[0]["age"], rtol=0, atol=1e-8
        ):
            rejected_metadata += 1
            continue
        ensemble_rows.append({
            "donor_key": donor_key,
            "study": rows[0]["study"],
            "age": rows[0]["age"],
            "predicted_age": float(np.mean([row["predicted_age"] for row in rows])),
            "n_cell_types": len(rows),
            "cell_types": ";".join(sorted(row["cell_type"] for row in rows)),
        })
        for row in rows:
            residual_rows.append({
                "donor_key": donor_key,
                "study": row["study"],
                "cell_type": row["cell_type"],
                "residual_years": row["predicted_age"] - row["age"],
            })
    if rejected_metadata:
        LOGGER.warning("Excluded %d donors with inconsistent age or study across cell types", rejected_metadata)
    ensemble = pd.DataFrame(ensemble_rows)
    if not ensemble.empty:
        ensemble["age_error_years"] = ensemble["predicted_age"] - ensemble["age"]
        ensemble = ensemble.set_index("donor_key")
        ensemble.index.name = "donor_key"
        ensemble["n_cell_types"] = ensemble["n_cell_types"].astype(int)
    residuals = pd.DataFrame(residual_rows)
    correlations = []
    if not residuals.empty:
        wide = residuals.pivot(index="donor_key", columns="cell_type", values="residual_years")
        for first_index, first in enumerate(wide.columns):
            for second in wide.columns[first_index + 1:]:
                pair = wide[[first, second]].dropna()
                correlations.append({
                    "cell_type_a": first,
                    "cell_type_b": second,
                    "n_donors": len(pair),
                    "pearson_r": (
                        float(pair[first].corr(pair[second]))
                        if len(pair) > 1 and pair[first].std() > 0 and pair[second].std() > 0
                        else None
                    ),
                })
    return ensemble, residuals, pd.DataFrame(correlations)


def aggregate_cross_type_clocks(
    output_dir: Path, cell_types: list[str] | None = None,
) -> dict:
    """Write per-method mean predictions and pairwise residual correlations across types."""
    manifests = []
    requested = set(cell_types) if cell_types is not None else None
    for path in sorted(output_dir.glob("*/ageing_clock/age_clock.json")):
        manifest = json.loads(path.read_text())
        if requested is not None and manifest.get("cell_type") not in requested:
            continue
        manifests.append((path.parent, manifest))
    cohorts = {
        json.dumps({
            key: manifest.get(key)
            for key in (
                "input", "configuration_sha256", "training_studies", "external_studies",
                "cross_validation", "normalization",
            )
        }, sort_keys=True)
        for _, manifest in manifests
    }
    if len(cohorts) > 1:
        raise ValueError(
            "Per-cell-type clock outputs have different input/configuration or study splits; "
            "rerun `make run-age-clocks` before aggregating them."
        )
    method_records: dict[str, list[tuple[Path, dict, str]]] = {}
    for clock_dir, manifest in manifests:
        for record in manifest.get("methods", []):
            if record.get("method"):
                method_records.setdefault(record["method"], []).append(
                    (clock_dir, record, str(manifest.get("cell_type", record.get("cell_type", ""))))
                )

    output_root = output_dir / "age_clock_cross_type"
    output_root.mkdir(parents=True, exist_ok=True)
    summary = {
        "kind": "cross_cell_type_age_clock_ensemble",
        "ensemble": "equal-weight mean of per-cell-type predictions for each donor",
        "minimum_cell_types_per_donor": MIN_CELL_TYPES,
        "cell_types": sorted(manifest.get("cell_type", "") for _, manifest in manifests),
        "methods": [],
        "residual_correlations_file": "residual_correlations.csv",
    }
    all_correlations = []
    for method_id, entries in sorted(method_records.items()):
        if len(entries) < MIN_CELL_TYPES:
            continue
        method_summary = {
            "id": method_id,
            "name": entries[0][1].get("method_name", method_id),
            "description": "Equal-weight mean of this clock method's predictions across available cell types for each donor.",
            "training": None,
            "held_out": None,
            "component_metrics": {},
        }
        for split, filename, destination in (
            ("training", "cv_oof_predictions.csv", "training"),
            ("held_out", "external_predictions.csv", "held_out"),
        ):
            frames = []
            for clock_dir, record, cell_type in entries:
                path = clock_dir / method_id / filename
                if path.is_file():
                    frames.append(_read_predictions(path, cell_type))
            if len(frames) < MIN_CELL_TYPES:
                continue
            ensemble, _residuals, correlations = _combine_predictions(frames)
            if ensemble.empty:
                continue
            out_dir = output_root / method_id
            out_dir.mkdir(parents=True, exist_ok=True)
            predictions_file = f"{method_id}/{destination}_predictions.csv"
            ensemble.to_csv(output_root / predictions_file)
            method_summary[destination] = {
                "predictions_file": predictions_file,
                "metrics": _metrics(ensemble),
                "cell_types": sorted(frame.cell_type.iloc[0] for frame in frames),
            }
            for row in correlations.to_dict("records"):
                all_correlations.append({"split": destination, "method": method_id, **row})
            for frame in frames:
                cell_type = str(frame.cell_type.iloc[0])
                common = ensemble.index.intersection(frame.donor_key)
                if not len(common):
                    continue
                component = frame.set_index("donor_key").loc[common]
                mean_rows = ensemble.loc[common]
                method_summary["component_metrics"].setdefault(destination, []).append({
                    "cell_type": cell_type,
                    "n_donors": len(common),
                    "clock_mae_years": float(np.abs(component.predicted_age - component.age).mean()),
                    "mean_ensemble_mae_years": float(np.abs(mean_rows.predicted_age - mean_rows.age).mean()),
                })
        if method_summary["training"] and method_summary["held_out"]:
            summary["methods"].append(method_summary)
    correlation_columns = ["split", "method", "cell_type_a", "cell_type_b", "n_donors", "pearson_r"]
    pd.DataFrame(all_correlations, columns=correlation_columns).to_csv(
        output_root / "residual_correlations.csv", index=False
    )
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    LOGGER.info(
        "Cross-cell-type clock analysis includes %d cell types and %d shared clock methods",
        len(summary["cell_types"]), len(summary["methods"]),
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cell-types", help="Optional comma-separated cell-type labels to include")
    parser.add_argument("--report-template", type=Path)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    requested = [value.strip() for value in args.cell_types.split(",") if value.strip()] if args.cell_types else None
    summary = aggregate_cross_type_clocks(args.output_dir, cell_types=requested)
    LOGGER.info("Wrote cross-cell-type clock summary with %d methods", len(summary["methods"]))
    if args.report_template and summary["methods"]:
        from .age_clock_report import generate_cross_type_age_clock_report

        report = generate_cross_type_age_clock_report(
            report_dir=args.output_dir / "age_clock_cross_type",
            template_path=args.report_template,
            project_root=args.project_root,
        )
        LOGGER.info("Wrote cross-cell-type clock diagnostics: %s", report)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()

"""Validate and write provenance for one downstream cell-type analysis run."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_cell_type_manifest(
    *,
    input_path: Path,
    config_path: Path,
    model_paths: list[Path],
    report_dirs: list[Path],
    output_path: Path,
    split_manifest_path: Path | None = None,
    analysed_slugs: list[str] | None = None,
    selected_slugs: list[str] | None = None,
) -> dict:
    config = json.loads(config_path.read_text())
    analysis = config["cell_type_analysis"]
    types = []
    for report_dir in sorted(report_dirs):
        required = [report_dir / name for name in ("analysis.h5ad", "report.json", "report.html", "executed.ipynb")]
        absent = [str(path) for path in required if not path.exists()]
        if absent:
            raise FileNotFoundError(f"{report_dir}: required report outputs are absent: {absent}")
        report = json.loads((report_dir / "report.json").read_text())
        if report.get("status") not in {"complete", "insufficient_cells"}:
            raise ValueError(f"{report_dir}: unsupported report status {report.get('status')!r}")
        trajectory_files = {}
        for name in ("age_trajectory_clusters.csv", "age_trajectory_cluster_means.csv"):
            path = report_dir / name
            if path.is_file():
                trajectory_files[name] = {"path": name, "sha256": _sha256(path)}
        types.append({
            "slug": report_dir.name,
            "cell_type": report["cell_type"],
            "aifi_l1_parent": report["aifi_l1_parent"],
            "status": report["status"],
            "n_cells": report["n_cells"],
            "n_clusters": report.get("n_clusters"),
            "analysis_h5ad_sha256": _sha256(report_dir / "analysis.h5ad"),
            "report_html_sha256": _sha256(report_dir / "report.html"),
            "age_trajectory_outputs": trajectory_files,
        })
    expected_types = types
    if split_manifest_path is not None:
        split_manifest = json.loads(split_manifest_path.read_text())
        expected_types = split_manifest.get("cell_types", [])
    expected_slugs = [record["slug"] for record in expected_types]
    if not expected_slugs or len(set(expected_slugs)) != len(expected_slugs):
        raise ValueError("Expected cell types must be non-empty and have unique slugs")
    if selected_slugs:
        selected = set(selected_slugs)
        missing = selected - set(expected_slugs)
        if missing:
            raise ValueError(
                f"Requested cell types are absent from the input: {sorted(missing)}. "
                f"Available: {sorted(expected_slugs)}"
            )
        expected_types = [record for record in expected_types if record.get("slug") in selected]
    successful_slugs = {record["slug"] for record in types}
    if len(successful_slugs) != len(types):
        raise ValueError("Report directories must have unique cell-type slugs")
    unexpected = successful_slugs - {record["slug"] for record in expected_types}
    if unexpected:
        raise ValueError(f"Reports are outside the expected cell-type selection: {sorted(unexpected)}")
    analysed = set(analysed_slugs or ())
    failed_types = [
        {
            "slug": str(record["slug"]),
            "cell_type": str(record.get("cell_type") or record["slug"]),
            "failed_stage": "report" if record["slug"] in analysed else "analysis",
        }
        for record in expected_types
        if record.get("slug") not in successful_slugs
    ]
    analysis_bytes = json.dumps(analysis, sort_keys=True, separators=(",", ":")).encode()
    document = {
        "schema_version": 2,
        "kind": "cell_type_analysis",
        "status": "partial" if failed_types else "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "merged_input": {"path": str(input_path), "sha256": _sha256(input_path)},
        "configuration": {"path": str(config_path), "sha256": _sha256(config_path)},
        "analysis_specification": analysis,
        "analysis_specification_sha256": hashlib.sha256(analysis_bytes).hexdigest(),
        "models": [{"path": str(path), "sha256": _sha256(path)} for path in model_paths],
        "cell_types": types,
        "failed_cell_types": failed_types,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(document, indent=2) + "\n")
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate and write a cell-type analysis manifest")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--model", type=Path, nargs="+", required=True)
    parser.add_argument("--split-manifest", type=Path)
    parser.add_argument("--analysed-slug", action="append", default=[])
    parser.add_argument("--selected-slug", action="append", default=[])
    parser.add_argument("--report-dir", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_cell_type_manifest(
        input_path=args.input, config_path=args.config, model_paths=args.model,
        report_dirs=args.report_dir, output_path=args.output,
        split_manifest_path=args.split_manifest, analysed_slugs=args.analysed_slug,
        selected_slugs=args.selected_slug,
    )


if __name__ == "__main__":
    main()

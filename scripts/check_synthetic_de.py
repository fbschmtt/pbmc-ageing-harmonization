"""Assert the synthetic DE workflow produced usable model and report outputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import anndata as ad


def _load_fixture(outdir: Path) -> tuple[list[str], list[str]]:
    fixture_path = outdir / "synthetic_de" / "pseudobulk_merged.h5ad"
    if not fixture_path.is_file():
        raise SystemExit(f"Missing generated synthetic fixture: {fixture_path}")
    fixture = ad.read_h5ad(fixture_path, backed="r")
    expected_genes = fixture.uns["synthetic_test_data"].get("expected_age_associated_genes")
    if expected_genes is None or len(expected_genes) == 0:
        raise SystemExit(f"Fixture does not declare planted age markers: {fixture_path}")
    cell_types = sorted(fixture.obs["aifi_l2_majority"].astype(str).unique())
    fixture.file.close()
    return cell_types, list(expected_genes)


def check_results(outdir: Path) -> None:
    cell_types, expected_genes = _load_fixture(outdir)
    de_dir = outdir / "differential_expression"
    manifest_path = de_dir / "differential_expression.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing DE manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise SystemExit(f"Expected complete synthetic DE run; got {manifest.get('status')!r}")
    if manifest.get("models_completed") != 8 or manifest.get("models_failed") != 0:
        raise SystemExit(
            "Expected eight completed fits and zero failed fits; got "
            f"{manifest.get('models_completed')} completed, {manifest.get('models_failed')} failed"
        )
    if not manifest.get("synthetic_test_data") or manifest.get("test_mode"):
        raise SystemExit("Expected synthetic labeling with production inclusion filters enabled")
    if not manifest.get("sample_inclusion_applied", {}).get("minimum_cells_filter_enabled"):
        raise SystemExit("The synthetic run did not apply the configured minimum-cell filter")

    records = {record["cell_type"]: record for record in manifest.get("cell_types", [])}
    if len(cell_types) != 2 or len(records) != len(cell_types):
        raise SystemExit("Expected exactly two synthetic cell types in the DE manifest")
    for cell_type in cell_types:
        record = records.get(cell_type)
        if record is None:
            raise SystemExit(f"Missing DE manifest record for {cell_type!r}")
        models = record.get("models", [])
        if len(models) != 4 or any(model.get("status") != "complete" for model in models):
            raise SystemExit(f"Expected four completed models for {cell_type!r}")
        merged = next((model for model in models if model.get("model") == "merged"), None)
        if merged is None or not merged.get("results"):
            raise SystemExit(f"Missing merged-model result for {cell_type!r}")
        result_path = de_dir / merged["results"]
        with result_path.open(newline="") as handle:
            results = {row["gene"]: row for row in csv.DictReader(handle)}
        for gene in expected_genes:
            row = results.get(gene)
            if row is None or not row.get("padj"):
                raise SystemExit(f"Missing planted marker {gene!r} in {result_path}")
            if float(row["padj"]) >= 0.05:
                raise SystemExit(f"Planted marker {gene!r} is not significant in {result_path}")
    print("Synthetic DE check passed: eight fits completed and planted age markers were detected.")


def check_reports(outdir: Path) -> None:
    cell_types, expected_genes = _load_fixture(outdir)
    for cell_type in cell_types:
        slug = "-".join(cell_type.lower().split())
        report_path = outdir / "cell_type_analysis" / slug / "report.html"
        if not report_path.is_file():
            raise SystemExit(f"Missing standard cell-type report: {report_path}")
        html = report_path.read_text(errors="replace").lower()
        required = ["synthetic test data", "differential expression", "volcano", *expected_genes]
        missing = [item for item in required if item.lower() not in html]
        de_section = html.split("pseudobulk differential expression", maxsplit=1)
        if len(de_section) < 2 or "<img" not in de_section[1]:
            missing.append("rendered DE plot image")
        if missing:
            raise SystemExit(f"{report_path} is missing expected DE content: {missing}")
    print("Synthetic DE report check passed: warnings, marker results, and volcano plots are present.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("results", "reports"))
    parser.add_argument("--outdir", type=Path, default=Path("output/test"))
    args = parser.parse_args()
    if args.phase == "results":
        check_results(args.outdir)
    else:
        check_reports(args.outdir)


if __name__ == "__main__":
    main()

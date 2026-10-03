"""Assert the synthetic DE workflow produced usable model and report outputs."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import anndata as ad


def _load_fixture(outdir: Path) -> tuple[list[str], dict[str, list[str]], dict[str, list[str]]]:
    fixture_path = outdir / "synthetic_de" / "pseudobulk_merged.h5ad"
    if not fixture_path.is_file():
        raise SystemExit(f"Missing generated synthetic fixture: {fixture_path}")
    fixture = ad.read_h5ad(fixture_path, backed="r")
    metadata = fixture.uns["synthetic_test_data"]
    expected_genes = {
        covariate: metadata.get(f"expected_{covariate}_associated_genes")
        for covariate in ("age", "sex", "bmi", "cmv")
    }
    if any(genes is None or len(genes) == 0 for genes in expected_genes.values()):
        raise SystemExit(f"Fixture does not declare planted covariate markers: {fixture_path}")
    covariate_studies = metadata.get("expected_covariate_studies")
    if not covariate_studies:
        raise SystemExit(f"Fixture does not declare expected covariate study coverage: {fixture_path}")
    cell_types = sorted(fixture.obs["aifi_l2_majority"].astype(str).unique())
    fixture.file.close()
    return (
        cell_types,
        {key: list(value) for key, value in expected_genes.items()},
        {key: list(value) for key, value in covariate_studies.items()},
    )


def check_results(outdir: Path) -> None:
    cell_types, expected_genes, covariate_studies = _load_fixture(outdir)
    de_dir = outdir / "differential_expression"
    manifest_path = de_dir / "differential_expression.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing DE manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise SystemExit(f"Expected complete synthetic DE run; got {manifest.get('status')!r}")
    expected_per_cell_type = 3 + 1 + len(covariate_studies)
    expected_completed = expected_per_cell_type * len(cell_types)
    if manifest.get("models_completed") != expected_completed or manifest.get("models_failed") != 0:
        raise SystemExit(
            f"Expected {expected_completed} completed fits and zero failed fits; got "
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
        completed = [model for model in models if model.get("status") == "complete"]
        if len(completed) != expected_per_cell_type or any(
            model.get("status") == "failed" for model in models
        ):
            raise SystemExit(
                f"Expected {expected_per_cell_type} completed and zero failed models for {cell_type!r}; "
                f"got {len(completed)} complete, "
                f"{sum(model.get('status') == 'failed' for model in models)} failed"
            )
        merged = next((model for model in models if model.get("model") == "merged"), None)
        if merged is None or not merged.get("results"):
            raise SystemExit(f"Missing merged-model result for {cell_type!r}")
        result_path = de_dir / merged["results"]
        with result_path.open(newline="") as handle:
            results = {row["gene"]: row for row in csv.DictReader(handle)}
        for gene in expected_genes["age"]:
            row = results.get(gene)
            if row is None or not row.get("padj"):
                raise SystemExit(f"Missing planted marker {gene!r} in {result_path}")
            if float(row["padj"]) >= 0.05:
                raise SystemExit(f"Planted marker {gene!r} is not significant in {result_path}")
        sex_path = de_dir / merged["results_by_covariate"]["sex"]
        with sex_path.open(newline="") as handle:
            sex_results = {row["gene"]: row for row in csv.DictReader(handle)}
        for gene in expected_genes["sex"]:
            row = sex_results.get(gene)
            if row is None or not row.get("padj") or float(row["padj"]) >= 0.05:
                raise SystemExit(f"Planted sex marker {gene!r} is not significant in {sex_path}")

        for covariate, studies in covariate_studies.items():
            combined = next(
                (model for model in models
                 if model.get("model") == "combined" and covariate in model.get("covariates", [])),
                None,
            )
            if combined is None or combined.get("status") != "complete":
                raise SystemExit(f"Missing completed combined {covariate} fit for {cell_type!r}")
            if combined.get("studies") != studies:
                raise SystemExit(
                    f"Combined {covariate} fit for {cell_type!r} used "
                    f"{combined.get('studies')}, expected {studies}"
                )
            if covariate == "cmv" and combined.get("design") != (
                "~ study_site + age + sex + cmv"
            ):
                raise SystemExit(
                    f"Single-study CMV fit for {cell_type!r} omitted its varying site term: "
                    f"{combined.get('design')}"
                )
            expected_n = 7
            if set(combined.get("study_sample_counts", {}).values()) != {expected_n}:
                raise SystemExit(
                    f"Combined {covariate} fit for {cell_type!r} did not select the "
                    f"expected {expected_n} complete cases per study"
                )
            covariate_path = de_dir / combined["results_by_covariate"][covariate]
            with covariate_path.open(newline="") as handle:
                covariate_results = {row["gene"]: row for row in csv.DictReader(handle)}
            for gene in expected_genes[covariate]:
                row = covariate_results.get(gene)
                if row is None or not row.get("padj") or float(row["padj"]) >= 0.05:
                    raise SystemExit(
                        f"Planted {covariate} marker {gene!r} is not significant in "
                        f"{covariate_path}"
                    )

            per_study_results = [
                model for model in models
                if model.get("model") == "per_study" and covariate in model.get("covariates", [])
            ]
            completed_studies = sorted(
                model.get("study") for model in per_study_results
                if model.get("status") == "complete"
            )
            if completed_studies != studies:
                raise SystemExit(
                    f"Per-study {covariate} fits for {cell_type!r} used "
                    f"{completed_studies}, expected {studies}"
                )
        per_study_models = [model for model in models if model.get("model") == "per_study"]
        expected_optional_by_study = {
            study: sorted(
                covariate for covariate, studies in covariate_studies.items() if study in studies
            )
            for study in sorted({
                study for studies in covariate_studies.values() for study in studies
            } | {"synthetic_study_b"})
        }
        for study, optional_covariates in expected_optional_by_study.items():
            matching = [model for model in per_study_models if model.get("study") == study]
            if len(matching) != 1 or matching[0].get("status") != "complete":
                raise SystemExit(
                    f"Expected one completed maximal per-study model for {study!r}, got {matching}"
                )
            expected_covariates = ["age", "sex", *optional_covariates]
            if matching[0].get("covariates") != expected_covariates:
                raise SystemExit(
                    f"Per-study model for {study!r} used {matching[0].get('covariates')}, "
                    f"expected {expected_covariates}"
                )
    print(
        "Synthetic DE check passed: age, sex, BMI, and CMV markers were detected; "
        "covariate study selection, complete-case counts, and single-study CMV fit matched."
    )


def check_reports(outdir: Path) -> None:
    cell_types, expected_genes, _ = _load_fixture(outdir)
    for cell_type in cell_types:
        slug = "-".join(cell_type.lower().split())
        report_path = outdir / "cell_type_analysis" / slug / "report.html"
        if not report_path.is_file():
            raise SystemExit(f"Missing standard cell-type report: {report_path}")
        html = report_path.read_text(errors="replace").lower()
        required = [
            "synthetic test data", "differential expression", "shared age-model diagnostics",
        ]
        missing = [item for item in required if item.lower() not in html]
        section_heading = re.search(
            r'<h2\b[^>]*id="pseudobulk-differential-expression"[^>]*>', html
        )
        if section_heading is None:
            missing.append("DE section heading")
        else:
            de_section = html[section_heading.start():]
            for covariate in expected_genes:
                if f"{covariate} contrast:" not in de_section:
                    missing.append(f"{covariate} contrast summary")
            # The age recurrence panel and one volcano per covariate are rendered
            # as images. Gene labels inside those plots are rasterized, so their
            # names are checked against the DE result CSVs in check_results().
            expected_plots = 2 + len(expected_genes)
            if de_section.count("<img") < expected_plots:
                missing.append(
                    f"{expected_plots} DE plot images; found {de_section.count('<img')}"
                )
        if missing:
            raise SystemExit(f"{report_path} is missing expected DE content: {missing}")
    print("Synthetic DE report check passed: all covariate summaries and volcano plots are present.")


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

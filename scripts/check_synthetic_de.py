"""Assert the synthetic DE workflow produced usable model and report outputs."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def _load_fixture(
    outdir: Path,
) -> tuple[
    list[str], dict[str, list[str]], dict[str, list[str]], int, dict[str, list[float]],
    dict[str, list[str]],
]:
    fixture_path = outdir / "synthetic_de" / "pseudobulk_merged.h5ad"
    if not fixture_path.is_file():
        raise SystemExit(f"Missing generated synthetic fixture: {fixture_path}")
    fixture = ad.read_h5ad(fixture_path, backed="r")
    metadata = fixture.uns["synthetic_test_data"]
    expected_genes = {
        covariate: metadata.get(f"expected_{covariate}_associated_genes")
        for covariate in ("age", "sex", "bmi", "cmv", "age_trajectory")
    }
    if any(genes is None or len(genes) == 0 for genes in expected_genes.values()):
        raise SystemExit(f"Fixture does not declare planted covariate markers: {fixture_path}")
    covariate_studies = metadata.get("expected_covariate_studies")
    if not covariate_studies:
        raise SystemExit(f"Fixture does not declare expected covariate study coverage: {fixture_path}")
    cell_types = sorted(fixture.obs["aifi_l2_majority"].astype(str).unique())
    samples_per_study = int(metadata["samples_per_study"])
    trajectory_profiles = metadata.get("expected_age_trajectory_profiles")
    if not trajectory_profiles:
        raise SystemExit(f"Fixture does not declare planted trajectory profiles: {fixture_path}")
    per_study_age_genes = metadata.get("expected_per_study_age_associated_genes")
    if not per_study_age_genes:
        raise SystemExit(f"Fixture does not declare planted per-study age markers: {fixture_path}")
    fixture.file.close()
    return (
        cell_types,
        {key: list(value) for key, value in expected_genes.items()},
        {key: list(value) for key, value in covariate_studies.items()},
        samples_per_study,
        {gene: list(profile) for gene, profile in trajectory_profiles.items()},
        {study: list(genes) for study, genes in per_study_age_genes.items()},
    )


def check_results(outdir: Path) -> None:
    (
        cell_types, expected_genes, covariate_studies, samples_per_study,
        trajectory_profiles, per_study_age_genes,
    ) = _load_fixture(outdir)
    de_dir = outdir / "differential_expression"
    manifest_path = de_dir / "differential_expression.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing DE manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise SystemExit(f"Expected complete synthetic DE run; got {manifest.get('status')!r}")
    expected_per_cell_type = 3 + 1 + len(covariate_studies) + 1
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
        if any("log10_total_counts" not in model.get("design", "") for model in completed):
            raise SystemExit(
                f"A completed DE fit for {cell_type!r} omitted the log10_total_counts adjustment"
            )
        merged = next((model for model in models if model.get("model") == "merged"), None)
        if merged is None or not merged.get("results"):
            raise SystemExit(f"Missing merged-model result for {cell_type!r}")
        if "log10_total_counts" not in merged.get("results_by_covariate", {}):
            raise SystemExit(f"Missing combined depth coefficient for {cell_type!r}")
        depth_path = de_dir / merged["results_by_covariate"]["log10_total_counts"]
        if not depth_path.is_file():
            raise SystemExit(f"Missing combined depth results: {depth_path}")
        result_path = de_dir / merged["results"]
        with result_path.open(newline="") as handle:
            results = {row["gene"]: row for row in csv.DictReader(handle)}
        for gene in expected_genes["age"]:
            row = results.get(gene)
            if row is None or not row.get("padj"):
                raise SystemExit(f"Missing planted marker {gene!r} in {result_path}")
            if float(row["padj"]) >= 0.05:
                raise SystemExit(f"Planted marker {gene!r} is not significant in {result_path}")
        for study, genes in per_study_age_genes.items():
            per_study_model = next(
                (model for model in models
                 if model.get("model") == "per_study" and model.get("study") == study),
                None,
            )
            if per_study_model is None or "age" not in per_study_model.get(
                "results_by_covariate", {}
            ):
                raise SystemExit(
                    f"Missing per-study age result for planted marker study {study!r}"
                )
            per_study_path = de_dir / per_study_model["results_by_covariate"]["age"]
            with per_study_path.open(newline="") as handle:
                per_study_results = {row["gene"]: row for row in csv.DictReader(handle)}
            for gene in genes:
                row = per_study_results.get(gene)
                if row is None or not row.get("padj") or float(row["padj"]) >= 0.05:
                    raise SystemExit(
                        f"Planted per-study age marker {gene!r} is not significant in "
                        f"{per_study_path}"
                    )
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
                "~ study_site + age + sex + log10_total_counts + cmv"
            ):
                raise SystemExit(
                    f"Single-study CMV fit for {cell_type!r} omitted its varying site term: "
                    f"{combined.get('design')}"
                )
            expected_n = samples_per_study - 1
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

        trajectory_model = next(
            (model for model in models if model.get("purpose") == "age_trajectory"), None
        )
        if trajectory_model is None or trajectory_model.get("status") != "complete":
            raise SystemExit(f"Missing completed age-trajectory fit for {cell_type!r}")
        trajectory_path = de_dir / trajectory_model["results_by_covariate"]["age_bin"]
        with trajectory_path.open(newline="") as handle:
            trajectory_rows = list(csv.DictReader(handle))
        trajectory_genes = {
            row["gene"] for row in trajectory_rows
            if row.get("omnibus_padj") and float(row["omnibus_padj"]) < 0.05
        }
        for gene in expected_genes["age_trajectory"]:
            if gene not in trajectory_genes:
                raise SystemExit(
                    f"Planted non-linear age marker {gene!r} is not omnibus-significant "
                    f"in {trajectory_path}"
                )
        for gene in expected_genes["age_trajectory"]:
            rows = [row for row in trajectory_rows if row["gene"] == gene]
            if not rows or len({row["age_bin"] for row in rows}) != len(
                trajectory_model.get("retained_bins", [])
            ):
                raise SystemExit(f"Age trajectory for {gene!r} is missing retained bins")
            reference = next(row for row in rows if row["age_bin"] == "20-30")
            if float(reference["trajectory_scaled"]) != 0.0:
                raise SystemExit(f"Age trajectory for {gene!r} is not zero at 20–30")
            profile = np.asarray(trajectory_profiles[gene], dtype=float)
            retained_indices = [
                (int(age_bin.split("-", maxsplit=1)[0])
                 - int(trajectory_model["reference_bin"].split("-", maxsplit=1)[0]))
                // int(trajectory_model["bin_width_years"])
                for age_bin in trajectory_model.get("retained_bins", [])
            ]
            profile = profile[retained_indices]
            expected_shape = np.log2(profile / profile[0])
            differences = np.diff(expected_shape)
            direction_changes = np.count_nonzero(
                np.sign(differences[1:]) != np.sign(differences[:-1])
            )
            if direction_changes < 2:
                raise SystemExit(
                    f"Planted trajectory for {gene!r} must contain multiple direction changes"
                )
            ordered_rows = sorted(rows, key=lambda row: int(row["age_bin"].split("-")[0]))
            observed_shape = np.asarray(
                [float(row["trajectory_scaled"]) for row in ordered_rows], dtype=float
            )
            expected_shape = expected_shape / np.std(expected_shape, ddof=0)
            shape_correlation = float(np.corrcoef(expected_shape, observed_shape)[0, 1])
            if not np.isfinite(shape_correlation) or shape_correlation < 0.75:
                raise SystemExit(
                    f"Fitted shape for {gene!r} does not recover its planted non-linear "
                    f"trajectory (correlation={shape_correlation:.3f})"
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
            expected_covariates = ["age", "sex", "log10_total_counts", *optional_covariates]
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
    cell_types, expected_genes, _, _, _, _ = _load_fixture(outdir)
    for cell_type in cell_types:
        slug = "-".join(cell_type.lower().split())
        cell_type_dir = outdir / "cell_type_analysis" / slug
        report_path = cell_type_dir / "report.html"
        if not report_path.is_file():
            raise SystemExit(f"Missing standard cell-type report: {report_path}")
        cluster_path = cell_type_dir / "age_trajectory_clusters.csv"
        cluster_means_path = cell_type_dir / "age_trajectory_cluster_means.csv"
        if not cluster_path.is_file() or not cluster_means_path.is_file():
            raise SystemExit(
                f"Missing materialized trajectory clusters or means for {cell_type!r}"
            )
        cluster_rows = pd.read_csv(cluster_path)
        cluster_means = pd.read_csv(cluster_means_path)
        if cluster_rows.empty or cluster_means.empty:
            raise SystemExit(f"Trajectory clustering produced no export rows for {cell_type!r}")
        if not cluster_rows["trajectory_cluster"].notna().all():
            raise SystemExit(f"Trajectory cluster labels are missing for {cell_type!r}")
        if cluster_rows[["umap_1", "umap_2"]].notna().all(axis=1).sum() < 4:
            raise SystemExit(f"Per-cell-type UMAP was not materialized for {cell_type!r}")
        html = report_path.read_text(errors="replace").lower()
        required = [
            "synthetic test data", "differential gene expression",
            "shared age-model diagnostics",
        ]
        missing = [item for item in required if item.lower() not in html]
        section_heading = re.search(
            r'<h2\b[^>]*id="differential-gene-expression"[^>]*>', html
        )
        if section_heading is None:
            missing.append("DE section heading")
        else:
            next_section = re.search(r"<h2\b", html[section_heading.end():])
            de_end = (
                section_heading.end() + next_section.start()
                if next_section else len(html)
            )
            de_section = html[section_heading.start():de_end]
            toc_link = html.find('href="#differential-gene-expression"')
            if toc_link < 0 or toc_link > section_heading.start():
                missing.append("DE link in the top report contents box")

            def heading_text(match: re.Match[str]) -> str:
                content = re.sub(
                    r"<a\b[^>]*>.*?</a>", "", match.group(1), flags=re.DOTALL
                )
                return re.sub(r"<[^>]+>", "", content).strip().lower()

            h3s = list(re.finditer(r"<h3\b[^>]*>(.*?)</h3>", de_section, re.DOTALL))
            h3_names = [heading_text(match) for match in h3s]
            for heading in (
                "shared age-model diagnostics", "covariate effects",
                "per-study covariate models", "per-study gene-intersection diagnostic",
            ):
                if heading not in h3_names:
                    missing.append(f"{heading} subsection")

            rendered_h4 = list(re.finditer(r"<h4\b[^>]*>(.*?)</h4>", de_section, re.DOTALL))
            h4_names = [heading_text(match) for match in rendered_h4]
            plotted_covariates = [
                covariate for covariate in expected_genes if covariate != "age_trajectory"
            ]
            display_names = {
                "age": "age", "sex": "sex", "bmi": "bmi", "cmv": "cmv",
                "log10_total_counts": "log10(total_counts)",
            }
            covariate_order = [display_names[name] for name in plotted_covariates]
            covariate_order.append("log10(total_counts)")
            heading_positions = {}
            for covariate in covariate_order:
                match_index = next(
                    (index for index, name in enumerate(h4_names) if name == covariate),
                    None,
                )
                if match_index is None:
                    missing.append(f"{covariate} covariate subsection")
                    continue
                heading = rendered_h4[match_index]
                end = (
                    rendered_h4[match_index + 1].start()
                    if match_index + 1 < len(rendered_h4)
                    else len(de_section)
                )
                subsection = de_section[heading.end():end]
                heading_positions[covariate] = heading.start()
                if "<img" not in subsection:
                    missing.append(f"{covariate} volcano plot")
                if "<th>contrast</th>" not in subsection:
                    missing.append(f"{covariate} model details contrast row")
                if "<th>studies included (samples per study)</th>" not in subsection:
                    missing.append(f"{covariate} model details per-study sample counts")
            if (
                "sex" in heading_positions
                and "log10(total_counts)" in heading_positions
                and heading_positions["log10(total_counts)"] < heading_positions["sex"]
            ):
                missing.append("depth volcano subsection after sex")
            if de_section.count("<img") < len(covariate_order):
                missing.append(
                    f"at least {len(covariate_order)} covariate volcano images; "
                    f"found {de_section.count('<img')} total DE images"
                )

        age_heading = re.search(r'<h2\b[^>]*id="age-bin-trajectories"[^>]*>', html)
        if age_heading is None:
            missing.append("age-bin trajectory section heading")
        else:
            next_section = re.search(r"<h2\b", html[age_heading.end():])
            age_end = age_heading.end() + next_section.start() if next_section else len(html)
            age_section = html[age_heading.start():age_end]
            for phrase, label in (
                ("show trajectory clusters, umap, and merge diagnostic", "trajectory clusters"),
                (
                    "../../differential_expression/trajectory_analysis/report.html",
                    "cross-cell-type trajectory report link",
                ),
            ):
                if phrase not in age_section:
                    missing.append(label)
        if missing:
            raise SystemExit(f"{report_path} is missing expected DE content: {missing}")
    trajectory_dir = outdir / "differential_expression" / "trajectory_analysis"
    required_outputs = [
        trajectory_dir / "report.html",
        trajectory_dir / "executed.ipynb",
        trajectory_dir / "cross_cell_type_trajectory_clusters.csv",
        trajectory_dir / "cross_cell_type_cluster_means.csv",
        trajectory_dir / "gene_recurrence.csv",
        trajectory_dir / "gene_pattern_concordance.csv",
        trajectory_dir / "expression_atlas_matrix.csv",
        trajectory_dir / "expression_atlas_gene_umap.csv",
        trajectory_dir / "expression_atlas_gene_clusters.csv",
        trajectory_dir / "expression_atlas_study_support.csv",
        trajectory_dir / "expression_atlas_gene_availability.csv",
        trajectory_dir / "expression_atlas_gene_intersection_accounting.csv",
        trajectory_dir / "expression_atlas_technology_contrast.csv",
        trajectory_dir / "expression_atlas_intronic_contrast.csv",
        trajectory_dir / "cross_cell_type_residual_cluster_markers.csv",
        trajectory_dir / "cross_cell_type_residual_cluster_feature_composition.csv",
    ]
    missing_outputs = [str(path) for path in required_outputs if not path.is_file()]
    if missing_outputs:
        raise SystemExit(f"Missing cross-cell-type trajectory report outputs: {missing_outputs}")
    obsolete_type_clusters = [
        trajectory_dir / name for name in (
            "cell_type_trajectory_clusters.csv", "cell_type_cluster_means.csv",
        ) if (trajectory_dir / name).exists()
    ]
    if obsolete_type_clusters:
        raise SystemExit(
            f"Per-cell-type clusters should be exported beside their reports: "
            f"{obsolete_type_clusters}"
        )
    metadata = json.loads((trajectory_dir / "analysis_metadata.json").read_text())
    expression_atlas = metadata.get("expression_atlas", {})
    if expression_atlas.get("status") != "complete":
        raise SystemExit(f"Synthetic report did not produce an expression atlas: {expression_atlas}")
    if not expression_atlas.get("technology_available"):
        raise SystemExit("Synthetic report did not retain its technology labels")
    if len(expression_atlas.get("technology_contrast_cell_types", [])) != 2:
        raise SystemExit(
            "Synthetic report did not calculate one-versus-two study 3′/5′ contrasts "
            "for both cell types"
        )
    if len(expression_atlas.get("intronic_contrast_cell_types", [])) != 2:
        raise SystemExit(
            "Synthetic report did not calculate one-versus-one intronic contrasts "
            "for both cell types"
        )
    if expression_atlas.get("n_gene_expression_umap_genes", 0) < 4:
        raise SystemExit("Synthetic report did not materialize a gene-expression UMAP")
    residual_support = metadata.get("pearson_residual_clustering", {}).get("cross_cell_type", {})
    if residual_support.get("residual_cluster_method") != (
        "K-means on UMAP coordinates; silhouette-selected with at least 3 samples per cluster"
    ):
        raise SystemExit("Synthetic report did not use the configured residual clustering method")
    cross_trajectories = pd.read_csv(
        trajectory_dir / "cross_cell_type_trajectory_clusters.csv"
    )
    finite_umap = cross_trajectories[["umap_1", "umap_2"]].notna().all(axis=1).sum()
    if finite_umap < 4:
        raise SystemExit(
            "Synthetic report fixture did not exercise cross-cell-type trajectory UMAP"
        )
    html = (trajectory_dir / "report.html").read_text(errors="replace").lower()
    for phrase in (
        "cross-cell-type age trajectories", "recurrence", "each point is one gene",
    ):
        if phrase not in html:
            raise SystemExit(
                f"Cross-cell-type trajectory report is missing expected content {phrase!r}"
            )
    if not metadata.get("common_bins") or len(metadata.get("cell_types_in_cross_cell_type_analysis", [])) < 2:
        raise SystemExit("Synthetic report did not combine multiple cell types on shared bins")
    print(
        "Synthetic DE report check passed: per-cell-type cluster exports, cross-cell-type "
        "trajectory clusters, UMAPs, recurrence summaries, and the expression atlas are present."
    )


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

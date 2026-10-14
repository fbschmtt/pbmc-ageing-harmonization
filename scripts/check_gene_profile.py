"""Check the published artifacts from the synthetic gene-profile smoke test."""
from __future__ import annotations

import argparse
import gzip
import html
import json
import re
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd


def check_profile(profile_dir: Path, gene: str) -> None:
    required = [
        "report.html",
        "executed.ipynb",
        "fit_summary.json",
        "gene_age_result.csv",
        "exploratory_infant_adult_comparison.csv",
        "gene_sample_counts.csv",
        "mean_variance.csv",
        "model_fit_diagnostics.csv",
        "gene_residual_correlations.csv",
        "top_correlated_gene_residuals.csv",
        "top_gene_residual_correlations.csv",
        "naive_linear_model_result.csv",
        "naive_age_bin_samples.csv",
        "naive_age_bin_model_result.csv",
        "naive_age_bin_predictions.csv",
        "age_bin_residual_correlations.csv",
        "age_bin_correlated_gene_residuals.csv",
        "naive_study_slope_model_result.csv",
        "naive_study_slopes.csv",
        "naive_study_slope_predictions.csv",
        "naive_mixed_model_fixed_effects.csv",
        "naive_mixed_model_study_slopes.csv",
        "naive_mixed_model_predictions.csv",
        "fit_sample_metadata.csv",
        "vst_normalized_counts.csv.gz",
    ]
    missing = [name for name in required if not (profile_dir / name).is_file()]
    if missing:
        raise AssertionError(f"Missing gene-profile outputs: {', '.join(missing)}")

    report = (profile_dir / "report.html").read_text()
    report_text = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", report)))
    if "REPORT_TOC_PLACEHOLDER" in report:
        raise AssertionError("The HTML report still has an unfilled TOC placeholder")
    expected_title = f"<title>{gene} · CD14 monocyte | PBMC ageing</title>"
    if expected_title not in report:
        raise AssertionError(f"The HTML report title does not identify the gene and type: {expected_title}")
    for heading in (
        "PyDESeq2", "Model and gene result", "Mean–dispersion position",
        "Exploratory infant-to-adult comparison",
        "Pearson residuals", "Against Age", "Against Depth",
        "Correlated genes", "Naive normalized expression versus VST",
        "Distribution of normalized counts", "Normalization procedure",
        "Naive linear age model", "Naive age-bin model",
        "Naive age model with study-specific slopes",
        "Naive mixed age model with study fixed intercepts and random slopes",
        "Sample-level expression and age-bin levels",
        "Sample-level expression and fitted slopes", "Per-study age slopes",
        "limma-voom", "DESeq2", "TODO: implement R execution",
    ):
        if heading not in report:
            raise AssertionError(f"The HTML report is missing the {heading!r} section")
    if (
        "before model filters; includes lower-age samples" not in report_text
        or "all samples in this cell type first" not in report_text
    ):
        raise AssertionError("The age-bin all-sample normalization scale is not explained")
    if "ordinary R² is not reported for the negative-binomial PyDESeq2 model" not in report_text:
        raise AssertionError("The report should explain why PyDESeq2 has no ordinary R²")
    if (
        report.count("<h1") != 8
        or '<h1 id="PyDESeq2">' not in report
        or '<h1 id="Naive-normalization">' not in report
        or '<h1 id="Naive-linear-age-model">' not in report
        or '<h1 id="Naive-age-bin-model">' not in report
        or '<h1 id="Naive-age-model-with-study-specific-slopes">' not in report
        or '<h1 id="Naive-mixed-age-model-with-study-fixed-intercepts-and-random-slopes">' not in report
        or '<h1 id="limma-voom">' not in report
        or '<h1 id="DESeq2">' not in report
    ):
        raise AssertionError("Each fitted model should have a level-one report heading")
    for anchor in (
        "#PyDESeq2", "#Model-and-gene-result",
        "#Mean%E2%80%93dispersion-position", "#Pearson-residuals",
        "#Against-Age", "#Against-Depth", "#Correlated-genes",
        "#Naive-normalization", "#Normalization-procedure",
        "#Naive-normalized-expression-versus-VST", "#Naive-linear-age-model",
        "#Naive-age-bin-model", "#Sample-level-expression-and-age-bin-levels",
        "#Naive-age-model-with-study-specific-slopes", "#Per-study-age-slopes",
        "#Sample-level-expression-and-fitted-slopes", "#limma-voom", "#DESeq2",
        "#Naive-mixed-age-model-with-study-fixed-intercepts-and-random-slopes",
    ):
        if f'href="{anchor}"' not in report:
            raise AssertionError(f"The report TOC is missing the {anchor!r} heading link")
    notebook = nbformat.read(profile_dir / "executed.ipynb", as_version=4)
    errors = [cell for cell in notebook.cells if cell.cell_type == "code" for output in cell.get("outputs", []) if output.output_type == "error"]
    if errors:
        raise AssertionError(f"The executed notebook contains {len(errors)} error output(s)")
    code = "\n".join(
        cell.source for cell in notebook.cells if cell.cell_type == "code"
    )
    if (
        '"facet_column": "study"' not in code
        or '"color_column": "study_site"' not in code
        or "def plot_residuals(resids, x_values, x_label)" not in code
        or 'gene_samples["age"]' not in code
    ):
        raise AssertionError("Residual plot metadata or the shared age plot call is missing")
    if 'hue="study"' not in code or 'genewise_dispersions' not in code:
        raise AssertionError("Correlation colors or the gene-wise dispersion marker are missing")
    if "log10(total_counts)" not in report or 'gene_samples["log10_total_counts"]' not in code:
        raise AssertionError("Residual diagnostics are missing the depth axis")
    if (
        'residual_correlations.loc[residual_correlations["pearson_r"] > 0].nlargest' not in code
        or 'residual_correlations.loc[residual_correlations["pearson_r"] < 0].nsmallest' not in code
    ):
        raise AssertionError("The residual scatter plots should include positive and negative genes")
    for residual_column in (
        'age_bin_samples["residual"]',
        'gene_samples["study_slope_residual"]',
        'gene_samples["mixed_model_residual"]',
    ):
        if residual_column not in code:
            raise AssertionError(f"A new model is missing its residual plots: {residual_column}")
    if (
        'gene_samples["naive_log_expression"]' not in code
        or 'gene_samples["vst_count"]' not in code
        or 'colorbar.set_label("log10(total_counts)")' not in code
    ):
        raise AssertionError("The naive-expression versus VST depth-colored plot is missing")

    result = pd.read_csv(profile_dir / "gene_age_result.csv")
    if result.empty or result.loc[0, "gene"] != gene:
        raise AssertionError(f"The gene result does not contain {gene!r}")
    samples = pd.read_csv(profile_dir / "gene_sample_counts.csv")
    metadata = pd.read_csv(profile_dir / "fit_sample_metadata.csv")
    if samples.empty or len(samples) != len(metadata):
        raise AssertionError("The gene counts and fit metadata do not contain matching samples")
    for column in (
        "normalized_count", "vst_count", "pearson_residual",
        "size_factor", "fitted_normalized_mean",
        "naive_size_factor", "naive_normalized_count", "naive_log_expression",
        "naive_fitted", "naive_residual",
        "study_slope_residual", "mixed_model_residual",
    ):
        if column not in samples or not samples[column].notna().all():
            raise AssertionError(f"Gene sample output has missing {column} values")
    summary = json.loads((profile_dir / "fit_summary.json").read_text())
    if not np.isfinite(summary["geometric_mean_total_counts"]) or summary["geometric_mean_total_counts"] <= 0:
        raise AssertionError("The naive normalization geometric mean is not positive and finite")
    age_comparison = pd.read_csv(
        profile_dir / "exploratory_infant_adult_comparison.csv"
    )
    age_comparison_summary = summary["exploratory_infant_adult_comparison"]
    if age_comparison_summary is None:
        if (
            not age_comparison.empty
            or "does not change the primary adult-only DESeq2 fit" not in report_text
            or "do not include both age bins centered near 5 and 25 years" not in report_text
        ):
            raise AssertionError("Unavailable infant/adult age bins should be reported as exploratory and omitted")
    else:
        required_age_comparison = {
            "age_bin_effect_log2_count_plus_1",
            "adult_age_log2_fold_change_per_year",
            "adult_age_log2_fold_change_20_year_projection",
        }
        if age_comparison.empty or not required_age_comparison.issubset(age_comparison.columns):
            raise AssertionError("The exploratory infant/adult comparison is missing its effect scales")
        if not np.isclose(
            age_comparison.iloc[0]["adult_age_log2_fold_change_20_year_projection"],
            20 * age_comparison.iloc[0]["adult_age_log2_fold_change_per_year"],
        ):
            raise AssertionError("The adult 20-year projection does not match the per-year DESeq2 effect")
    mixed_summary = summary["mixed_model"]
    mixed_slopes = pd.read_csv(profile_dir / "naive_mixed_model_study_slopes.csv")
    mixed_effects = pd.read_csv(profile_dir / "naive_mixed_model_fixed_effects.csv")
    if (
        mixed_slopes.empty
        or len(mixed_slopes) != mixed_summary["n_studies"]
        or "age_centered" not in mixed_effects["term"].tolist()
        or not np.isfinite(mixed_summary["random_slope_sd_per_year"])
    ):
        raise AssertionError("The mixed model is missing fixed effects, study slopes, or diagnostics")
    mixed_predictions = pd.read_csv(profile_dir / "naive_mixed_model_predictions.csv")
    prediction_slopes = mixed_predictions.groupby("study", sort=True).apply(
        lambda rows: np.polyfit(rows["age"], rows["fitted_expression"], 1)[0],
        include_groups=False,
    )
    expected_slopes = mixed_slopes.set_index("study")["conditional_age_slope"].sort_index()
    if (
        not np.allclose(prediction_slopes.sort_index(), expected_slopes, rtol=1e-4, atol=1e-6)
    ):
        raise AssertionError("Mixed-model prediction lines disagree with the reported study slopes")
    diagnostics = pd.read_csv(profile_dir / "model_fit_diagnostics.csv")
    deseq2_diagnostics = diagnostics.loc[diagnostics["model"] == "PyDESeq2"].iloc[0]
    if (
        pd.notna(deseq2_diagnostics["r_squared"])
        or deseq2_diagnostics["r_squared_label"] != "Not reported for negative-binomial GLM"
    ):
        raise AssertionError("PyDESeq2 diagnostics should not report ordinary raw-count R²")
    naive_result = pd.read_csv(profile_dir / "naive_linear_model_result.csv")
    if naive_result.empty or "age" not in naive_result["term"].tolist():
        raise AssertionError("The naive linear model is missing its age coefficient")
    mean_variance = pd.read_csv(profile_dir / "mean_variance.csv")
    target = mean_variance.loc[mean_variance["gene"] == gene]
    if target.empty or not bool(target.iloc[0]["is_target_gene"]):
        raise AssertionError(f"The mean–variance output does not identify {gene!r}")
    correlations = pd.read_csv(profile_dir / "gene_residual_correlations.csv")
    if correlations.empty or gene in set(correlations["gene"]):
        raise AssertionError("Residual correlations must include other genes and exclude the target")
    if not correlations["pearson_r"].between(-1, 1).all():
        raise AssertionError("Residual correlations contain values outside [-1, 1]")
    top_summary = pd.read_csv(profile_dir / "top_gene_residual_correlations.csv")
    if len(top_summary) > min(11, len(correlations) + 1):
        raise AssertionError("The correlation table should contain the target plus up to ten genes")
    if top_summary.loc[0, "gene"] != gene or top_summary.loc[0, "pearson_r"] != 1.0:
        raise AssertionError("The correlation table must start with the target gene at r=1")
    if not top_summary[["mean_normalized_count", "dispersions"]].notna().all().all():
        raise AssertionError("The correlation table is missing gene means or dispersions")
    top_residuals = pd.read_csv(profile_dir / "top_correlated_gene_residuals.csv")
    if (
        top_residuals.empty
        or top_residuals["correlated_gene"].nunique() != len(top_summary) - 1
        or not top_summary["correlation_direction"].isin({"self", "positive", "negative"}).all()
    ):
        raise AssertionError("Sample residual output does not cover the signed correlation selections")
    age_bin_correlations = pd.read_csv(profile_dir / "age_bin_residual_correlations.csv")
    if not age_bin_correlations.empty:
        if not {"gene", "pearson_r", "correlation_direction", "padj"}.issubset(
            age_bin_correlations.columns
        ):
            raise AssertionError("Age-bin residual correlations need correlation and DE results")
        if not age_bin_correlations["pearson_r"].between(-1, 1).all():
            raise AssertionError("Age-bin residual correlations contain values outside [-1, 1]")
    with gzip.open(profile_dir / "vst_normalized_counts.csv.gz", "rt") as handle:
        if len(handle.readline().rstrip("\n").split(",")) < 2:
            raise AssertionError("The VST matrix has no gene columns")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, required=True)
    parser.add_argument("--gene", required=True)
    args = parser.parse_args()
    check_profile(args.profile_dir, args.gene)
    print(f"Gene-profile smoke outputs are valid: {args.profile_dir}")


if __name__ == "__main__":
    main()

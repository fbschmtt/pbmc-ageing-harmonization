# %%
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display

from pbmc_pipeline.report_theme import display_accessible_figure, render_report_header

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams.update({"figure.dpi": 115, "savefig.dpi": 115})
report_dir = Path(os.environ["AGE_CLOCK_REPORT_DIR"])
manifest = json.loads((report_dir / "age_clock.json").read_text())
methods = manifest.get("methods", [])
if not methods:
    raise ValueError(f"{report_dir}: age-clock experiment has no fitted methods")

display(HTML(render_report_header(
    title="Age-clock method comparison",
    eyebrow="PBMC ageing · shared external validation",
    subtitle=(
        f"{manifest['cell_type']} clocks trained on identical study partitions "
        "and evaluated on the same untouched external studies."
    ),
    metrics=[
        ("Methods", f"{len(methods):,}"),
        ("Training studies", f"{len(manifest['training_studies']):,}"),
        ("External studies", f"{len(manifest['external_studies']):,}"),
        ("Cell type", manifest["cell_type"]),
        ("Validation", "SYNTHETIC SMOKE TEST" if manifest.get("synthetic_test_data") else "EXTERNAL STUDY HOLDOUT"),
    ],
)))

if manifest.get("synthetic_test_data"):
    display(Markdown(
        "> **SYNTHETIC TEST OUTPUT — NOT BIOLOGICAL EVIDENCE.** "
        "This report checks the experiment pipeline only."
    ))

summary = []
study_rows = []
calibration_rows = []
prediction_frames = []
for method in methods:
    summary.append({
        "method": method["method_name"],
        "feature representation": method.get("feature_representation"),
        "training weighting": method["training_weighting"],
        "shared genes": method["n_genes"],
        "selected parameters": str(method["best_parameters"]),
        "LOSO mean study MAE (years)": method["best_training_cv_mae_years"],
        "external pooled MAE (years)": method["external_metrics_pooled"]["mae_years"],
        "external pooled RMSE (years)": method["external_metrics_pooled"]["rmse_years"],
        "external mean study MAE (years)": method["external_metrics_macro_study_mae_years"],
        "nonzero coefficients": method["n_nonzero_coefficients"],
        "final fit reached iteration limit": method.get("fit_diagnostics", {}).get(
            "iteration_limit_reached"
        ),
    })
    method_dir = report_dir / method["method"]
    predictions = pd.read_csv(method_dir / "external_predictions.csv", index_col=0)
    predictions["method"] = method["method_name"]
    predictions["age_error_years"] = predictions["predicted_age"] - predictions["age"]
    prediction_frames.append(predictions)
    for study, metrics in method["external_metrics_by_study"].items():
        study_rows.append({
            "method": method["method_name"],
            "study": study,
            "n": metrics["n_samples"],
            "MAE (years)": metrics["mae_years"],
            "RMSE (years)": metrics["rmse_years"],
            "bias (years)": predictions.loc[
                predictions["study"].astype(str) == study, "age_error_years"
            ].mean(),
        })
    for study, calibration in method.get("target_expression_calibration", {}).items():
        metrics = calibration.get("metrics", {})
        raw = metrics.get("none", {})
        corrected = metrics.get("age_protected_mean_shift", {})
        calibration_rows.append({
            "method": method["method_name"],
            "study": study,
            "status": calibration.get("status"),
            "calibration donors": calibration.get("n_calibration_samples"),
            "evaluation donors per split": calibration.get("n_evaluation_samples_per_repeat"),
            "repeats": calibration.get("n_repeats"),
            "raw held-out MAE (years)": raw.get("mean_mae_years"),
            "corrected held-out MAE (years)": corrected.get("mean_mae_years"),
            "MAE change (years)": (
                corrected.get("mean_mae_years", float("nan"))
                - raw.get("mean_mae_years", float("nan"))
            ),
            "corrected MAE SD across splits (years)": corrected.get(
                "sd_mae_across_repeats_years"
            ),
            "corrected mean bias (years)": corrected.get("mean_bias_years"),
        })

tuning_description = (
    "Hyperparameters are selected by leave-one-training-study-out validation; "
    "each study contributes one fold score to the tuning average."
    if manifest.get("cross_validation") != "disabled"
    else "Cross-validation is disabled; every method uses its fixed default parameters."
)
display(Markdown(
    f"Training studies: **{', '.join(manifest['training_studies'])}**. External "
    f"studies: **{', '.join(manifest['external_studies'])}**. {tuning_description} "
    "Choose a method using development-study CV before inspecting external results; "
    "using this external comparison for selection makes these studies part of model development."
))
display(pd.DataFrame(summary))
coverage_rows = [
    {"method": method["method_name"], **record}
    for method in methods
    for record in method.get("external_gene_coverage_by_study", [])
]
if coverage_rows:
    display(Markdown("### External assay coverage by method"))
    display(pd.DataFrame(coverage_rows))
display(Markdown("### External performance by study"))
study_table = pd.DataFrame(study_rows)
display(study_table)

figure, axis = plt.subplots(figsize=(8, 4.5))
sns.barplot(data=study_table, x="study", y="MAE (years)", hue="method", ax=axis)
axis.set_title("External MAE by held-out study")
axis.tick_params(axis="x", rotation=25)
figure.tight_layout()
display_accessible_figure(
    figure,
    alt_text="Bar chart comparing external mean absolute age error by study and clock method.",
)

if calibration_rows:
    display(Markdown(
        "### Target-site gene-expression calibration\n\n"
        "Each repeat randomly selects the configured number of known-age donors to estimate "
        "an additive per-gene site shift after accounting for a shared linear age "
        "trend and study-specific training intercepts. Calibration donors are "
        "excluded from these scores; raw and corrected predictions use the same "
        "remaining donors. The mean-shift correction does not estimate gene-wise "
        "scale changes."
    ))
    display(pd.DataFrame(calibration_rows))

all_predictions = pd.concat(prediction_frames, ignore_index=True)
figure, axis = plt.subplots(figsize=(7, 6))
sns.scatterplot(
    data=all_predictions, x="age", y="predicted_age", hue="study",
    style="method", alpha=0.75, s=40, ax=axis,
)
limits = [
    min(all_predictions["age"].min(), all_predictions["predicted_age"].min()),
    max(all_predictions["age"].max(), all_predictions["predicted_age"].max()),
]
axis.plot(limits, limits, color="#333333", linestyle="--", linewidth=1)
axis.set(xlim=limits, ylim=limits, xlabel="Observed age (years)",
         ylabel="Predicted age (years)", title="External predictions")
figure.tight_layout()
display_accessible_figure(
    figure,
    alt_text="External observed versus predicted age for all clock methods, with points colored by study.",
)

display(Markdown(
    "The study-balanced model gives each training study equal total weight in "
    "both gene scaling and regression fitting. This changes the training "
    "objective; it does not calibrate predictions by external study. External "
    "pooled metrics are sample-weighted, while mean-study MAE gives each "
    "external study equal weight."
))
display(Markdown(
    "The rank-based model replaces log-CPM magnitudes with each gene's "
    "within-sample percentile rank over the training-defined gene universe. "
    "Ranks are centered and scaled using training samples only. This removes "
    "expression magnitude information while preserving relative ordering; tied "
    "values (including zeros) receive average ranks."
))
if any(method.get("within_study_zscore") for method in methods):
    display(Markdown(
        "The within-study z-score model centers and scales each gene across the "
        "unlabeled samples of each study independently, then applies the fitted "
        "elastic net. Held-out predictions therefore use the expression cohort "
        "from that study to estimate its normalization. This is cohort-level "
        "prediction; a single-sample input has no within-study variation and "
        "normalizes to zero. The recorded normalization diagnostics include "
        "study sample counts and genes with zero within-study variance."
    ))

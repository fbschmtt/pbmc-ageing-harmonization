# %%
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display

from pbmc_pipeline.report_theme import display_accessible_figure, render_report_header

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams.update({"figure.dpi": 115, "savefig.dpi": 115})
report_dir = Path(os.environ["AGE_CLOCK_REPORT_DIR"])
manifest = json.loads((report_dir / "age_clock.json").read_text())
cell_type_records = manifest.get("cell_types", [])
if not cell_type_records:
    raise ValueError(f"{report_dir}: age-clock run has no completed cell-type models")

sample_count = sum(item["n_training_samples"] + item["n_external_samples"]
                   for item in cell_type_records)
gene_counts = [item["n_genes"] for item in cell_type_records]
external_studies = manifest["external_studies"]
training_studies = manifest["training_studies"]
synthetic = manifest.get("synthetic_test_data", False)
cv_enabled = manifest.get("cross_validation") != "disabled"
validation_coefficient_coverages = []
for item in cell_type_records:
    validation_coefficient_coverages.extend(
        item.get("external_gene_coverage_by_study", [])
    )
    validation_coefficient_coverages.extend(item.get("cv_mae_by_held_out_study", []))
missing_weight_fractions = [
    record["missing_coefficient_weight_fraction"]
    for record in validation_coefficient_coverages
    if record.get("missing_coefficient_weight_fraction") is not None
]
header_metrics = [
    ("Cell types", f"{len(cell_type_records):,}"),
    ("Training studies", f"{len(training_studies):,}"),
    ("External studies", f"{len(external_studies):,}"),
    ("Pseudobulks", f"{sample_count:,}"),
    ("Shared genes", f"{min(gene_counts):,}–{max(gene_counts):,}"),
    ("Validation", "SYNTHETIC SMOKE TEST" if synthetic else "EXTERNAL STUDY HOLDOUT"),
]
if missing_weight_fractions:
    header_metrics.append((
        "Mean missing coefficient mass",
        f"{np.mean(missing_weight_fractions):.1%}",
    ))
display(HTML(render_report_header(
    title=f"{cell_type_records[0]['method_name']} age clock",
    eyebrow="PBMC ageing · cross-study age prediction",
    subtitle=(
        "Continuous-age models trained on pseudobulk expression. "
        "External studies were kept out of feature selection, scaling, and tuning."
    ),
    metrics=header_metrics,
)))

if synthetic:
    display(Markdown(
        "> **SYNTHETIC TEST OUTPUT — NOT BIOLOGICAL EVIDENCE.** "
        "This report checks the training and reporting pipeline only."
    ))

missing_weight_warnings = []
for item in cell_type_records:
    for coverage in item.get("external_gene_coverage_by_study", []):
        if coverage.get("missing_coefficient_weight_warning", False):
            missing_weight_warnings.append({
                "model": item["method_name"],
                "evaluation": "Final external validation",
                "study": coverage["study"],
                "missing coefficient-weight fraction": coverage[
                    "missing_coefficient_weight_fraction"
                ],
            })
    for fold in item.get("cv_mae_by_held_out_study", []):
        if fold.get("missing_coefficient_weight_warning", False):
            missing_weight_warnings.append({
                "model": item["method_name"],
                "evaluation": "LOSO held-out training study",
                "study": fold["held_out_study"],
                "missing coefficient-weight fraction": fold[
                    "missing_coefficient_weight_fraction"
                ],
            })
if missing_weight_warnings:
    warning_threshold = manifest.get("missing_coefficient_weight_warning_fraction", 0.05)
    display(Markdown(
        "> **ASSAY-COVERAGE WARNING.** One or more validation assays are missing "
        f"more than {warning_threshold:.0%} of the fitted model's absolute coefficient "
        "weight. Their predictions set these features to zero and should be interpreted "
        "with particular caution."
    ))
    display(pd.DataFrame(missing_weight_warnings))

display(Markdown(
    f"Training studies: **{', '.join(training_studies)}**. External studies: "
    f"**{', '.join(external_studies)}**. CPM is computed from the available genes "
    "in each fold's training-defined model universe, followed by `log1p`; this "
    "prevents extra assay-only genes from changing shared-gene expression. Missing "
    "model genes are zero-filled after normalization, and cached `total_counts` "
    "metadata is ignored. This method then uses **"
    f"{cell_type_records[0].get('feature_representation', 'log-CPM expression')}**. "
    "Feature scaling is learned on training samples inside each cross-validation "
    "fold. The gene universe is selected from that fold's training studies. "
    f"{cell_type_records[0]['method_description']} "
    "Predictions are not calibrated by study."
))

# %% [markdown]
# ## Validation design and run summary
#
# When enabled, hyperparameters are selected by leave-one-training-study-out
# cross-validation. This tuning score is not an independent estimate after tuning.
# The external studies are evaluated once after the final refit.

# %%
summary_rows = []
cv_rows = []
for item in cell_type_records:
    summary_rows.append({
        "cell type": item["cell_type"],
        "genes": item["n_genes"],
        "training pseudobulks": item["n_training_samples"],
        "external pseudobulks": item["n_external_samples"],
        "nonzero coefficients": item.get("n_nonzero_coefficients"),
        "selected parameters": json.dumps(item["best_parameters"], sort_keys=True),
        "LOSO tuning MAE (years)": item.get("best_training_cv_mae_years"),
        "external MAE (years)": item["external_metrics_pooled"]["mae_years"],
        "external RMSE (years)": item["external_metrics_pooled"]["rmse_years"],
        "final fit hit iteration limit": item.get("fit_diagnostics", {}).get(
            "iteration_limit_reached", False
        ),
        "feature representation": item.get("feature_representation"),
    })
    for fold in item.get("cv_mae_by_held_out_study", []):
        cv_rows.append({"cell type": item["cell_type"], **fold})

display(pd.DataFrame(summary_rows))
if cv_enabled:
    display(Markdown(
        f"The current tuning procedure uses **{manifest.get('cross_validation', 'unknown')}** "
        "on the training studies. For each candidate parameter configuration, the "
        "pipeline is refit using only that fold's training studies. Candidate scores "
        "are mean absolute error across held-out training studies. The selected "
        "configuration is then refit using all training studies."
    ))
else:
    display(Markdown(
        "**Cross-validation is disabled.** The model uses fixed parameters; "
        "no CV scores or held-out-training-study predictions were generated for this run."
    ))
if cv_rows:
    cv_table = pd.DataFrame(cv_rows).rename(columns={
        "held_out_study": "held-out training study", "mae_years": "MAE (years)",
    })
    display(cv_table)
    figure, axis = plt.subplots(figsize=(max(8, len(cv_table) * 0.7), 4.5))
    sns.barplot(data=cv_table, x="held-out training study", y="MAE (years)",
                hue="cell type", ax=axis)
    axis.set_title("Leave-one-training-study-out tuning MAE at selected settings")
    axis.tick_params(axis="x", rotation=35)
    figure.tight_layout()
    display_accessible_figure(
        figure,
        alt_text="Bar chart of selected elastic-net tuning MAE for each held-out training study.",
    )

    convergence_rows = []
    for item in cell_type_records:
        grid_path = report_dir / "cv_grid_search.csv"
        if not grid_path.is_file():
            continue
        grid = pd.read_csv(grid_path)
        if grid.empty or "n_iteration_limit_folds" not in grid:
            continue
        grid["cell type"] = item["cell_type"]
        convergence_rows.append(grid)
    if convergence_rows:
        display(Markdown("### Convergence diagnostics"))
        convergence_grid = pd.concat(convergence_rows, ignore_index=True)
        affected = convergence_grid.loc[
            convergence_grid["n_iteration_limit_folds"].gt(0)
        ]
        if not affected.empty:
            display(Markdown(
                "> **CONVERGENCE WARNING.** Some candidate fits reached the estimator's "
                "iteration limit. Their tuning scores are shown for transparency, but may be "
                "less reliable. The plots show whether the selected setting was affected."
            ))
        if {"regression__alpha", "regression__l1_ratio", "mean_mae_years"}.issubset(
            convergence_grid.columns
        ):
            figure, axis = plt.subplots(figsize=(8.5, 5.2))
            sns.scatterplot(
                data=convergence_grid,
                x="regression__alpha", y="mean_mae_years",
                hue="n_iteration_limit_folds", style="regression__l1_ratio",
                palette="magma_r", s=85, ax=axis,
            )
            axis.set_xscale("log")
            axis.set(
                xlabel="Elastic-net alpha", ylabel="Mean LOSO MAE (years)",
                title="CV tuning performance and convergence burden",
            )
            axis.legend(
                title="Iteration-limit folds / L1 ratio", bbox_to_anchor=(1.02, 1),
                loc="upper left",
            )
            figure.tight_layout()
            display_accessible_figure(
                figure,
                alt_text=(
                    "Scatter plot of mean LOSO tuning error by elastic-net alpha, with color "
                    "showing how many held-out-study folds reached the iteration limit."
                ),
            )

    if (
        {"n_iterations", "iteration_limit_reached"}.issubset(cv_table.columns)
        and cv_table["n_iterations"].notna().any()
    ):
        figure, axis = plt.subplots(figsize=(max(8, len(cv_table) * 0.7), 4.5))
        colors = np.where(cv_table["iteration_limit_reached"], "#c44e52", "#4c72b0")
        axis.bar(cv_table["held-out training study"], cv_table["n_iterations"], color=colors)
        max_iter = max(
            (item.get("fit_diagnostics", {}).get("max_iter", 0) or 0)
            for item in cell_type_records
        )
        if max_iter:
            axis.axhline(max_iter, color="#c44e52", linestyle="--", linewidth=1,
                        label=f"Iteration limit ({max_iter:g})")
            axis.legend()
        axis.set(
            title="Selected-setting optimizer iterations by held-out study",
            xlabel="Held-out training study", ylabel="Iterations",
        )
        axis.tick_params(axis="x", rotation=35)
        figure.tight_layout()
        display_accessible_figure(
            figure,
            alt_text=(
                "Bar chart of optimizer iterations for selected-setting LOSO fits; red bars "
                "reached the iteration limit."
            ),
        )

# %% [markdown]
# ## Assay coverage: cross-validation versus external validation
#
# Orange rows are **inner LOSO cross-validation**: each development study is
# temporarily held out while selecting hyperparameters. Blue rows are **final
# external validation**: those studies were never used in feature selection,
# tuning, scaling, or fitting. Each LOSO fold has its own training-defined gene
# universe; final external validation uses the final universe from every
# development study. A missing gene is zero-filled at prediction time. The
# absolute coefficient-weight fraction identifies whether missing genes are
# important to that fitted model; values above 5% trigger the report warning.

# %%
coverage_rows = []
for item in cell_type_records:
    if cv_enabled:
        for fold in item.get("cv_mae_by_held_out_study", []):
            coverage_rows.append({
                "cell type": item["cell_type"],
                "evaluation": "Inner LOSO cross-validation",
                "study": fold["held_out_study"],
                "samples": fold["n_samples"],
                "model genes": fold["n_model_genes"],
                "zero-filled genes": fold["n_missing_genes"],
                "available genes": fold["n_available_genes"],
                "available fraction": fold["fraction_available"],
                "missing absolute coefficient weight": fold.get(
                    "missing_coefficient_absolute_weight"
                ),
                "total absolute coefficient weight": fold.get(
                    "total_coefficient_absolute_weight"
                ),
                "missing coefficient-weight fraction": fold.get(
                    "missing_coefficient_weight_fraction"
                ),
                "coefficient-weight warning": fold.get(
                    "missing_coefficient_weight_warning", False
                ),
            })
    for coverage in item.get("external_gene_coverage_by_study", []):
        coverage_rows.append({
            "cell type": item["cell_type"],
            "evaluation": "Final external validation",
            "study": coverage["study"],
            "samples": coverage["n_samples"],
            "model genes": coverage["n_model_genes"],
            "zero-filled genes": coverage["n_missing_genes"],
            "available genes": coverage["n_available_genes"],
            "available fraction": coverage["fraction_available"],
            "missing absolute coefficient weight": coverage.get(
                "missing_coefficient_absolute_weight"
            ),
            "total absolute coefficient weight": coverage.get(
                "total_coefficient_absolute_weight"
            ),
            "missing coefficient-weight fraction": coverage.get(
                "missing_coefficient_weight_fraction"
            ),
            "coefficient-weight warning": coverage.get(
                "missing_coefficient_weight_warning", False
            ),
        })

coverage_table = pd.DataFrame(coverage_rows)
if not coverage_table.empty:
    coverage_table = coverage_table.sort_values(
        ["cell type", "evaluation", "study"], kind="stable"
    )
    display(Markdown(
        "A nonzero count means that the corresponding prediction received one or "
        "more unavailable model genes as zero. The denominator is shown because a "
        "LOSO fold and the final external model can have different gene universes. "
        "The missing and total absolute coefficient weights are reported explicitly. "
        "Their ratio, missing coefficient-weight fraction, is the fraction of the "
        "fitted model's coefficient mass forced to zero."
    ))
    display(coverage_table)
    colors = {"Inner LOSO cross-validation": "#bd6a29", "Final external validation": "#176b87"}
    evaluations = [
        evaluation for evaluation in colors
        if evaluation in set(coverage_table["evaluation"])
    ]
    figure, axes = plt.subplots(
        len(evaluations), len(cell_type_records),
        figsize=(7 * len(cell_type_records), 3.8 * len(evaluations)), squeeze=False,
    )
    for row_index, evaluation in enumerate(evaluations):
        for column_index, item in enumerate(cell_type_records):
            axis = axes[row_index, column_index]
            rows = coverage_table.loc[
                (coverage_table["cell type"] == item["cell_type"])
                & (coverage_table["evaluation"] == evaluation)
            ]
            labels = rows["study"].tolist()
            values = rows["zero-filled genes"].to_numpy()
            bars = axis.bar(labels, values, color=colors[evaluation])
            for bar, missing, total in zip(
                bars, rows["zero-filled genes"], rows["model genes"], strict=True
            ):
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    f"{missing} / {total}",
                    ha="center", va="bottom", fontsize=8,
                )
            axis.set(
                title=f"{item['cell_type']}: {evaluation}", xlabel="Study",
                ylabel="Zero-filled model genes",
                ylim=(0, max(1, values.max() + 1)),
            )
            axis.tick_params(axis="x", rotation=35)
    figure.suptitle("Missing model genes by evaluation stage", y=1.01)
    figure.tight_layout()
    display_accessible_figure(
        figure,
        alt_text=(
            "Separate bar charts showing zero-filled model genes for inner LOSO "
            "cross-validation and final external validation. Labels show "
            "zero-filled genes divided by the applicable model gene universe."
        ),
    )

# %% [markdown]
# ## Training-gene intersection
#
# The final model is restricted using development studies only: a gene must be
# available in every development study. This diagnostic makes the cost of that
# portability constraint explicit. It is distinct from validation assay
# coverage above, which concerns genes unavailable only after a model universe
# has already been chosen.

# %%
intersection_rows = []
for item in cell_type_records:
    diagnostic = item.get("training_gene_intersection", {})
    for record in diagnostic.get("by_study", []):
        intersection_rows.append({
            "cell type": item["cell_type"],
            "study": record["study"],
            "genes available before intersection": record[
                "n_genes_available_before_intersection"
            ],
            "final shared model genes": diagnostic.get(
                "n_genes_in_final_model_universe"
            ),
            "genes dropped by shared-universe restriction": record[
                "n_genes_dropped_from_study_universe"
            ],
        })
intersection_table = pd.DataFrame(intersection_rows)
if not intersection_table.empty:
    display(Markdown(
        "The final model retains only genes shared across all development studies. "
        "`genes dropped` is relative to each study's own assay; it can differ by "
        "study because their initial gene universes differ."
    ))
    display(intersection_table)
    figure, axes = plt.subplots(
        1, len(cell_type_records), figsize=(7 * len(cell_type_records), 4.5), squeeze=False,
    )
    for axis, item in zip(axes.flat, cell_type_records):
        rows = intersection_table.loc[
            intersection_table["cell type"] == item["cell_type"]
        ]
        available = rows["genes available before intersection"].to_numpy()
        shared = rows["final shared model genes"].to_numpy()
        dropped = rows["genes dropped by shared-universe restriction"].to_numpy()
        axis.bar(rows["study"], shared, label="final shared model genes", color="#176b87")
        axis.bar(rows["study"], dropped, bottom=shared,
                 label="dropped by shared-universe restriction", color="#bd6a29")
        for position, total in enumerate(available):
            axis.text(position, total, str(total), ha="center", va="bottom", fontsize=8)
        axis.set(
            title=item["cell_type"], xlabel="Development study", ylabel="Genes",
            ylim=(0, max(1, available.max() * 1.08)),
        )
        axis.tick_params(axis="x", rotation=35)
        axis.legend()
    figure.suptitle("Training-gene intersection cost", y=1.02)
    figure.tight_layout()
    display_accessible_figure(
        figure,
        alt_text=(
            "Stacked bars for each development study: final shared model genes and "
            "genes removed when restricting to the shared training-gene universe."
        ),
    )

cv_oof_path = report_dir / "cv_oof_predictions.csv"
if cv_enabled and cv_oof_path.is_file():
    cv_oof = pd.read_csv(cv_oof_path, index_col="pseudobulk_id")
    cv_oof["age"] = pd.to_numeric(cv_oof["age"], errors="coerce")
    cv_oof["predicted_age"] = pd.to_numeric(cv_oof["predicted_age"], errors="coerce")
    cv_oof["absolute_error_years"] = (cv_oof["predicted_age"] - cv_oof["age"]).abs()
    cv_oof_metrics = cv_oof.groupby("study", observed=True, sort=True).agg(
        n_samples=("age", "size"),
        mae_years=("absolute_error_years", "mean"),
        bias_years=("age_error_years", "mean"),
    ).reset_index()
    display(Markdown(
        "### Leave-one-training-study-out predictions\n"
        "Each point is predicted by a fit that held out its entire study. "
        "These are tuning diagnostics: the same folds were used to select "
        "hyperparameters, so external-study performance remains the primary "
        "generalization estimate."
    ))
    display(cv_oof_metrics.rename(columns={
        "study": "held-out training study", "n_samples": "samples",
        "mae_years": "MAE (years)", "bias_years": "bias (years)",
    }))
    figure, axis = plt.subplots(figsize=(7.2, 6))
    sns.scatterplot(
        data=cv_oof, x="age", y="predicted_age", hue="study",
        alpha=0.78, s=42, ax=axis,
    )
    limits = [
        min(cv_oof["age"].min(), cv_oof["predicted_age"].min()),
        max(cv_oof["age"].max(), cv_oof["predicted_age"].max()),
    ]
    axis.plot(limits, limits, color="#333333", linestyle="--", linewidth=1,
              label="perfect prediction")
    axis.set(
        xlim=limits, ylim=limits, xlabel="Observed age (years)",
        ylabel="Predicted age (years)",
        title=f"{cell_type_records[0]['cell_type']}: held-out training-study predictions",
    )
    axis.legend(title="Held-out study", fontsize=8)
    figure.tight_layout()
    display_accessible_figure(
        figure,
        alt_text="Scatter plot of observed versus predicted age for training samples, colored by held-out study.",
    )

# %% [markdown]
# ## External validation by study
#
# Error metrics are reported separately per external study and pooled. Bias is
# mean predicted minus observed age; calibration slope and intercept describe
# the observed-age-on-predicted-age relationship. They are diagnostics only
# and are not applied to predictions.

# %%
external_rows = []
age_coverage_rows = []
prediction_tables = {}
metadata_tables = {}
for item in cell_type_records:
    type_dir = report_dir
    for coverage in item.get("training_age_coverage_by_study", []):
        age_coverage_rows.append({
            "cell type": item["cell_type"], "partition": "training",
            "study": coverage["study"], "n": coverage["n_samples"],
            "minimum age": coverage["minimum_age"],
            "median age": coverage["median_age"],
            "maximum age": coverage["maximum_age"],
        })
    predictions = pd.read_csv(type_dir / "external_predictions.csv", index_col=0)
    predictions["age"] = pd.to_numeric(predictions["age"], errors="coerce")
    predictions["predicted_age"] = pd.to_numeric(predictions["predicted_age"], errors="coerce")
    predictions["age_error_years"] = predictions["predicted_age"] - predictions["age"]
    predictions["absolute_error_years"] = predictions["age_error_years"].abs()
    prediction_tables[item["cell_type"]] = predictions
    metadata_tables[item["cell_type"]] = json.loads((type_dir / "metadata.json").read_text())
    for study, group in predictions.groupby("study", observed=True, sort=True):
        metrics = item["external_metrics_by_study"][str(study)]
        external_rows.append({
            "cell type": item["cell_type"], "study": study,
            "n": len(group),
            "age range": f"{group['age'].min():.1f}–{group['age'].max():.1f}",
            "observed mean age": group["age"].mean(),
            "predicted mean age": group["predicted_age"].mean(),
            "bias (years)": group["age_error_years"].mean(),
            "MAE (years)": group["absolute_error_years"].mean(),
            "RMSE (years)": np.sqrt(np.mean(group["age_error_years"] ** 2)),
            "R²": metrics["r_squared"],
            "calibration slope": metrics["calibration_slope"],
            "calibration intercept": metrics["calibration_intercept_years"],
        })
        age_coverage_rows.append({
            "cell type": item["cell_type"], "partition": "external",
            "study": study, "n": len(group),
            "minimum age": group["age"].min(), "median age": group["age"].median(),
            "maximum age": group["age"].max(),
        })

external_table = pd.DataFrame(external_rows)
display(external_table)
display(pd.DataFrame(age_coverage_rows))
coverage_rows = [
    {"cell type": item["cell_type"], **record}
    for item in cell_type_records
    for record in item.get("external_gene_coverage_by_study", [])
]
if coverage_rows:
    display(Markdown(
        "### External assay coverage\n"
        "The training-defined model gene universe is fixed. Missing assay genes "
        "are zero-filled and partial coverage is reported here; zero overlap is rejected. "
        "Availability metadata is used only for this check and never to select model genes."
    ))
    display(pd.DataFrame(coverage_rows))

figure, axes = plt.subplots(1, len(cell_type_records), figsize=(7 * len(cell_type_records), 5),
                            squeeze=False)
for axis, item in zip(axes.flat, cell_type_records):
    predictions = prediction_tables[item["cell_type"]]
    sns.scatterplot(data=predictions, x="age", y="predicted_age", hue="study",
                    alpha=0.75, s=38, ax=axis)
    limits = [min(predictions["age"].min(), predictions["predicted_age"].min()),
              max(predictions["age"].max(), predictions["predicted_age"].max())]
    axis.plot(limits, limits, color="#333333", linestyle="--", linewidth=1,
              label="perfect prediction")
    axis.set(xlim=limits, ylim=limits, title=item["cell_type"],
             xlabel="Observed age (years)", ylabel="Predicted age (years)")
    axis.legend(fontsize=8)
figure.suptitle("External predictions vs observed age", y=1.02)
figure.tight_layout()
display_accessible_figure(
    figure,
    alt_text="External observed age versus predicted age, with separate panels for cell types and points colored by study.",
)

figure, axes = plt.subplots(1, len(cell_type_records), figsize=(7 * len(cell_type_records), 4.5),
                            squeeze=False)
for axis, item in zip(axes.flat, cell_type_records):
    predictions = prediction_tables[item["cell_type"]]
    sns.scatterplot(data=predictions, x="age", y="age_error_years", hue="study",
                    alpha=0.75, s=38, ax=axis)
    axis.axhline(0, color="#333333", linestyle="--", linewidth=1)
    axis.set(title=item["cell_type"], xlabel="Observed age (years)",
             ylabel="Prediction error (years)")
figure.suptitle("External prediction error across age", y=1.02)
figure.tight_layout()
display_accessible_figure(
    figure,
    alt_text="External prediction error versus observed age, with separate panels for cell types and points colored by study.",
)

# %% [markdown]
# ## Error across age and sample support
#
# Decade summaries make age-dependent errors and uneven support visible.
# Bins are descriptive and do not alter the continuous-age model.

# %%
age_bin_rows = []
for cell_type, predictions in prediction_tables.items():
    lower = int(np.floor(predictions["age"].min() / 10) * 10)
    upper = int(np.ceil(predictions["age"].max() / 10) * 10) + 10
    edges = np.arange(lower, upper + 1, 10)
    if len(edges) < 2:
        edges = np.array([lower, lower + 10])
    bins = predictions.assign(age_bin=pd.cut(predictions["age"], bins=edges, right=False))
    for (study, age_bin), group in bins.groupby(["study", "age_bin"], observed=True, sort=True):
        age_bin_rows.append({
            "cell type": cell_type, "study": study, "age band": str(age_bin),
            "n": len(group), "MAE (years)": group["absolute_error_years"].mean(),
            "bias (years)": group["age_error_years"].mean(),
        })
age_bin_table = pd.DataFrame(age_bin_rows)
display(age_bin_table)
if not age_bin_table.empty:
    figure, axes = plt.subplots(1, len(cell_type_records),
                                figsize=(7 * len(cell_type_records), 4.5), squeeze=False)
    for axis, item in zip(axes.flat, cell_type_records):
        rows = age_bin_table.loc[age_bin_table["cell type"] == item["cell_type"]]
        sns.barplot(data=rows, x="age band", y="MAE (years)", hue="study", ax=axis)
        axis.set_title(item["cell_type"])
        axis.tick_params(axis="x", rotation=30)
    figure.suptitle("External MAE by age band", y=1.02)
    figure.tight_layout()
    display_accessible_figure(
        figure,
        alt_text="External mean absolute prediction error across age bands, separated by study.",
    )

# %% [markdown]
# ## Model and fitting diagnostics
#
# Coefficients are on the method's training-standardized feature scale. Their magnitude
# indicates model contribution under this representation, not causal effects.

# %%
for item in cell_type_records:
    type_dir = report_dir
    info = metadata_tables[item["cell_type"]]
    display(Markdown(f"### {item['cell_type']}"))
    design = pd.DataFrame([
        *[
            {"diagnostic": f"Selected {name}", "value": value}
            for name, value in item["best_parameters"].items()
        ],
        {"diagnostic": "Nonzero coefficients", "value": info.get("n_nonzero_coefficients")},
        {"diagnostic": "Model intercept (years)", "value": info.get("model_intercept_years")},
        {"diagnostic": "LOSO tuning MAE (years)", "value": item.get("best_training_cv_mae_years")},
    ])
    display(design)
    diagnostics = item.get("fit_diagnostics", {})
    if diagnostics:
        display(Markdown("#### Final-fit convergence diagnostics"))
        display(pd.DataFrame([
            {"diagnostic": key, "value": value} for key, value in diagnostics.items()
        ]))
        if diagnostics.get("iteration_limit_reached"):
            display(Markdown(
                "> **Iteration limit reached:** inspect the solver warnings and dual gap "
                "before interpreting the fitted model."
            ))
    inclusion = info.get("sample_inclusion_counts", {})
    if inclusion:
        display(pd.DataFrame([{"inclusion diagnostic": key, "pseudobulks": value}
                              for key, value in inclusion.items() if key != "by_study"]))
        if inclusion.get("by_study"):
            display(Markdown("#### Sample inclusion by study"))
            display(pd.DataFrame.from_dict(inclusion["by_study"], orient="index"))
    if info.get("software_versions"):
        display(Markdown("#### Software and configuration provenance"))
        display(pd.DataFrame([
            {"package": key, "version": value}
            for key, value in info["software_versions"].items()
        ]))
        display(Markdown(f"Configuration SHA-256: `{info['configuration_sha256']}`"))
    grid_path = type_dir / "cv_grid_search.csv"
    if cv_enabled and grid_path.is_file():
        grid = pd.read_csv(grid_path)
        display(Markdown(
            "#### Hyperparameter search: leave-one-training-study-out MAE\n"
            "Points show mean MAE across held-out training studies; bars show "
            "the standard deviation across studies. Lower is better."
        ))
        alpha_column = "regression__alpha"
        ratio_column = "regression__l1_ratio"
        if {alpha_column, ratio_column, "mean_mae_years", "sd_mae_years"}.issubset(grid):
            figure, axis = plt.subplots(figsize=(7.5, 4.8))
            for ratio, candidates in grid.groupby(ratio_column, sort=True):
                candidates = candidates.sort_values(alpha_column)
                axis.errorbar(
                    candidates[alpha_column], candidates["mean_mae_years"],
                    yerr=candidates["sd_mae_years"], marker="o", capsize=3,
                    label=f"l1 ratio = {ratio:g}",
                )
            selected = item["best_parameters"]
            axis.scatter(
                [selected[alpha_column]], [item["best_training_cv_mae_years"]],
                marker="*", s=180, color="black", zorder=5, label="selected setting",
            )
            axis.set_xscale("log")
            axis.set(
                xlabel="Elastic-net alpha (log scale)",
                ylabel="Mean held-out-study MAE (years)",
                title=f"{item['cell_type']}: CV score across parameter grid",
            )
            axis.legend(title="Mixing ratio", fontsize=8)
            figure.tight_layout()
            display_accessible_figure(
                figure,
                alt_text=(
                    f"Cross-validation MAE across elastic-net alpha values and mixing ratios "
                    f"for {item['cell_type']}; the selected setting is marked with a star."
                ),
            )
        show_columns = [column for column in (
            "candidate_order", *item["best_parameters"], "mean_mae_years",
            "sd_mae_years", "n_iteration_limit_folds",
        ) if column in grid]
        display(grid[show_columns].head(12))
    coefficients_path = type_dir / "coefficients.csv"
    if coefficients_path.is_file():
        coefficients = pd.read_csv(coefficients_path)
        coefficients["coefficient"] = pd.to_numeric(
            coefficients["coefficient"], errors="coerce"
        )
        coefficients = coefficients.loc[np.isfinite(coefficients["coefficient"])].copy()
        coefficients["absolute_coefficient"] = coefficients["coefficient"].abs()
        selected = coefficients.loc[coefficients["coefficient"] != 0].copy()
        n_features = len(coefficients)
        n_selected = len(selected)
        display(Markdown("#### Fitted feature selection and coefficient sizes"))
        display(Markdown(
            "A feature is counted as selected when its final-fit coefficient is "
            "exactly nonzero. L1 regularization can produce exact zeros, but this "
            "count is a property of the fitted solution, not a count of validated "
            "biological signals. Correlated genes can remain jointly selected, "
            "and their individual coefficients may be unstable across training sets."
        ))
        display(pd.DataFrame([
            {"diagnostic": "Features in final model", "value": n_features},
            {"diagnostic": "Features with nonzero coefficient", "value": n_selected},
            {"diagnostic": "Selected fraction", "value": (
                n_selected / n_features if n_features else np.nan
            )},
            {"diagnostic": "Positive coefficients", "value": int((selected["coefficient"] > 0).sum())},
            {"diagnostic": "Negative coefficients", "value": int((selected["coefficient"] < 0).sum())},
        ]))
        if n_selected:
            selected["coefficient sign"] = np.where(
                selected["coefficient"] > 0, "positive", "negative"
            )
            figure, axis = plt.subplots(figsize=(8, 4.8))
            sns.histplot(
                data=selected,
                x="absolute_coefficient",
                hue="coefficient sign",
                bins=min(30, max(8, int(np.sqrt(n_selected)))),
                multiple="layer",
                element="step",
                log_scale=True,
                ax=axis,
            )
            axis.set(
                xlabel="Absolute coefficient (years per training SD; log scale)",
                ylabel="Number of selected features",
                title=f"{item['cell_type']}: nonzero coefficient magnitudes",
            )
            figure.tight_layout()
            display_accessible_figure(
                figure,
                alt_text=(
                    f"Log-scale histogram of absolute nonzero coefficients for "
                    f"{item['cell_type']}, separated by coefficient sign."
                ),
            )
        display(Markdown(
            "#### Strongest fitted predictors\n"
            "Ranked by maximum absolute coefficient in the final full-training fit. "
            "Coefficients are years per training standard deviation of the feature; "
            "the sign gives the direction conditional on the other fitted features. "
            "These are predictive associations, not causal effects."
        ))
        display(coefficients.sort_values("absolute_coefficient", ascending=False).head(20))
        display(Markdown(
            f"Coefficient artifact: [{coefficients_path.name}]({coefficients_path.name}). "
            f"Prediction artifact: [external_predictions.csv](external_predictions.csv)."
        ))

# %% [markdown]
# ## Interpretation and limits
#
# External errors describe transfer for the selected studies and their observed
# age ranges. They do not prove generalization to every future site. Calibration
# values are shown as diagnostics only; applying them would change the stated
# no-calibration objective. The synthetic fixture is explicitly marked above.

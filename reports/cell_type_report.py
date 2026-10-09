# %%
import base64
import json
import os
import re
import warnings
from html import escape
from io import BytesIO
from pathlib import Path

# The report runs in a notebook kernel where ipywidgets is intentionally not a
# runtime dependency; tqdm's optional-progress warning is not analytically
# relevant.
warnings.filterwarnings("ignore", message="IProgress not found.*")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.api as sm
from IPython.display import HTML, Markdown, display
from matplotlib.colors import Normalize

from pbmc_pipeline.report_theme import render_report_header

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 110})
from pbmc_pipeline.config import AgeTrajectorySettings, read_json
from pbmc_pipeline.covariates import (
    categorical_contrast_label,
    categorical_levels,
    normalize_cmv_status,
)
from pbmc_pipeline.reporting import (
    plot_trajectory_merge_diagnostic,
    sample_cell_type_fractions,
    sample_cluster_fractions,
    signed_value_at_largest_absolute_magnitude,
)
from pbmc_pipeline.trajectory_analysis import (
    cluster_sample_residuals,
    cluster_trajectory_profiles,
    summarize_residual_cluster_markers,
    trajectory_merge_diagnostics,
)

primitive = sc.read_h5ad(Path(os.environ["CELL_TYPE_H5AD"]))
analysis_path = Path(os.environ["CELL_TYPE_ANALYSIS_H5AD"])
output_dir = analysis_path.parent
pipeline = read_json(Path(os.environ["CELL_TYPE_CONFIG"]))
trajectory_settings = AgeTrajectorySettings.from_mapping(
    pipeline["differential_expression"]["age_trajectory"]
)
sc.settings.verbosity = 0
report = read_json(output_dir / "report.json")
adata = sc.read_h5ad(analysis_path)
de_dir = os.environ.get("CELL_TYPE_DIFFERENTIAL_EXPRESSION_DIR")
de_run_metadata = {}
if de_dir:
    de_metadata_path = Path(de_dir) / "run_metadata.json"
    if de_metadata_path.is_file():
        de_run_metadata = read_json(de_metadata_path)
cell_type_name = str(primitive.obs["aifi_l2_majority"].iloc[0])
status_label = str(report.get("status", "unknown")).replace("_", " ").title()
hero_samples = (
    primitive.obs.groupby(["study", "sample"], observed=True)
    .agg(age=("age", "first"))
    .reset_index()
)
hero_ages = pd.to_numeric(hero_samples["age"], errors="coerce").dropna()
hero_age_range = (
    f"{hero_ages.min():g}–{hero_ages.max():g} years" if not hero_ages.empty else "Not provided"
)
summary_metrics = [
    ("Cells", f"{primitive.n_obs:,}"),
    ("Studies", f"{hero_samples['study'].nunique():,}"),
    ("Retained samples", f"{len(hero_samples):,}"),
    ("Age support", hero_age_range),
    ("Local clusters", str(report.get("n_clusters", "—"))),
    ("Analysis status", status_label),
]
display(HTML(render_report_header(
    title=cell_type_name,
    eyebrow="PBMC ageing · per-cell-type report",
    subtitle=(
        "Descriptive sample-level and single-cell diagnostics. "
        "Cell-level associations are not independent-sample inference."
    ),
    metrics=summary_metrics,
)))
provenance = pd.Series(report, name="value").to_frame()
display(HTML(
    "<details class=\"report-provenance\"><summary>Analysis provenance</summary>"
    + provenance.to_html(escape=True)
    + "</details>"
))
if de_run_metadata.get("analysis_mode") == "test_only":
    display(Markdown(
        "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
        + de_run_metadata["interpretation_warning"]
    ))


def display_collapsible_table(data, summary, *, index=False, max_height=None):
    """Render a table inside a closed-by-default details panel."""
    table_html = data.to_html(index=index, escape=True, border=0)
    if max_height is not None:
        table_html = (
            f'<div style="max-height: {int(max_height)}px; overflow: auto;" '
            'tabindex="0">' + table_html + "</div>"
        )
    display(HTML(
        f'<details class="report-details"><summary>{escape(summary)}</summary>'
        f"{table_html}</details>"
    ))


def display_collapsible_figure(figure, summary, *, open_by_default=False):
    """Render a figure in a details panel, optionally open on load."""
    image_buffer = BytesIO()
    figure.savefig(image_buffer, format="png", bbox_inches="tight")
    image_data = base64.b64encode(image_buffer.getvalue()).decode("ascii")
    plt.close(figure)
    escaped_summary = escape(summary, quote=True)
    display(HTML(
        '<details class="report-details report-details-figure"'
        + (" open" if open_by_default else "")
        + ">"
        f"<summary>{escape(summary)}</summary>"
        f'<img alt="{escaped_summary}" src="data:image/png;base64,{image_data}">'
        "</details>"
    ))


def load_residual_sample_clusters(de_dir, fit, settings):
    """Cluster one cell type's fitted age-bin residuals for its report."""
    if fit.get("status") != "complete":
        return None, None, 0, None
    residual_path = Path(de_dir) / fit.get("pearson_residuals", "")
    metadata_path = Path(de_dir) / fit.get("pearson_residual_sample_metadata", "")
    if not residual_path.is_file() or not metadata_path.is_file():
        return None, None, 0, None
    with np.load(residual_path, allow_pickle=False) as stored:
        residuals = pd.DataFrame(
            stored["residuals"],
            index=stored["pseudobulk_id"].astype(str),
            columns=stored["genes"].astype(str),
        )
    samples = pd.read_csv(metadata_path, dtype={"pseudobulk_id": str})
    samples["sample_key"] = samples["study"].astype(str) + "|" + samples["sample"].astype(str)
    sample_by_id = samples.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id")
    shared_ids = residuals.index.intersection(sample_by_id.index, sort=False)
    residuals = residuals.loc[shared_ids]
    sample_by_id = sample_by_id.loc[shared_ids]
    residuals.index = sample_by_id["sample_key"].to_numpy()
    keep = ~residuals.index.duplicated(keep="first")
    residuals = residuals.loc[keep]
    sample_by_key = sample_by_id.loc[~sample_by_id["sample_key"].duplicated()].set_index("sample_key")
    annotations, pc_scores, n_components = cluster_sample_residuals(
        residuals, sample_by_key, settings
    )
    residuals = residuals.loc[annotations["sample_key"].astype(str)].copy()
    residuals.index = annotations["sample_key"].astype(str)
    return annotations, pc_scores, n_components, residuals


def plot_residual_cluster_markers(residuals, annotations):
    """Show the residual UMAP clusters and their largest mean-residual contrasts."""
    display(Markdown("### Residual Cluster Markers"))
    coordinates = annotations.dropna(subset=["umap_1", "umap_2"]).copy()
    if coordinates.empty:
        display(Markdown("_Residual UMAP coordinates are unavailable for this cell type._"))
    else:
        cluster_counts = annotations["residual_cluster"].value_counts()
        coordinates["residual_cluster_label"] = coordinates["residual_cluster"].map(
            lambda cluster: f"Cluster {int(cluster)} (n={cluster_counts[cluster]})"
        )
        figure, axis = plt.subplots(figsize=(6.6, 4.8))
        sns.scatterplot(
            data=coordinates, x="umap_1", y="umap_2", hue="residual_cluster_label",
            palette="tab20", s=40, alpha=0.84, linewidth=0, ax=axis,
        )
        axis.set_title("Residual UMAP: sample clusters")
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        if axis.legend_ is not None:
            axis.legend_.set_title("Residual cluster")
            axis.legend_.set_bbox_to_anchor((1.02, 1))
        figure.tight_layout()
        plt.show()

    markers, _, reference_cluster = summarize_residual_cluster_markers(
        residuals, annotations
    )
    if markers.empty:
        display(Markdown(
            "_No between-cluster marker contrasts are available; at least two residual "
            "clusters are needed._"
        ))
        return

    display(Markdown(
        f"For each cluster, these are the genes with the largest lower and higher mean "
        f"Pearson residuals versus the largest cluster (cluster {reference_cluster}). "
        "The top 10 in each direction are shown. These are descriptive residual contrasts, "
        "not differential-expression significance tests."
    ))
    marker_clusters = sorted(markers["residual_cluster"].unique())
    figure, axes = plt.subplots(
        len(marker_clusters), 2, figsize=(13, 3.4 * len(marker_clusters)),
        squeeze=False,
    )
    for row_index, cluster in enumerate(marker_clusters):
        for column_index, direction in enumerate(("lower", "higher")):
            axis = axes[row_index, column_index]
            selected = markers.loc[
                (markers["residual_cluster"] == cluster)
                & (markers["direction"] == direction)
            ].sort_values("rank").head(10)
            axis.barh(
                selected["feature"].iloc[::-1],
                selected["mean_residual_difference_vs_reference"].iloc[::-1],
                color="#4c72b0",
            )
            axis.axvline(0, color="#555555", linewidth=0.8)
            axis.set_title(f"Cluster {cluster}: {direction} vs. largest cluster")
            axis.set_xlabel("Mean Pearson-residual difference")
    figure.tight_layout()
    plt.show()


def plot_residual_sample_umaps(table, *, minimum_umap_samples):
    """Show residual clusters and sample metadata over the residual UMAP."""
    coordinates = table.dropna(subset=["umap_1", "umap_2"])
    if coordinates.empty:
        display(Markdown(
            f"_At least {minimum_umap_samples} samples are required for residual UMAPs; "
            "the cluster assignments are still available below._"
        ))
        return
    def plot_covariate_grid(
        covariates, title, *, categorical_covariates, figsize, open_by_default=False
    ):
        figure, axes = plt.subplots(2, 2, figsize=figsize)
        for axis, covariate in zip(axes.flat, covariates):
            if covariate not in coordinates:
                axis.text(0.5, 0.5, "Not recorded", ha="center", va="center")
                axis.set_axis_off()
                continue
            values = coordinates.dropna(subset=[covariate]).copy()
            present = len(values)
            if values.empty:
                axis.text(0.5, 0.5, "No recorded values", ha="center", va="center")
                axis.set_axis_off()
                continue
            categorical = covariate in categorical_covariates
            sns.scatterplot(
                data=values, x="umap_1", y="umap_2", hue=covariate,
                palette="tab20" if categorical else "viridis", s=50, alpha=0.85,
                linewidth=0, ax=axis, legend="brief",
            )
            axis.set_title(f"{covariate.replace('_', ' ')} · {present}/{len(coordinates)} present")
            axis.set_xlabel("UMAP 1")
            axis.set_ylabel("UMAP 2")
            if axis.legend_ is not None:
                axis.legend_.set_bbox_to_anchor((1.02, 1))
        for axis in axes.flat[len(covariates):]:
            axis.set_axis_off()
        figure.suptitle(title, y=1.01)
        figure.tight_layout()
        display_collapsible_figure(
            figure, title, open_by_default=open_by_default,
        )

    optional_covariates = [
        covariate for covariate in ("bmi", "cmv") if covariate in coordinates
    ]
    if optional_covariates or "residual_cluster" in coordinates:
        plot_covariate_grid(
            ["residual_cluster", *optional_covariates],
            "Residual clusters and non-model covariates over residual UMAP",
            categorical_covariates={"residual_cluster", "cmv"},
            figsize=(18, 10.2),
            open_by_default=True,
        )
    plot_covariate_grid(
        ["sex", "age", "log10_total_counts", "study"],
        "Model covariates over residual UMAP",
        categorical_covariates={"sex", "study"},
        figsize=(18, 10.2),
    )


def natural_sort_key(value):
    """Sort labels naturally so numbered clusters follow numeric order."""
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(value))
    )


def plot_umap(adata, *, color, title=None, label_clusters=False, marker_area_scale=1.0):
    """Render comparable low-resolution, equal-aspect UMAP panels."""
    # Scanpy's size controls marker area, so 4x area gives about 2x diameter.
    point_size = max(2, 100_000 / adata.n_obs) * marker_area_scale
    plot_kwargs = {"color": color, "size": point_size, "show": False}
    if label_clusters:
        plot_kwargs["legend_loc"] = "on data"
    sc.pl.umap(adata, **plot_kwargs)
    figure = plt.gcf()
    figure.set_size_inches(9.5, 6.75)
    figure.set_dpi(110)
    for axis in figure.axes:
        if axis.has_data():
            axis.set_aspect("equal", adjustable="box")
    if title:
        figure.axes[0].set_title(title)
    figure.tight_layout()
    plt.show()


def has_age_span(data, *, min_samples=2):
    return len(data) >= min_samples and data["age"].nunique() >= 2


def facet_grid_shape(n_panels, *, max_columns=4):
    """Keep study panels compact while growing to any number of studies."""
    n_columns = min(max_columns, n_panels)
    return (n_panels + n_columns - 1) // n_columns, n_columns


def plot_study_linear_trends(
    data, *, fraction_column, title, ylabel, palette, show_fit_legend=True,
    open_by_default=True,
):
    """Show within-study trends with and without samples younger than 20."""
    if data.empty:
        display(Markdown(f"_No samples meet the denominator cutoff for {title}._"))
        return
    studies = sorted(data["study"].unique())
    n_rows, n_columns = facet_grid_shape(len(studies))
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(4.8 * n_columns, 3.9 * n_rows),
        squeeze=False, sharex=True, sharey=True,
    )
    for axis, study in zip(axes.flat, studies):
        subset = data.loc[data["study"] == study]
        sns.scatterplot(
            data=subset, x="age", y=fraction_column, color=palette[study],
            s=24, alpha=0.7, edgecolor="none", ax=axis,
        )
        for include_all_ages, color, linestyle in [
            (False, "#111111", "-"),
            (True, "#777777", "--"),
        ]:
            fit_data = subset if include_all_ages else subset.loc[subset["age"] >= 20]
            if has_age_span(fit_data):
                sns.regplot(
                    data=fit_data, x="age", y=fraction_column, scatter=False, ci=None,
                    color=color,
                    line_kws={"linewidth": 2.4, "alpha": 0.95, "linestyle": linestyle},
                    ax=axis,
                )
        if not has_age_span(subset):
            axis.text(0.5, 0.08, "No age span for trend", transform=axis.transAxes,
                      ha="center", va="bottom", fontsize=9)
        axis.set_title(f"{study} (n={len(subset)})")
        axis.set_ylabel(ylabel)
    for axis in axes.flat[len(studies):]:
        axis.remove()
    figure.suptitle(title, y=0.99)
    if show_fit_legend:
        fit_handles = [
            plt.Line2D([], [], color="#111111", linewidth=2.4,
                       label="Linear fit: age ≥20"),
            plt.Line2D([], [], color="#777777", linewidth=2.4, linestyle="--",
                       label="Linear fit: all samples"),
        ]
        figure.legend(handles=fit_handles, loc="upper center", ncol=2,
                      bbox_to_anchor=(0.5, 0.94))
        figure.tight_layout(rect=(0, 0, 1, 0.87))
    else:
        figure.tight_layout(rect=(0, 0, 1, 0.94))
    display_collapsible_figure(
        figure, f"Within-study trends: {title}", open_by_default=open_by_default,
    )


def plot_study_fraction_curves(fraction_specs, *, palette):
    """Plot each study's descriptive lowess curve without a pooled fit."""
    figure, axes = plt.subplots(
        1, len(fraction_specs), figsize=(7.8 * len(fraction_specs), 5.5), squeeze=False,
    )
    for axis, (data, fraction_column, label) in zip(axes.flat, fraction_specs):
        plot_data = data[["study", "age", fraction_column]].copy()
        plot_data["study"] = plot_data["study"].astype(str)
        plot_data["age"] = pd.to_numeric(plot_data["age"], errors="coerce")
        plot_data[fraction_column] = pd.to_numeric(plot_data[fraction_column], errors="coerce")
        plot_data = plot_data.dropna().reset_index(drop=True)
        if plot_data.empty:
            axis.text(
                0.5, 0.5, "No samples available for descriptive curves",
                transform=axis.transAxes, ha="center", va="center",
            )
        else:
            studies = sorted(plot_data["study"].unique())
            sns.scatterplot(
                data=plot_data, x="age", y=fraction_column, hue="study", palette=palette,
                s=26, alpha=0.2, edgecolor="none", ax=axis,
            )
            for study in studies:
                study_data = plot_data.loc[plot_data["study"] == study]
                if has_age_span(study_data):
                    sns.regplot(
                        data=study_data, x="age", y=fraction_column, scatter=False,
                        lowess=True, ci=None, color=palette[study],
                        line_kws={"linewidth": 2.4, "alpha": 0.95}, ax=axis,
                    )
            axis.legend(title="Study", loc="best", fontsize=8)
        axis.set_title(label)
        axis.set_xlabel("Age (years)")
        axis.set_ylabel(label)
    figure.suptitle("Study-specific descriptive fraction curves", y=0.98)
    figure.text(
        0.5, 0.015,
        "Points are study × sample (alpha = 0.2). Each coloured curve is a separate study lowess fit; no combined model is fitted.",
        ha="center", fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.06, 1, 0.92))
    display_collapsible_figure(figure, "Show study-specific fraction curves")


def plot_cluster_fractions(data, *, palette):
    """Keep study identity in points while showing one uncluttered pooled trend."""
    if data.empty:
        display(Markdown("_No samples meet the denominator cutoff for cluster fractions._"))
        return
    plot_data = data.assign(study=data["study"].astype(str))
    clusters = sorted(plot_data["cluster"].astype(str).unique(), key=natural_sort_key)
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(5.5 * n_columns, 4.4 * n_rows), squeeze=False,
        sharex=True, sharey=True,
    )
    studies = sorted(plot_data["study"].unique())
    for axis, cluster in zip(axes.flat, clusters):
        subset = plot_data.loc[plot_data["cluster"].astype(str) == cluster]
        sns.scatterplot(
            data=subset, x="age", y="fraction_within_cell_type", hue="study", palette=palette,
            s=20, alpha=0.65, edgecolor="none", ax=axis, legend=False,
        )
        if has_age_span(subset):
            sns.regplot(
                data=subset, x="age", y="fraction_within_cell_type", scatter=False, ci=None,
                color="#1f1f1f", line_kws={"linewidth": 2.6}, ax=axis,
            )
        axis.set_title(f"Cluster {cluster}")
        axis.set_ylabel("Fraction within cell type")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study, markersize=5)
        for study in studies
    ]
    figure.legend(handles=handles, title="Study", loc="upper center", ncol=min(4, len(studies)))
    figure.suptitle("Cluster fractions: pooled linear trends (study-unadjusted)", y=0.98)
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    display_collapsible_figure(figure, "Show cluster-composition curves")


def _bootstrap_binomial_residual_sd(design, denominators, fitted_probabilities, *, rng, n_bootstrap):
    """Estimate post-model residual SD expected from binomial cell sampling alone.

    The fitted probabilities and design are fixed inputs.  Thus this is a
    conditional measurement-noise reference: only finite cell sampling varies
    between replicates, while the observed covariate structure is retained.
    """
    n_samples, n_terms = design.shape
    if n_samples <= n_terms:
        return np.nan
    denominators = np.asarray(denominators, dtype=np.int64)
    fitted_probabilities = np.asarray(fitted_probabilities, dtype=float)
    simulated = rng.binomial(
        denominators, fitted_probabilities, size=(n_bootstrap, n_samples),
    ) / denominators
    orthonormal_design, _ = np.linalg.qr(design, mode="reduced")
    # Refit the same fixed design because the observed quantity is also a
    # post-covariate-model residual SD.  This removes chance alignment between
    # simulated binomial noise and the covariates, just as in the observed fit.
    residuals = simulated - (simulated @ orthonormal_design) @ orthonormal_design.T
    return float(np.sqrt(np.mean(np.sum(residuals**2, axis=1) / (n_samples - n_terms))))


def fit_parent_fraction_models(
    fractions, adata, *, fraction_column, bootstrap_replicates, bootstrap_seed,
):
    """Fit adult-only, study-specific OLS models and residual diagnostics."""
    candidate_covariates = ["age", "sex", "bmi", "cmv"]
    present = [column for column in candidate_covariates if column in adata.obs]
    metadata_columns = [column for column in present if column != "age"]
    metadata = adata.obs.groupby(["study", "sample"], observed=True)[metadata_columns].first().reset_index()
    # Repeated cell metadata must describe one sample consistently.
    for column in present:
        n_unique = adata.obs.groupby(["study", "sample"], observed=True)[column].nunique(dropna=False)
        if (n_unique > 1).any():
            raise ValueError(f"{column} must be constant within each study × sample")
    data = fractions.merge(metadata, on=["study", "sample"], how="left", validate="one_to_one")
    data["age"] = pd.to_numeric(data["age"], errors="coerce")
    if "cmv" in present:
        data["cmv"] = normalize_cmv_status(data["cmv"])
    data = data.loc[data["age"] >= 20].copy()
    results = []
    residual_diagnostics = []
    unavailable = {"", "not_provided", "unknown", "nan", "none"}
    bootstrap_rng = np.random.default_rng(bootstrap_seed)

    def valid_values(values, column):
        if column in {"age", "bmi"}:
            numeric = pd.to_numeric(values, errors="coerce")
            return numeric.notna() & np.isfinite(numeric)
        normalized = values.astype("string").str.strip().str.lower()
        return normalized.notna() & ~normalized.isin(unavailable)

    for study, subset in data.groupby("study", sort=True, observed=True):
        eligible = subset.copy()
        selected = []
        for column in candidate_covariates:
            if column not in present:
                continue
            valid = valid_values(eligible[column], column)
            coverage = float(valid.mean()) if len(eligible) else np.nan
            if coverage >= 0.95:
                if column in {"sex", "cmv"}:
                    levels = categorical_levels(column, eligible.loc[valid, column])
                    variable = len(levels) >= 2
                else:
                    variable = pd.to_numeric(eligible.loc[valid, column], errors="coerce").nunique() >= 2
                if variable:
                    selected.append(column)

        if not selected:
            continue
        valid_rows = pd.Series(True, index=eligible.index)
        for column in selected:
            valid_rows &= valid_values(eligible[column], column)
        model_data = eligible.loc[valid_rows].copy()
        if model_data.empty:
            continue
        pieces = [pd.Series(1.0, index=model_data.index, name="intercept")]
        term_covariates = {"intercept": "intercept"}
        for column in selected:
            if column in {"age", "bmi"}:
                pieces.append(pd.to_numeric(model_data[column], errors="coerce").rename(column))
                term_covariates[column] = column
            else:
                values = model_data[column].astype("string").str.strip().str.lower()
                levels = categorical_levels(column, values)
                values = pd.Series(pd.Categorical(values, categories=levels), index=model_data.index)
                dummies = pd.get_dummies(values, prefix=column, drop_first=True, dtype=float)
                for term in dummies:
                    pieces.append(dummies[term].rename(term))
                    term_covariates[term] = column
        design = pd.concat(pieces, axis=1).astype(float)
        if len(model_data) <= design.shape[1] or np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
            continue
        design_array = design.to_numpy()
        fit = sm.OLS(model_data[fraction_column].to_numpy(dtype=float), design_array).fit()
        confidence = fit.conf_int()
        mean_fraction = float(model_data[fraction_column].mean())
        unbounded_fitted_probabilities = np.asarray(fit.predict(design_array), dtype=float)
        n_clipped_probabilities = int(
            ((unbounded_fitted_probabilities < 0) | (unbounded_fitted_probabilities > 1)).sum()
        )
        # OLS is retained as the reported fraction model.  Its fitted values
        # provide donor-specific conditional means for the binomial null; clip
        # only because binomial probabilities must lie in [0, 1].  The count is
        # published so any material boundary correction is visible.
        fitted_probabilities = np.clip(unbounded_fitted_probabilities, 0, 1)
        residual_diagnostics.append({
            "study": str(study),
            "residual_sd": float(np.sqrt(fit.mse_resid)),
            "binomial_sampling_sd": _bootstrap_binomial_residual_sd(
                design_array,
                model_data["denominator_cells"].to_numpy(),
                fitted_probabilities,
                rng=bootstrap_rng,
                n_bootstrap=bootstrap_replicates,
            ),
            "n_samples": len(model_data),
            "mean_fraction": mean_fraction,
            "n_clipped_binomial_probabilities": n_clipped_probabilities,
            "fraction_clipped_binomial_probabilities": n_clipped_probabilities / len(model_data),
            "covariates": ";".join(selected),
            "model_terms": ";".join(design.columns.drop("intercept")),
        })
        for index, term in enumerate(design.columns):
            if term == "intercept":
                continue
            covariate = term_covariates[term]
            if covariate in {"sex", "cmv"}:
                levels = categorical_levels(
                    covariate,
                    model_data[covariate].astype("string").str.strip().str.lower(),
                )
                contrast = f"{term.removeprefix(f'{covariate}_')} vs {levels[0]}"
            else:
                contrast = "per decade" if covariate == "age" else "per BMI unit"
            effect_scale = 10 if covariate == "age" else 1
            results.append({
                "study": str(study), "covariate": covariate, "term": term,
                "contrast": contrast, "estimate": fit.params[index] * effect_scale,
                "ci_low": confidence[index, 0] * effect_scale,
                "ci_high": confidence[index, 1] * effect_scale,
                "n_samples": len(model_data),
            })
    coefficient_columns = [
        "study", "covariate", "term", "contrast", "estimate", "ci_low", "ci_high",
        "n_samples",
    ]
    diagnostic_columns = [
        "study", "residual_sd", "binomial_sampling_sd", "n_samples", "mean_fraction",
        "n_clipped_binomial_probabilities", "fraction_clipped_binomial_probabilities",
        "covariates", "model_terms",
    ]
    # Keep the report-facing schema stable even when no study has enough adult
    # samples or covariate variation for an estimable model.
    return (
        pd.DataFrame(results, columns=coefficient_columns),
        pd.DataFrame(residual_diagnostics, columns=diagnostic_columns),
    )


def plot_parent_fraction_coefficients(
    coefficients, residual_diagnostics, *, palette, cell_type, parent,
):
    """Plot adjusted coefficients plus observed and binomial residual variation."""
    if coefficients.empty:
        display(Markdown("_No estimable adult-only study models for this fraction._"))
        return
    covariate_order = [name for name in ["age", "sex", "bmi", "cmv"]
                       if name in coefficients["covariate"].unique()]
    studies = sorted(coefficients["study"].unique())
    residual_label = "residual_sd"
    fraction_to_percentage_points = 100.0
    coefficient_rows = []
    for covariate in covariate_order:
        subset = coefficients.loc[coefficients["covariate"] == covariate]
        for contrast in sorted(subset["contrast"].astype(str).unique()):
            if covariate in {"sex", "cmv"}:
                label = categorical_contrast_label(
                    covariate, contrast.split(" vs ", maxsplit=1)[0],
                    contrast.split(" vs ", maxsplit=1)[1],
                )
            elif covariate == "age":
                label = "AGE (per decade)"
            else:
                label = "BMI (per BMI unit)"
            coefficient_rows.append((covariate, contrast, label))
    y_positions = {
        (covariate, contrast): index
        for index, (covariate, contrast, _) in enumerate(coefficient_rows)
    }
    residual_y = len(coefficient_rows)
    y_positions[residual_label] = residual_y

    # Keep broad confidence intervals from setting the scale. Estimates and
    # residual variation define the visible range; truncated CI ends are marked
    # explicitly at the panel boundary below.
    scale_values = [0.0]
    scale_values.extend(
        pd.to_numeric(coefficients["estimate"], errors="coerce")
        * fraction_to_percentage_points
    )
    for column in ["residual_sd", "binomial_sampling_sd"]:
        if column in residual_diagnostics:
            scale_values.extend(
                pd.to_numeric(residual_diagnostics[column], errors="coerce")
                * fraction_to_percentage_points
            )
    scale_values = np.asarray(scale_values, dtype=float)
    scale_values = scale_values[np.isfinite(scale_values)]
    data_min, data_max = float(scale_values.min()), float(scale_values.max())
    padding = max((data_max - data_min) * 0.08, max(abs(data_min), abs(data_max), 1.0) * 0.02)
    x_limits = (data_min - padding, data_max + padding)

    figure, axis = plt.subplots(figsize=(10, max(4.8, 1.15 * (len(coefficient_rows) + 1))))
    clipped_ci_ends = []
    for study_index, study in enumerate(studies):
        subset = coefficients.loc[coefficients["study"] == study]
        for _, row in subset.iterrows():
            offset = (study_index - (len(studies) - 1) / 2) * 0.055
            y = y_positions[(row["covariate"], row["contrast"])] + offset
            estimate = row["estimate"] * fraction_to_percentage_points
            ci_low = row["ci_low"] * fraction_to_percentage_points
            ci_high = row["ci_high"] * fraction_to_percentage_points
            axis.errorbar(
                estimate, y, xerr=[[estimate - ci_low], [ci_high - estimate]],
                fmt="o", color=palette[study], capsize=3, markersize=5,
            )
            if ci_low < x_limits[0]:
                clipped_ci_ends.append((x_limits[0], y, "<", palette[study]))
            if ci_high > x_limits[1]:
                clipped_ci_ends.append((x_limits[1], y, ">", palette[study]))
    for study_index, study in enumerate(studies):
        subset = residual_diagnostics.loc[residual_diagnostics["study"] == study]
        if subset.empty:
            continue
        row = subset.iloc[0]
        # Interpretation check: assess whether residual variation is too large
        # for the coefficient estimates above to support useful interpretation.
        offset = (study_index - (len(studies) - 1) / 2) * 0.055
        y = residual_y + offset
        axis.plot(
            row["residual_sd"] * fraction_to_percentage_points,
            y, "o", color=palette[study], markersize=6,
        )
        axis.plot(
            row["binomial_sampling_sd"] * fraction_to_percentage_points,
            y, marker="D", linestyle="", markersize=5,
            markerfacecolor="none", markeredgecolor=palette[study], markeredgewidth=1.25,
        )
    axis.axvline(0, color="#333333", linewidth=1, linestyle="--")
    axis.set_xlim(*x_limits)
    for boundary, y, direction, color in clipped_ci_ends:
        axis.plot(
            boundary, y, marker=direction, linestyle="", color=color,
            markersize=5, clip_on=False, zorder=5,
        )
    axis.set_yticks(
        list(range(residual_y + 1)),
        [*(label for _, _, label in coefficient_rows), "RESIDUAL SD"],
    )
    axis.set_xlabel("Fraction change (percentage points; effects with 95% CI)")
    axis.set_ylabel("Covariate / diagnostic")
    axis.invert_yaxis()
    study_legend = axis.legend(
        handles=[plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study)
                 for study in studies],
        title="Study", bbox_to_anchor=(1.02, 1), loc="upper left",
    )
    axis.add_artist(study_legend)
    axis.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="", color="#333333", label="Observed residual SD"),
            plt.Line2D([], [], marker="D", linestyle="", color="#333333", markerfacecolor="none",
                       label="Simulated sampling SD"),
        ],
        title="Residual diagnostic", bbox_to_anchor=(1.02, 0.53), loc="upper left",
    )
    axis.set_title(
        f"Adult-only covariate effects on fraction of {cell_type}\nwithin {parent}"
    )
    figure.tight_layout()
    plt.show()


def directional_cluster_genes(markers, *, cluster, n_genes=8):
    """Return the strongest positive and negative genes for one cluster contrast."""
    subset = markers.loc[markers["group"].astype(str) == str(cluster)].copy()
    positive = subset.loc[subset["logfoldchanges"] > 0].nlargest(n_genes, "logfoldchanges")
    negative = subset.loc[subset["logfoldchanges"] < 0].nsmallest(n_genes, "logfoldchanges")
    result = pd.concat([positive, negative]).copy()
    result["direction"] = result["logfoldchanges"].ge(0).map({
        True: "Higher in cluster", False: "Lower in cluster",
    })
    return result.reset_index(drop=True)


def plot_pca_study_technology_and_introns(adata):
    """Show native PC scores with detailed and collapsed protocol covariates."""
    scores = adata.obsm["X_pca"]
    if scores.shape[1] < 2:
        display(Markdown("_Fewer than two local PCs are available to plot._"))
        return

    pairs = [(0, 1)]
    if scores.shape[1] >= 4:
        pairs.append((2, 3))
    score_data = adata.obs[["study", "technology", "include_intronic"]].copy()
    score_data["study"] = score_data["study"].astype(str)
    score_data["technology"] = score_data["technology"].astype(str)
    score_data["include_intronic"] = score_data["include_intronic"].astype(str)

    def technology_end(value):
        normalized = value.lower().replace("′", "'").replace(" ", "")
        if "3'" in normalized:
            return "3′"
        if "5'" in normalized:
            return "5′"
        return "other"

    intronic_labels = {
        "yes": "Yes",
        "no": "No",
        "not_provided": "Not provided",
    }
    score_data["technology_end"] = score_data["technology"].map(technology_end)
    score_data["include_intronic_label"] = score_data["include_intronic"].map(
        lambda value: intronic_labels.get(value.lower(), value)
    )
    for component in {component for pair in pairs for component in pair}:
        score_data[f"PC {component + 1}"] = scores[:, component]

    studies = sorted(score_data["study"].unique())
    technologies = sorted(score_data["technology"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    intronic_levels = [
        label for label in ["Yes", "No", "Not provided"]
        if label in score_data["include_intronic_label"].unique()
    ]
    intronic_palette = {
        level: {"Yes": "#0072B2", "No": "#D55E00", "Not provided": "#7F7F7F"}[level]
        for level in intronic_levels
    }
    marker_cycle = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "h"]
    technology_markers = {
        technology: marker_cycle[index % len(marker_cycle)]
        for index, technology in enumerate(technologies)
    }
    technology_end_order = [
        end for end in ["3′", "5′", "other"] if end in score_data["technology_end"].unique()
    ]
    technology_end_markers = {"3′": "o", "5′": "s", "other": "X"}
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])

    figure, axes = plt.subplots(2, len(pairs), figsize=(9 * len(pairs), 12), squeeze=False)
    for axis, (x_component, y_component) in zip(axes[0], pairs):
        sns.scatterplot(
            data=score_data,
            x=f"PC {x_component + 1}",
            y=f"PC {y_component + 1}",
            hue="study",
            style="technology",
            palette=palette,
            markers=technology_markers,
            s=max(8, 120_000 / adata.n_obs),
            alpha=0.65,
            edgecolor="none",
            legend=False,
            rasterized=True,
            ax=axis,
        )
        axis.set_title(f"PC {x_component + 1} vs. PC {y_component + 1}: study and technology")
        axis.set_xlabel(f"PC {x_component + 1} ({variance_ratio[x_component]:.1%} variance)")
        axis.set_ylabel(f"PC {y_component + 1} ({variance_ratio[y_component]:.1%} variance)")
    for axis, (x_component, y_component) in zip(axes[1], pairs):
        sns.scatterplot(
            data=score_data,
            x=f"PC {x_component + 1}",
            y=f"PC {y_component + 1}",
            hue="include_intronic_label",
            hue_order=intronic_levels,
            style="technology_end",
            style_order=technology_end_order,
            palette=intronic_palette,
            markers=technology_end_markers,
            s=max(8, 120_000 / adata.n_obs),
            alpha=0.65,
            edgecolor="none",
            legend=False,
            rasterized=True,
            ax=axis,
        )
        axis.set_title(f"PC {x_component + 1} vs. PC {y_component + 1}: intronic reads and 3′/5′")
        axis.set_xlabel(f"PC {x_component + 1} ({variance_ratio[x_component]:.1%} variance)")
        axis.set_ylabel(f"PC {y_component + 1} ({variance_ratio[y_component]:.1%} variance)")

    study_handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study, markersize=6)
        for study in studies
    ]
    technology_handles = [
        plt.Line2D(
            [], [], marker=technology_markers[technology], linestyle="", color="0.35",
            label=technology, markersize=6,
        )
        for technology in technologies
    ]
    intronic_handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=intronic_palette[level], label=level, markersize=6)
        for level in intronic_levels
    ]
    technology_end_handles = [
        plt.Line2D(
            [], [], marker=technology_end_markers[end], linestyle="", color="0.35",
            label=end, markersize=6,
        )
        for end in technology_end_order
    ]
    study_legend = figure.legend(
        handles=study_handles, title="Study", loc="upper left",
        bbox_to_anchor=(0.81, 0.94),
    )
    figure.add_artist(study_legend)
    technology_legend = figure.legend(
        handles=technology_handles, title="Technology", loc="upper left",
        bbox_to_anchor=(0.81, 0.67),
    )
    figure.add_artist(technology_legend)
    intronic_legend = figure.legend(
        handles=intronic_handles, title="Intronic reads used\nin alignment", loc="upper left",
        bbox_to_anchor=(0.81, 0.39),
    )
    figure.add_artist(intronic_legend)
    figure.legend(
        handles=technology_end_handles, title="Technology end", loc="upper left",
        bbox_to_anchor=(0.81, 0.17),
    )
    figure.suptitle("Native local PCA scores by study and protocol covariates", y=0.99)
    figure.tight_layout(rect=(0, 0, 0.8, 0.96))
    display_collapsible_figure(figure, "Show PCA scores by study and protocol")

# %% [markdown]
# <a id="sample-coverage"></a>
# ## Sample Coverage Among Retained Samples
#
# One row represents one study × sample. The compact table retains the metadata
# context for the fraction plots below without repeating separate count plots
# for fields that are usually constant within a study.

# %%
sample_covariates = (
    primitive.obs.groupby(["study", "sample"], observed=True)
    .agg(
        age=("age", "first"),
        sex=("sex", "first"),
        technology=("technology", "first"),
        include_intronic=("include_intronic", "first"),
        disease_status=("disease_status", "first"),
    )
    .reset_index()
)
study_coverage = (
    sample_covariates.groupby("study", observed=True)
    .agg(
        n_samples=("sample", "nunique"),
        age_min=("age", "min"),
        age_max=("age", "max"),
        sex=("sex", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        technology=("technology", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        include_intronic=("include_intronic", lambda values: ", ".join(sorted(values.astype(str).unique()))),
        disease_status=("disease_status", lambda values: ", ".join(sorted(values.astype(str).unique()))),
    )
    .reset_index()
)
display_collapsible_table(study_coverage, "Show study coverage and protocol details")
figure, axis = plt.subplots(figsize=(10, 4))
sns.stripplot(data=sample_covariates, x="study", y="age", hue="sex", dodge=True, ax=axis)
axis.set_title("Age coverage by study and sex")
axis.tick_params(axis="x", rotation=45)
figure.tight_layout()
plt.show()

# %% [markdown]
# <a id="sample-fractions"></a>
# ## Sample-Level Fraction Across Age
#
# Each point is one study × sample. The left panel uses all retained PBMCs as
# denominator; the right panel uses the manually configured AIFI-L1 parent
# derived from every cell's merged L2 call, never from a separate L1 classifier.
# Samples with fewer than the configured denominator-cell cutoff are excluded.
# Facets show descriptive within-study linear fits for all samples and for age
# ≥20. The all-study view places all studies in one panel per denominator, with
# semi-transparent points and a separate lowess curve per study. It is purely
# descriptive: no combined model or pooled age effect is fitted.

# %%
min_fraction_denominator = pipeline["cell_type_analysis"]["min_fraction_denominator_cells"]
fraction_specs = []
for denominator, fraction_column, title, ylabel in [
    ("n_cells_in_sample", "fraction_of_retained_pbmc", "Fraction of all retained PBMCs", "Fraction of retained PBMCs"),
    ("n_cells_in_sample_l1_parent", "fraction_within_aifi_l1_parent", None, None),
]:
    fraction_specs.append((
        sample_cell_type_fractions(
            adata, denominator=denominator, fraction_name=fraction_column,
            min_denominator_cells=min_fraction_denominator,
        ),
        fraction_column,
        title,
        ylabel,
    ))
parent = adata.obs["aifi_l1_parent_for_l2"].iloc[0]
fraction_specs[1] = (
    fraction_specs[1][0], fraction_specs[1][1], f"Fraction within {parent}", f"Fraction of {parent}"
)
all_fraction_samples = [
    len(sample_cell_type_fractions(adata, denominator=denominator, fraction_name=fraction_column))
    for denominator, fraction_column, _, _ in [
        ("n_cells_in_sample", "fraction_of_retained_pbmc", None, None),
        ("n_cells_in_sample_l1_parent", "fraction_within_aifi_l1_parent", None, None),
    ]
]
display_collapsible_table(pd.DataFrame({
    "fraction": [spec[2] for spec in fraction_specs],
    "eligible_samples": [len(spec[0]) for spec in fraction_specs],
    "excluded_below_denominator_cutoff": [
        total - len(spec[0]) for total, spec in zip(all_fraction_samples, fraction_specs)
    ],
    "minimum_denominator_cells": min_fraction_denominator,
}), "Show sample eligibility by fraction denominator")

studies = sorted(pd.concat([spec[0] for spec in fraction_specs])["study"].unique())
palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))

# %%
display(Markdown(
    f"### Adjusted Fraction Effects: {cell_type_name} Within {parent}"
))

# %% [markdown]
# Each study has a separate ordinary least-squares model with one row per
# eligible sample aged **20 years or older**. The response is this cell type's
# fraction of its derived AIFI-L1 parent compartment. Age and BMI are numeric;
# sex and CMV are categorical treatment-coded terms. A covariate is included
# when it is present and observed in at least 95% of adult samples in that
# study. Models use complete cases for the included covariates, and categorical
# contrasts use female as the sex reference and no CMV (negative, when that is
# the source label) as the CMV reference. Age coefficients and intervals are
# scaled to a 10-year change. The bottom row compares each model's
# residual SD with a conditional parametric bootstrap of binomial cell-sampling
# noise. For donor *i*, the reported OLS model supplies its fitted fraction
# \(p_i\), clipped to [0, 1], and its observed parent-cell count is \(n_i\).
# Each replicate draws \(K_i \sim Binomial(n_i, p_i)\), analyzes
# \(K_i / n_i\) with the same covariate design, and records the residual SD.
# This holds the fitted covariate structure fixed and varies only finite
# cell-count sampling; it is therefore directly comparable with the observed
# post-model residual SD. The number and fraction of clipped OLS predictions
# are published in the diagnostics TSV. More bootstrap replicates improve Monte-Carlo
# precision, but do not change the target quantity. The residuals need review
# before interpreting the other terms: if they are too large relative to their
# scale, those effects may not support a reasonable interpretation. In the forest
# plot, “Observed residual SD” is the fitted OLS residual SD; “Simulated sampling
# SD” is the mean residual SD across the conditional binomial bootstrap replicates.
# The response is a fraction on the 0–1 scale. The plot multiplies effects, their
# 95% confidence limits, and residual SDs by 100 for display in percentage points
# (a 0.01 fraction change is 1 percentage point); these are absolute percentage-
# point changes, not relative percent changes.

# %%
display(Markdown("**Model sample restriction:** adults aged ≥20 years only."))
parent_fraction_data = fraction_specs[1][0]
parent_fraction_column = fraction_specs[1][1]
fraction_model_settings = pipeline["cell_type_analysis"]["fraction_model"]
fraction_coefficients, fraction_residual_diagnostics = fit_parent_fraction_models(
    parent_fraction_data, adata, fraction_column=parent_fraction_column,
    bootstrap_replicates=fraction_model_settings["binomial_bootstrap_replicates"],
    bootstrap_seed=fraction_model_settings["binomial_bootstrap_seed"],
)
fraction_residual_diagnostics.insert(0, "cell_type", cell_type_name)
fraction_residual_diagnostics.insert(1, "aifi_l1_parent", parent)
fraction_residual_diagnostics.insert(2, "response", parent_fraction_column)
fraction_residual_diagnostics["binomial_bootstrap_replicates"] = (
    fraction_model_settings["binomial_bootstrap_replicates"]
)
fraction_residual_diagnostics["binomial_bootstrap_seed"] = fraction_model_settings["binomial_bootstrap_seed"]
fraction_diagnostics_path = output_dir / "fraction_model_diagnostics.tsv"
fraction_residual_diagnostics.to_csv(fraction_diagnostics_path, sep="\t", index=False)
display(Markdown("#### Fraction-Model Residual and Sampling Diagnostics"))
display_collapsible_table(
    fraction_residual_diagnostics, "Show full per-study residual diagnostics",
)

age_effects = fraction_coefficients.loc[fraction_coefficients["covariate"] == "age"]
n_positive_age_effects = int((age_effects["estimate"] > 0).sum())
n_negative_age_effects = int((age_effects["estimate"] < 0).sum())
if age_effects.empty:
    age_direction = "No estimable age coefficient"
elif len(age_effects) == 1:
    age_direction = "1 positive estimate" if n_positive_age_effects else "1 negative estimate"
elif n_positive_age_effects and n_negative_age_effects:
    age_direction = f"Mixed: {n_positive_age_effects} positive, {n_negative_age_effects} negative"
else:
    age_direction = (
        f"{n_positive_age_effects} positive estimates"
        if n_positive_age_effects else f"{n_negative_age_effects} negative estimates"
    )
sampling_ratio = (
    fraction_residual_diagnostics["residual_sd"]
    / fraction_residual_diagnostics["binomial_sampling_sd"].replace(0, np.nan)
)
sampling_ratio = sampling_ratio.replace([np.inf, -np.inf], np.nan).dropna()
residual_summary = (
    f"Median {sampling_ratio.median():.1f}× the conditional sampling reference"
    if not sampling_ratio.empty else "Sampling comparison unavailable"
)
model_covariates = sorted({
    covariate
    for terms in fraction_residual_diagnostics["covariates"].dropna()
    for covariate in str(terms).split(";") if covariate
})
display(HTML(
    "<section id=\"fraction-model-evidence\" class=\"report-evidence\">"
    "<h3>Fraction-Model Evidence at a Glance</h3>"
    "<p>These are study-specific adult-sample OLS fits, not a pooled effect estimate. "
    "Use them to assess direction, precision, and the residual diagnostic before interpreting individual coefficients.</p>"
    "<div class=\"report-evidence-grid\">"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Estimable study models</p><p class=\"report-evidence-value\">{len(fraction_residual_diagnostics)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Age-effect direction</p><p class=\"report-evidence-value\">{escape(age_direction)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Residual check</p><p class=\"report-evidence-value\">{escape(residual_summary)}</p></div>"
    f"<div class=\"report-evidence-item\"><p class=\"report-evidence-label\">Covariates represented</p><p class=\"report-evidence-value\">{escape(', '.join(model_covariates) or 'None')}</p></div>"
    "</div></section>"
))
plot_parent_fraction_coefficients(
    fraction_coefficients, fraction_residual_diagnostics, palette=palette,
    cell_type=cell_type_name, parent=parent,
)
display(Markdown(
    "_How to read this plot: each colour is a study and intervals are its 95% OLS confidence "
    "intervals. Fraction effects and both residual SD markers are shown in percentage points. "
    "The filled dot on the residual row is the observed residual SD; the open diamond is the "
    "simulated sampling SD described above. Compare studies for consistency and precision, "
    "rather than treating their collection as a pooled estimate. Arrowheads mark 95% CIs "
    "extending beyond the displayed range._"
))

display(Markdown(
    "<a id=\"additional-fraction-trends\"></a>\n### Additional Descriptive Fraction Trends"
))
display(Markdown(
    "These age-trend views are secondary to the adult-only adjusted estimates above. "
    "The within-L1 within-study panel opens on load; the all-retained-PBMC panel starts "
    "closed. Both show all-age and adult-only linear fits. "
    "The all-study panels show one descriptive lowess curve per study."
))
for index, (fractions, fraction_column, title, ylabel) in enumerate(fraction_specs):
    plot_study_linear_trends(
        fractions, fraction_column=fraction_column, title=title, ylabel=ylabel, palette=palette,
        show_fit_legend=index == 0,
        open_by_default=index != 0,
    )
plot_study_fraction_curves(
    [(fractions, fraction_column, ylabel) for fractions, fraction_column, _, ylabel in fraction_specs],
    palette=palette,
)

# %% [markdown]
# <a id="local-embedding"></a>
# ## Local Embedding and Clusters
#
# This is a type-specific exploratory embedding. PCA is calculated from the
# local highly variable genes after normalization, log transformation, and
# scaling; configured V(D)J genes are excluded. When more than one study is
# represented, Harmony corrects these PCs by study. The configured local
# neighbor graph, UMAP, and Leiden clusters are calculated from the
# Harmony-adjusted PCs; native PCs are retained below for diagnostic plots.
# Integration can also change biological structure, so study colouring and
# sample-level cluster composition should be reviewed before interpreting a
# cluster as a biological state.

# %%
if report["status"] != "complete":
    display(Markdown("_Too few cells for a type-specific embedding and clustering._"))
else:
    plot_umap(adata, color="cluster", title="Local UMAP: cluster IDs", label_clusters=True)
    colors = ["study", "age"]
    if "aifi_l3_majority" in adata.obs:
        colors.append("aifi_l3_majority")
    for color in colors:
        plot_umap(adata, color=color, marker_area_scale=4.0)

# %% [markdown]
# <a id="cluster-markers"></a>
# ### Cluster Markers
#
# Genes are ranked by a Wilcoxon comparison of each local cluster against all
# other local clusters. Signed log fold changes therefore show whether a gene
# is higher or lower in that cluster; these are descriptive contrasts.

# %%
if report["status"] == "complete":
    markers = pd.read_csv(output_dir / "markers.tsv", sep="\t")
    clusters = sorted(markers["group"].astype(str).unique(), key=natural_sort_key)
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5.5 * n_columns, 4.4 * n_rows), squeeze=False)
    for axis, cluster in zip(axes.flat, clusters):
        subset = directional_cluster_genes(markers, cluster=cluster)
        sns.barplot(
            data=subset, x="logfoldchanges", y="names", hue="direction", dodge=False,
            order=subset["names"].astype(str).tolist(),
            palette={"Higher in cluster": "#c44e52", "Lower in cluster": "#4c72b0"}, ax=axis,
        )
        axis.axvline(0, color="black", linewidth=0.8)
        axis.set_title(f"Cluster {cluster}: differential genes")
        axis.set_xlabel("Log fold change versus all other clusters")
        axis.set_ylabel("")
        axis.legend(title="Direction", fontsize=8, title_fontsize=8, loc="best")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    figure.tight_layout()
    display_collapsible_figure(
        figure, f"Show marker genes across {len(clusters)} local clusters", open_by_default=True,
    )

# %% [markdown]
# <a id="cluster-composition"></a>
# ### Cluster Composition
#
# The plots use one study × sample fraction per point, so they describe
# sample-level cluster composition rather than treating cells as independent
# observations. Curves are minimized on load.

# %%
if report["status"] == "complete":
    all_cluster_fractions = sample_cluster_fractions(adata)
    cluster_fractions = sample_cluster_fractions(
        adata, min_denominator_cells=min_fraction_denominator,
    )
    excluded_cluster_samples = (
        all_cluster_fractions[["study", "sample"]].drop_duplicates().shape[0]
        - cluster_fractions[["study", "sample"]].drop_duplicates().shape[0]
    )
    display(Markdown(
        f"_{excluded_cluster_samples} sample(s) excluded from cluster fractions because the "
        f"split cell-type denominator has fewer than {min_fraction_denominator} cells._"
    ))
    studies = sorted(cluster_fractions["study"].unique())
    palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
    plot_cluster_fractions(cluster_fractions, palette=palette)
    display(Markdown(
        "_Each point is a sample-level cluster fraction. Differences among coloured study series "
        "are descriptive and should be read alongside the coverage and protocol context above._"
    ))

# %% [markdown]
# <a id="pca-diagnostics"></a>
# ### PCA Diagnostics
#
# These are the native, local PCA scores before Harmony adjustment. Colour uses
# the report's standard study palette and marker shape indicates the recorded
# single-cell technology in the top row. The lower row collapses technology to
# 3′ or 5′ marker shapes and maps intronic read inclusion to colour. The
# Harmony-adjusted representation is used only for the local neighbour graph,
# UMAP, and clustering.

# %%
if report["status"] == "complete":
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
    figure, axis = plt.subplots(figsize=(9, 4.5))
    components = np.arange(1, len(variance_ratio) + 1)
    axis.bar(components, variance_ratio * 100, color="#4c72b0")
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Variance explained (%)")
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    plt.show()

# %%
if report["status"] == "complete":
    plot_pca_study_technology_and_introns(adata)

# %%
if report["status"] == "complete":
    n_components_to_plot = min(10, adata.obsm["X_pca"].shape[1])
    sc.pl.pca_loadings(
        adata, components=list(range(1, n_components_to_plot + 1)), n_points=12, show=False,
    )
    figure = plt.gcf()
    n_columns = min(4, n_components_to_plot)
    n_rows = (n_components_to_plot + n_columns - 1) // n_columns
    figure.set_size_inches(14, 3.0 * n_rows)
    figure.subplots_adjust(left=0.06, right=0.98, bottom=0.08, top=0.94, hspace=0.75, wspace=0.35)
    display_collapsible_figure(figure, f"Show gene loadings for the first {n_components_to_plot} PCs")

# %%
if report["status"] == "complete":
    pc_age = pd.read_csv(output_dir / "pc_age_correlations.tsv", sep="\t")
    figure, axis = plt.subplots(figsize=(9, 4.5))
    sns.barplot(data=pc_age, x="pc", y="spearman_r", color="#4c72b0", ax=axis)
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_ylabel("Spearman correlation with age (cells; descriptive)")
    axis.tick_params(axis="x", rotation=45)
    figure.tight_layout()
    display_collapsible_figure(figure, "Show PC–age correlations")

# %%
if de_dir:
    de_dir = Path(de_dir)
    indexed_result_paths = {
        (de_dir / relative_path).resolve()
        for model in de_run_metadata.get("models", [])
        if model.get("status") == "complete"
        for relative_path in model.get("results_by_covariate", {}).values()
    }
    if "models" in de_run_metadata:
        per_study_files = sorted(
            path for path in (de_dir / "per_study").rglob("*.csv")
            if path.resolve() in indexed_result_paths
        )
        merged_file = de_dir / "merged.csv"
        if merged_file.resolve() not in indexed_result_paths:
            merged_file = None
        combined_files = sorted(
            path for path in (de_dir / "combined").glob("*.csv")
            if path.resolve() in indexed_result_paths
        )
    else:
        per_study_files = sorted((de_dir / "per_study").rglob("*.csv"))
        merged_file = de_dir / "merged.csv"
        combined_files = sorted((de_dir / "combined").glob("*.csv"))
    per_study_results = []
    for result_path in per_study_files:
        result = pd.read_csv(result_path)
        if "covariate" not in result:
            result["covariate"] = "age"
        if "contrast" not in result:
            result["contrast"] = "per year"
        if "study_sample_counts" not in result:
            result["study_sample_counts"] = "{}"
        per_study_results.append(result)
    has_per_study = bool(per_study_results)
    combined_results_by_covariate = {}
    age_trajectory_results = None
    if merged_file is not None and merged_file.is_file():
        result = pd.read_csv(merged_file)
        if "covariate" not in result:
            result["covariate"] = "age"
        if "contrast" not in result:
            result["contrast"] = "per year"
        combined_results_by_covariate["age"] = result
    for result_path in combined_files:
        result = pd.read_csv(result_path)
        if result_path.stem == "age_bin":
            age_trajectory_results = result
            continue
        if "covariate" not in result:
            result["covariate"] = result_path.stem
        if "contrast" not in result:
            result["contrast"] = "per BMI unit" if result_path.stem == "bmi" else result_path.stem
        combined_results_by_covariate[result_path.stem] = result
    has_combined = bool(combined_results_by_covariate)
    if has_per_study or has_combined or de_run_metadata:
        display(Markdown("<a id=\"differential-expression\"></a>\n## Differential Gene Expression"))
        display(Markdown(
            "Pseudobulk models aggregate counts at the study × sample × cell-type level. "
            "They provide the sample-level complement to the descriptive single-cell views above."
        ))
        display(Markdown(
            "**Gene universes:** per-study fits use genes observed in that study; each combined "
            "fit uses the intersection across its included studies, excluding outer-join zeros for "
            "study-absent genes. The cross-cell-type expression atlas uses the all-study intersection."
        ))
        if de_run_metadata.get("analysis_mode") == "test_only":
            display(Markdown(
                "> **TEST OUTPUT — NOT FOR BIOLOGICAL INTERPRETATION.** "
                + de_run_metadata["interpretation_warning"]
            ))
        if not has_per_study and not has_combined:
            display(Markdown(
                "_No estimable PyDESeq2 result files were produced for this cell type._"
            ))
    age_trajectory_fit = next(
        (model for model in de_run_metadata.get("models", [])
         if model.get("purpose") == "age_trajectory"),
        {},
    )
    age_trajectory_info = de_run_metadata.get("age_trajectory", {})
    if has_per_study or has_combined or age_trajectory_fit or age_trajectory_results is not None:
        de_alpha = pipeline["differential_expression"]["alpha"]
        per_study_by_covariate = {}
        for result in per_study_results:
            for covariate, subset in result.groupby("covariate", observed=True):
                per_study_by_covariate.setdefault(str(covariate), []).append(subset)
        age_per_study = (
            pd.concat(per_study_by_covariate.get("age", []), ignore_index=True)
            if per_study_by_covariate.get("age") else pd.DataFrame()
        )
        age_study_results = (
            [group for _, group in age_per_study.groupby("study", observed=True)]
            if not age_per_study.empty else []
        )
        has_age_per_study = bool(age_study_results)
        if has_age_per_study:
            gene_universes = [set(result["gene"].astype(str)) for result in age_study_results]
            common_per_study_genes = set.intersection(*gene_universes)
            per_study_significant = age_per_study.loc[
                pd.to_numeric(age_per_study["padj"], errors="coerce") < de_alpha
            ].copy()
            per_study_hits = per_study_significant.loc[
                per_study_significant["gene"].astype(str).isin(common_per_study_genes)
            ].copy()
            recurrence = (
                per_study_hits.groupby("gene", as_index=False)
                .agg(
                    studies_associated=("study", "nunique"),
                    studies=("study", lambda values: ", ".join(sorted(set(values)))),
                    maximum_absolute_log2_fold_change=(
                        "log2FoldChange", signed_value_at_largest_absolute_magnitude
                    ),
                )
                .assign(_absolute_effect=lambda frame: frame[
                    "maximum_absolute_log2_fold_change"
                ].abs())
                .sort_values(
                    ["studies_associated", "_absolute_effect", "gene"],
                    ascending=[False, False, True],
                )
                .drop(columns="_absolute_effect")
            )
        else:
            recurrence = pd.DataFrame()

        diagnostic_rows = []
        for covariate, model_results in per_study_by_covariate.items():
            covariate_results = pd.concat(model_results, ignore_index=True)
            study_gene_sets = [
                set(group["gene"].astype(str))
                for _, group in covariate_results.groupby("study", observed=True)
            ]
            if not study_gene_sets:
                continue
            common_genes = set.intersection(*study_gene_sets)
            for (study, contrast), study_results in covariate_results.groupby(
                ["study", "contrast"], sort=True, observed=True
            ):
                significant = study_results.loc[
                    pd.to_numeric(study_results["padj"], errors="coerce") < de_alpha
                ]
                in_intersection = significant.loc[
                    significant["gene"].astype(str).isin(common_genes), "gene"
                ].nunique()
                total_significant = int(significant["gene"].nunique())
                outside_intersection = total_significant - int(in_intersection)
                diagnostic_rows.append({
                    "covariate": covariate,
                    "contrast": contrast,
                    "study": study,
                    "design": (
                        study_results["design"].iloc[0]
                        if "design" in study_results else "—"
                    ),
                    "complete-case samples": (
                        int(study_results["n_samples_model"].iloc[0])
                        if "n_samples_model" in study_results else "—"
                    ),
                    "excluded for missing design values": (
                        int(study_results["n_samples_excluded_missing_design_covariates"].iloc[0])
                        if "n_samples_excluded_missing_design_covariates" in study_results else "—"
                    ),
                    "significant genes": total_significant,
                    "significant genes in study intersection": int(in_intersection),
                    "significant genes outside intersection": outside_intersection,
                    "outside intersection (%)": (
                        100 * outside_intersection / total_significant if total_significant else 0.0
                    ),
                })
        intersection_diagnostics = pd.DataFrame(diagnostic_rows)

        summary_rows = []
        if has_age_per_study:
            summary_rows.append(
                {
                    "model": "Per-study age significant-gene union",
                    "FDR-significant age-associated genes": int(
                        per_study_significant["gene"].nunique()
                    ),
                }
            )
            summary_rows.append(
                {
                    "model": "Per-study age significant genes in all-study intersection",
                    "FDR-significant age-associated genes": int(per_study_hits["gene"].nunique()),
                }
            )
        combined_fit_records = [
            model for model in de_run_metadata.get("models", [])
            if model.get("model") in {"merged", "combined"}
        ]
        age_diagnostics_path = de_dir / "age_model_diagnostics.csv"
        age_combined_results = combined_results_by_covariate.get("age")
        if age_diagnostics_path.is_file() and age_combined_results is not None:
            age_diagnostics = pd.read_csv(age_diagnostics_path)
            expected_columns = {
                "study", "age", "sex", "counts_in_age_gene_intersection",
                "genes_in_age_intersection",
            }
            if expected_columns.issubset(age_diagnostics):
                display(Markdown(
                    "<a id=\"shared-age-model-diagnostics\"></a>\n### Shared Age-Model Diagnostics"
                ))
                display(Markdown(
                    "This single figure checks study-site/age support, the relation between age and "
                    "the raw pseudobulk counts used in the shared age gene intersection, and "
                    "the unadjusted age-model p-value calibration. Q-Q points are the raw "
                    "p-values (`pvalue`, before multiple-testing correction) from the combined "
                    "age model fitted across all eligible studies for this cell type, adjusted "
                    "for study site, sex, age, and log10(total_counts), and using the model's shared "
                    "gene intersection. "
                    "The dashed y=x line is the null reference; "
                    "the axes retain independent scales so departures remain readable. This is "
                    "descriptive: it does not diagnose a specific gene or replace model checks."
                ))
                figure, axes = plt.subplots(1, 3, figsize=(18, 4.8))
                support_column = "study_site" if "study_site" in age_diagnostics else "study"
                sns.stripplot(
                    data=age_diagnostics, x=support_column, y="age", hue="sex", dodge=True,
                    jitter=0.18, alpha=0.8, ax=axes[0],
                )
                axes[0].set_xlabel("Study site" if support_column == "study_site" else "Study")
                axes[0].set_ylabel("Age (years)")
                axes[0].set_title("Age support by study site")
                axes[0].tick_params(axis="x", rotation=75, labelsize=7)
                axes[0].legend(title="Sex", fontsize=8, title_fontsize=8)

                plot_samples = age_diagnostics.loc[
                    age_diagnostics["counts_in_age_gene_intersection"] > 0
                ].copy()
                plot_samples["log10_counts"] = np.log10(
                    plot_samples["counts_in_age_gene_intersection"]
                )
                sns.scatterplot(
                    data=plot_samples, x="age", y="log10_counts", hue=support_column, style="sex",
                    s=42, alpha=0.85, ax=axes[1],
                )
                axes[1].set_xlabel("Age (years)")
                axes[1].set_ylabel("log10 raw counts in age gene intersection")
                axes[1].set_title("Age and pseudobulk depth")
                axes[1].legend(fontsize=7, title_fontsize=8)

                pvalue_source = (
                    age_combined_results["pvalue"]
                    if "pvalue" in age_combined_results else pd.Series(dtype=float)
                )
                pvalues = pd.to_numeric(pvalue_source, errors="coerce")
                pvalues = pvalues.loc[(pvalues > 0) & (pvalues <= 1)].sort_values().to_numpy()
                if len(pvalues):
                    expected = (np.arange(1, len(pvalues) + 1) - 0.5) / len(pvalues)
                    axes[2].scatter(
                        -np.log10(expected), -np.log10(pvalues), s=9, alpha=0.65,
                        color="#4c72b0", linewidths=0,
                    )
                    axes[2].axline(
                        (0, 0), slope=1, color="black", linestyle="--",
                        linewidth=1.5, zorder=4, label="Null (y = x)",
                    )
                    axes[2].set_xlabel("Expected −log10(p)")
                    axes[2].set_ylabel("Observed −log10(p)")
                    axes[2].set_title("Age-model p-value QQ plot")
                    axes[2].legend(fontsize=8)
                else:
                    axes[2].text(
                        0.5, 0.5, "No finite unadjusted p-values", ha="center", va="center",
                        transform=axes[2].transAxes,
                    )
                    axes[2].set_axis_off()
                figure.tight_layout()
                plt.show()

        display(Markdown("### Covariate Effects"))
        display(Markdown(
            f"Each combined fit estimates one covariate effect from eligible studies using "
            f"the complete cases in its displayed design. **FDR threshold:** {de_alpha:g}. "
            "Model details and sample counts are listed below each volcano plot."
        ))
        if not combined_results_by_covariate:
            display(Markdown("_No combined covariate results are available for volcano plots._"))
        covariate_order = {
            "age": 0,
            "sex": 1,
            "bmi": 2,
            "cmv": 3,
            "log10_total_counts": 4,
        }
        ordered_combined_covariates = sorted(
            combined_results_by_covariate.items(),
            key=lambda item: (covariate_order.get(item[0], 99), item[0]),
        )
        for covariate, combined_results in ordered_combined_covariates:
            covariate_label = {
                "age": "Age",
                "sex": "Sex",
                "bmi": "BMI",
                "cmv": "CMV",
                "log10_total_counts": "log10(total_counts)",
            }.get(covariate, covariate.title())
            covariate_anchor_suffix = re.sub(
                r"[^a-z0-9]+", "-", covariate_label.lower()
            ).strip("-")
            display(Markdown(
                f'<a id="combined-{covariate_anchor_suffix}"></a>\n#### {covariate_label}'
            ))
            fit_record = next(
                (model for model in combined_fit_records
                 if covariate in model.get("results_by_covariate", {})),
                {},
            )
            for contrast, contrast_results in combined_results.groupby("contrast", sort=True):
                volcano = contrast_results.copy()
                volcano["padj"] = pd.to_numeric(volcano["padj"], errors="coerce")
                volcano["log2FoldChange"] = pd.to_numeric(
                    volcano["log2FoldChange"], errors="coerce"
                )
                volcano = volcano.dropna(subset=["padj", "log2FoldChange"]).copy()
                if volcano.empty:
                    display(Markdown(
                        f"_{covariate} ({contrast}) has no finite adjusted p-values to plot._"
                    ))
                    continue
                volcano["minus_log10_padj"] = -np.log10(
                    volcano["padj"].clip(lower=np.finfo(float).tiny)
                )
                study_results = per_study_by_covariate.get(covariate, [])
                matching_study_results = pd.DataFrame()
                if study_results:
                    matching_study_results = pd.concat(study_results, ignore_index=True)
                    matching_study_results = matching_study_results.loc[
                        matching_study_results["contrast"].astype(str) == str(contrast)
                    ].copy()
                    matching_study_results["gene"] = matching_study_results[
                        "gene"
                    ].astype(str)
                    matching_study_results["padj"] = pd.to_numeric(
                        matching_study_results["padj"], errors="coerce"
                    )
                    study_hits = (
                        matching_study_results.loc[
                            matching_study_results["padj"] < de_alpha
                        ]
                        .groupby("gene")["study"]
                        .nunique()
                    )
                    volcano["studies_associated"] = (
                        volcano["gene"].astype(str).map(study_hits).fillna(0).astype(int)
                    )
                else:
                    volcano["studies_associated"] = 0
                studies_used = (
                    str(volcano["studies_included"].iloc[0]).split(";")
                    if "studies_included" in volcano else fit_record.get("studies", [])
                )
                sample_counts = (
                    json.loads(str(volcano["study_sample_counts"].iloc[0]))
                    if "study_sample_counts" in volcano
                    else fit_record.get("study_sample_counts", {})
                )
                study_summary = ", ".join(
                    f"{study} (n={sample_counts.get(study, 0)})" for study in studies_used if study
                )
                if not studies_used:
                    studies_used = sorted(sample_counts)
                if not study_summary and sample_counts:
                    study_summary = ", ".join(
                        f"{study} (n={count})" for study, count in sample_counts.items()
                    )
                nonintersection_age_markers = pd.DataFrame()
                if covariate == "age" and not matching_study_results.empty:
                    marker_results = matching_study_results.copy()
                    if studies_used:
                        marker_results = marker_results.loc[
                            marker_results["study"].astype(str).isin(studies_used)
                        ].copy()
                    marker_results["log2FoldChange"] = pd.to_numeric(
                        marker_results["log2FoldChange"], errors="coerce"
                    )
                    marker_results = marker_results.dropna(
                        subset=["gene", "padj", "log2FoldChange"]
                    )
                    study_gene_sets = [
                        set(group["gene"].astype(str))
                        for _, group in marker_results.groupby("study", observed=True)
                    ]
                    shared_age_genes = (
                        set.intersection(*study_gene_sets) if study_gene_sets else set()
                    )
                    nonintersection_age_markers = marker_results.loc[
                        (marker_results["padj"] < de_alpha)
                        & ~marker_results["gene"].astype(str).isin(shared_age_genes)
                    ].copy()
                    if not nonintersection_age_markers.empty:
                        nonintersection_age_markers["minus_log10_padj"] = -np.log10(
                            nonintersection_age_markers["padj"].clip(
                                lower=np.finfo(float).tiny
                            )
                        )
                        nonintersection_age_markers = (
                            nonintersection_age_markers.sort_values(
                                ["study", "padj", "gene"], kind="stable"
                            )
                            .groupby("study", as_index=False, sort=True)
                            .head(1)
                        )
                max_associated_studies = max(
                    int(volcano["studies_associated"].max()), 1
                )
                # Leave room for the study-count legend while keeping these panels
                # comparable in plot-area width to the baseMean-colored panels.
                figure, axis = plt.subplots(figsize=(12.6, 6.5))
                sns.scatterplot(
                    data=volcano, x="log2FoldChange", y="minus_log10_padj",
                    hue="studies_associated", palette="viridis",
                    hue_norm=(0, max_associated_studies), legend="brief",
                    s=16, alpha=0.7, linewidth=0, ax=axis,
                )
                study_count_legend = axis.legend(
                    title="Per-study FDR-significant\nassociations", fontsize=8,
                    bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
                )
                if not nonintersection_age_markers.empty:
                    study_colors = dict(zip(
                        sorted(nonintersection_age_markers["study"].astype(str).unique()),
                        sns.color_palette("tab10", n_colors=nonintersection_age_markers["study"].nunique()),
                    ))
                    for _, marker in nonintersection_age_markers.iterrows():
                        study = str(marker["study"])
                        axis.scatter(
                            marker["log2FoldChange"], marker["minus_log10_padj"],
                            marker="D", s=72, facecolors="none",
                            edgecolors=study_colors[study], linewidths=1.8,
                            label=study, zorder=4,
                        )
                        axis.annotate(
                            f"{marker['gene']} ({study})",
                            (marker["log2FoldChange"], marker["minus_log10_padj"]),
                            xytext=(4, -10), textcoords="offset points", fontsize=7,
                            color=study_colors[study],
                        )
                    axis.legend(
                        title="Top FDR-significant\nnon-shared age gene",
                        fontsize=8, bbox_to_anchor=(1.02, 0.46),
                        loc="upper left", borderaxespad=0,
                    )
                    axis.add_artist(study_count_legend)
                axis.axhline(-np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1)
                axis.axvline(0, color="#555555", linewidth=0.8)
                if covariate == "age":
                    effect_label = "log2 fold change per year"
                elif covariate == "bmi":
                    effect_label = "log2 fold change per BMI unit"
                elif covariate == "log10_total_counts":
                    effect_label = "log2 fold change per log10(total_counts) unit"
                else:
                    effect_label = f"log2 fold change ({contrast})"
                axis.set_xlabel(f"{covariate_label} effect ({effect_label})")
                axis.set_ylabel("−log10(adjusted p-value)")
                if len(studies_used) > 1:
                    title = f"Combined site-adjusted {covariate_label} association: {contrast}"
                else:
                    study_name = studies_used[0] if studies_used else "single study"
                    title = f"{covariate_label} association in {study_name}: {contrast}"
                axis.set_title(title)
                finite_label_rows = volcano.loc[
                    np.isfinite(volcano["padj"])
                    & np.isfinite(volcano["log2FoldChange"])
                ].copy()
                top_significance = finite_label_rows.sort_values(
                    ["padj", "gene"], kind="stable"
                ).head(20)
                top_effect = finite_label_rows.assign(
                    _absolute_log2_fold_change=(
                        finite_label_rows["log2FoldChange"].abs()
                    )
                ).sort_values(
                    ["_absolute_log2_fold_change", "padj", "gene"],
                    ascending=[False, True, True],
                    kind="stable",
                ).head(20)
                labels = pd.concat([top_significance, top_effect]).drop_duplicates(
                    "gene", keep="first"
                )
                for _, row in labels.iterrows():
                    axis.annotate(
                        str(row["gene"]),
                        (row["log2FoldChange"], row["minus_log10_padj"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7,
                    )
                figure.tight_layout(rect=(0, 0, 0.68 if not nonintersection_age_markers.empty else 0.78, 1))
                plt.show()
                if not nonintersection_age_markers.empty:
                    display(Markdown(
                        "Outlined diamonds use each study's own age-model estimate and adjusted "
                        "p-value. They select the smallest adjusted p-value among FDR-significant "
                        "genes outside the shared per-study gene universe."
                    ))
                    display(nonintersection_age_markers[[
                        "study", "gene", "log2FoldChange", "padj",
                    ]].rename(columns={
                        "log2FoldChange": "per-study log2 fold change per year",
                        "padj": "per-study adjusted p-value",
                    }).sort_values("study", kind="stable"))
                if covariate == "log10_total_counts" and "baseMean" in volcano.columns:
                    base_mean = pd.to_numeric(volcano["baseMean"], errors="coerce")
                    finite = np.isfinite(base_mean) & (base_mean > 0)
                    if finite.any():
                        log_base_mean = np.log10(base_mean.loc[finite])
                        base_mean_volcano = volcano.loc[finite]
                        low = float(log_base_mean.min())
                        high = float(log_base_mean.max())
                        norm = Normalize(vmin=low, vmax=high if high > low else low + 0.01)
                        figure, axis = plt.subplots(figsize=(10.5, 6.5))
                        points = axis.scatter(
                            base_mean_volcano["log2FoldChange"],
                            base_mean_volcano["minus_log10_padj"],
                            c=log_base_mean,
                            cmap="viridis",
                            norm=norm,
                            s=16,
                            alpha=0.7,
                            linewidths=0,
                        )
                        colorbar = figure.colorbar(points, ax=axis, pad=0.02)
                        colorbar.set_label("log10(combined-fit baseMean)")
                        axis.axhline(
                            -np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1
                        )
                        axis.axvline(0, color="#555555", linewidth=0.8)
                        axis.set_xlabel(f"{covariate_label} effect ({effect_label})")
                        axis.set_ylabel("−log10(adjusted p-value)")
                        axis.set_title(f"{title} (colored by log10(baseMean))")
                        for _, row in labels.iterrows():
                            axis.annotate(
                                str(row["gene"]),
                                (row["log2FoldChange"], row["minus_log10_padj"]),
                                xytext=(3, 3), textcoords="offset points", fontsize=7,
                            )
                        figure.tight_layout()
                        plt.show()
                fit_design = fit_record.get(
                    "design",
                    str(volcano["design"].iloc[0]) if "design" in volcano else "—",
                )
                study_sites = (
                    str(volcano["study_sites_included"].iloc[0]).split(";")
                    if "study_sites_included" in volcano
                    else fit_record.get("study_sites", [])
                )
                details = pd.DataFrame(
                    {
                        "Value": [
                            contrast,
                            fit_design,
                            f"{de_alpha:g}",
                            int((volcano["padj"] < de_alpha).sum()),
                            study_summary or "unavailable in this legacy result set",
                            ", ".join(site for site in study_sites if site) or "—",
                            sum(sample_counts.values()) if sample_counts else "unavailable",
                            fit_record.get(
                                "n_samples_excluded_missing_design_covariates", "—"
                            ),
                        ]
                    },
                    index=[
                        "Contrast",
                        "Model design",
                        "FDR threshold",
                        "FDR-significant genes",
                        "Studies included (samples per study)",
                        "Study sites included",
                        "Total samples used",
                        "Samples excluded for missing design values",
                    ],
                )
                display(details)
        display(Markdown(
            "_DE interpretation: the recurrence panels emphasize agreement among available per-study "
            "models; each combined volcano summarizes eligible samples with study-site, age, sex, "
            "and log-total-count adjustment. "
            "In the primary volcano, point color gives the number of available per-study fits "
            "for the same covariate and contrast with FDR-significant association for that gene "
            "(using the displayed FDR threshold); the plotted effect and adjusted p-value come "
            "from the combined fit. The log10(total_counts) volcano also has a companion colored "
            "by combined-fit log10(baseMean). "
            "For age, outlined diamonds mark the smallest per-study adjusted p-value among "
            "FDR-significant genes outside the shared per-study gene universe; their position instead uses that study's "
            "own effect and adjusted p-value, so they are a visibility check rather than combined-fit points. "
            "The site term is omitted when only one site remains estimable. Combined fits do not "
            "by themselves establish replication across studies._"
        ))
        per_study_covariate_summary = []
        for result in per_study_results:
            for (covariate, contrast), subset in result.groupby(
                ["covariate", "contrast"], observed=True
            ):
                sample_counts = json.loads(str(subset["study_sample_counts"].iloc[0]))
                study = str(subset["study"].iloc[0])
                fit_record = next(
                    (model for model in de_run_metadata.get("models", [])
                     if model.get("model") == "per_study"
                     and model.get("study") == study
                     and covariate in model.get("results_by_covariate", {})),
                    {},
                )
                per_study_covariate_summary.append({
                    "study": study,
                    "covariate": covariate,
                    "contrast": contrast,
                    "design": fit_record.get("design", "—"),
                    "complete-case samples used": (
                        sum(sample_counts.values()) if sample_counts else "—"
                    ),
                    "excluded for missing design values": fit_record.get(
                        "n_samples_excluded_missing_design_covariates", "—"
                    ),
                    "genes tested": subset["gene"].nunique(),
                    "FDR-significant genes": int(
                        (pd.to_numeric(subset["padj"], errors="coerce") < de_alpha).sum()
                    ),
                })
        if per_study_covariate_summary:
            display(Markdown(
                "<a id=\"per-study-covariate-models\"></a>\n### Per-Study Covariate Models"
            ))
            display(Markdown(
                "Each study is fitted once with its maximal available design. Every coefficient for "
                "that study comes from the same complete-case sample set and formula."
            ))
            per_study_summary = pd.DataFrame(per_study_covariate_summary)
            per_study_summary["_covariate_order"] = per_study_summary["covariate"].map({
                "age": 0, "sex": 1, "bmi": 2, "cmv": 3,
                "log10_total_counts": 4,
            }).fillna(99)
            per_study_summary = (
                per_study_summary
                .sort_values(
                    ["study", "_covariate_order", "contrast"], kind="stable",
                )
                .drop(columns="_covariate_order")
            )
            display_collapsible_table(
                per_study_summary,
                "Show per-study fits, sorted by study then coefficient",
            )

        display(Markdown(
            "Per-study fits use one maximal available design and also adjust for site when "
            "multiple sites contribute. Missing values exclude a sample only from models whose "
            "design includes that variable. Categorical effects use female as the sex reference "
            "and no CMV (negative where that is the source label) as the CMV reference."
        ))
        if summary_rows:
            display_collapsible_table(
                pd.DataFrame(summary_rows), "Show per-study age gene-count summary",
            )

        if has_age_per_study:
            figure, axes = plt.subplots(1, 2, figsize=(14.5, 5.8))
            if recurrence.empty:
                empty_message = (
                    "No genes were tested in every per-study model"
                    if not common_per_study_genes
                    else "No common-universe genes pass the FDR threshold"
                )
                axes[0].text(0.5, 0.5, empty_message,
                             transform=axes[0].transAxes, ha="center", va="center")
                axes[1].text(0.5, 0.5, empty_message,
                             transform=axes[1].transAxes, ha="center", va="center")
            else:
                top_genes = recurrence.head(50).sort_values("studies_associated")
                sns.barplot(
                    data=top_genes, x="studies_associated", y="gene", color="#4c72b0",
                    ax=axes[0],
                )
                axes[0].set_xlabel("Studies with FDR-significant age association")
                axes[0].set_ylabel("")
                axes[0].set_title("Genes recurring across studies (top 50)")
                recurrence_counts = recurrence["studies_associated"].value_counts().sort_index()
                sns.barplot(
                    x=recurrence_counts.index.astype(str), y=recurrence_counts.values,
                    color="#55a868", ax=axes[1],
                )
                axes[1].set_xlabel("Number of studies")
                axes[1].set_ylabel("Unique associated genes")
                axes[1].set_title(
                    f"Common-universe union: {recurrence['gene'].nunique():,} genes"
                )
            figure.tight_layout()
            plt.show()
            if not recurrence.empty:
                display(Markdown(
                    "The effect shown for each gene is the **signed** age log2 fold change "
                    "with the largest absolute magnitude among its per-study estimates."
                ))
                display_collapsible_table(
                    recurrence.head(50),
                    "Show top 50 recurring age-associated genes and per-study estimates",
                    max_height=480,
                )

        if not intersection_diagnostics.empty:
            display(Markdown(
                "<a id=\"per-study-gene-intersection\"></a>\n"
                "### Per-Study Gene-Intersection Diagnostic"
            ))
            display(Markdown(
                "This table shows how much restricting a cross-study comparison to genes tested "
                "by every available study would remove from each study's significant-gene list. "
                "It is a gene-universe sensitivity check, not a replication result. Each row uses "
                "the displayed single-study design and complete-case sample count; age, sex, BMI, "
                "and CMV coefficients from the same study therefore come from that one maximal model."
            ))
            display_collapsible_table(intersection_diagnostics.assign(
                **{"outside intersection (%)": intersection_diagnostics[
                    "outside intersection (%)"
                ].map(lambda value: f"{value:.1f}%")}
            ), "Show per-study gene-intersection sensitivity details")

        if age_trajectory_fit or age_trajectory_results is not None:
            display(Markdown(
                "<a id=\"age-bin-trajectories\"></a>\n## Age-Bin Trajectories"
            ))
            if de_dir:
                cross_type_report = (
                    Path(de_dir).resolve().parent / "trajectory_analysis" / "report.html"
                )
                if cross_type_report.is_file():
                    display(Markdown(
                        "See the [cross-cell-type trajectory report](../../differential_expression/"
                        "trajectory_analysis/report.html) for "
                        "shared-pattern clusters and gene recurrence across cell types."
                    ))
            width = trajectory_settings.bin_width_years
            minimum_bin_samples = trajectory_settings.minimum_samples_per_bin
            minimum_bins = trajectory_settings.minimum_bins
            strict_age_cutoff = trajectory_settings.strict_age_cutoff_exclusive
            age_cutoff_excluded = age_trajectory_info.get(
                "samples_excluded_by_age_cutoff", 0
            )
            manually_excluded_bins = age_trajectory_info.get(
                "manually_excluded_age_bins", []
            )
            raw_bin_counts = age_trajectory_info.get(
                "age_bin_counts_before_filtering", {}
            )
            manual_bin_counts_text = ", ".join(
                f"{age_bin}: {raw_bin_counts.get(age_bin, 0)} eligible pseudobulks"
                for age_bin in manually_excluded_bins
            ) or "none"
            reference_bin = age_trajectory_info.get(
                "reference_bin",
                f"{trajectory_settings.reference_bin_start_age}-"
                f"{trajectory_settings.reference_bin_start_age + width}",
            )
            retained_bins = age_trajectory_info.get("retained_bins", [])
            display(Markdown(
                f"Trajectory sample filtering: manually excluded age bin(s) "
                f"{manual_bin_counts_text}. The strict age cutoff is < {strict_age_cutoff} years; "
                f"{age_cutoff_excluded} eligible pseudobulks aged at or above it were excluded "
                "(this total includes the manually excluded 90–100 bin). Each retained "
                f"{width}-year bin must contain at least {minimum_bin_samples} eligible pseudobulks."
            ))
            if age_trajectory_fit.get("status") != "complete" or age_trajectory_results is None:
                display(Markdown(
                    f"_Age-bin trajectory fit skipped: "
                    f"{age_trajectory_fit.get('reason', 'no trajectory result was produced')}._"
                ))
            elif len(retained_bins) < trajectory_settings.minimum_bins:
                display(Markdown(
                    f"_Age-bin trajectory results omitted: only {len(retained_bins)} bins "
                    f"meet the sample-count threshold; at least "
                    f"{trajectory_settings.minimum_bins} "
                    "are required._"
                ))
            else:
                dropped_bins = age_trajectory_info.get("dropped_bins", [])
                top_gene_count = trajectory_settings.report_top_n_genes
                trajectory_de_fdr = trajectory_settings.de_fdr_threshold
                display(Markdown(
                    f"The combined model uses {width}-year age bins, adjusts for study and sex, "
                    f"and uses **{reference_bin}** as the zero point. Bins with fewer than "
                    f"{minimum_bin_samples} eligible pseudobulks are excluded; at least "
                    f"{minimum_bins} supported bins are required. Retained bins: "
                    f"{', '.join(retained_bins)}. "
                    f"Dropped bins: {', '.join(dropped_bins) if dropped_bins else 'none'}. "
                    "The omnibus joint Wald test asks whether any retained age bin differs from "
                    "the reference, with FDR controlled across genes. Each trajectory is scaled "
                    "by its population SD across all retained bins, including the zero reference; "
                    "a trajectory with incomplete estimates or zero SD has no standardized shape."
                ))
                omnibus = age_trajectory_results.drop_duplicates("gene").copy()
                omnibus["omnibus_padj"] = pd.to_numeric(
                    omnibus["omnibus_padj"], errors="coerce"
                )
                selected_genes = omnibus.loc[
                    omnibus["omnibus_padj"] < trajectory_de_fdr
                ].sort_values("omnibus_padj").head(top_gene_count)
                if selected_genes.empty:
                    display(Markdown(
                        f"_No genes pass the omnibus age-bin FDR threshold "
                        f"(< {trajectory_de_fdr:g}) for this cell type._"
                    ))
                else:
                    selected_names = selected_genes["gene"].astype(str).tolist()
                    selected_shape = age_trajectory_results.loc[
                        age_trajectory_results["gene"].astype(str).isin(selected_names)
                    ].copy()
                    shape_matrix = selected_shape.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(index=selected_names, columns=retained_bins)
                    figure, axis = plt.subplots(
                        figsize=(max(7.0, 1.1 * len(retained_bins)),
                                 max(3.5, 0.28 * len(selected_names) + 1.5))
                    )
                    sns.heatmap(
                        shape_matrix, cmap="vlag", center=0, ax=axis,
                        cbar_kws={"label": "Trajectory shape (SD units)"},
                    )
                    axis.set_xlabel("Age bin (years)")
                    axis.set_ylabel("Gene")
                    axis.set_title("Top omnibus-significant standardized trajectories")
                    figure.tight_layout()
                    display_collapsible_figure(figure, "Show omnibus-significant age trajectories")

                    selected_table = selected_shape.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(index=selected_names, columns=retained_bins)
                    selected_table.insert(
                        0, "Omnibus FDR", selected_genes.set_index("gene").loc[
                            selected_table.index, "omnibus_padj"
                        ].to_numpy(),
                    )
                    selected_table.insert(
                        1, "Trajectory SD (log2 fold change)", selected_genes.set_index(
                            "gene"
                        ).loc[selected_table.index, "trajectory_sd"].to_numpy(),
                    )
                    selected_table.columns = [
                        str(column) if column in {"Omnibus FDR", "Trajectory SD (log2 fold change)"}
                        else f"Shape {column}"
                        for column in selected_table.columns
                    ]
                    display_collapsible_table(
                        selected_table.reset_index(),
                        "Show top omnibus-significant age-bin trajectories",
                    )

                cluster_fdr = trajectory_settings.cluster_fdr_threshold
                cluster_significant = omnibus.loc[
                    omnibus["omnibus_padj"] <= cluster_fdr
                ]
                if cluster_significant.empty:
                    display(Markdown(
                        f"_No genes pass the trajectory clustering FDR threshold "
                        f"({cluster_fdr:g})._"
                    ))
                else:
                    profile_data = age_trajectory_results.loc[
                        age_trajectory_results["gene"].astype(str).isin(
                            cluster_significant["gene"].astype(str)
                        )
                    ].copy()
                    profiles = profile_data.pivot(
                        index="gene", columns="age_bin", values="trajectory_scaled"
                    ).reindex(columns=retained_bins).dropna()
                    cluster_labels, cluster_means = cluster_trajectory_profiles(
                        profiles,
                        trajectory_settings,
                    )
                    merge_diagnostics = trajectory_merge_diagnostics(
                        profiles, trajectory_settings, maximum_merges=20
                    )
                    cluster_rows = profiles.join(cluster_labels).reset_index()
                    cluster_rows.insert(1, "cell_type", cell_type_name)
                    cluster_rows["omnibus_padj"] = cluster_rows["gene"].map(
                        omnibus.set_index("gene")["omnibus_padj"]
                    )
                    cluster_rows.to_csv(output_dir / "age_trajectory_clusters.csv", index=False)
                    cluster_means = cluster_means.reset_index()
                    cluster_means.insert(0, "cell_type", cell_type_name)
                    cluster_means.to_csv(
                        output_dir / "age_trajectory_cluster_means.csv", index=False
                    )
                    figure, axes = plt.subplots(1, 3, figsize=(18, 4.8))
                    palette = sns.color_palette(
                        "husl", n_colors=max(len(cluster_means), 1)
                    )
                    cluster_colors = {
                        str(int(row["trajectory_cluster"])): palette[index]
                        for index, row in cluster_means.iterrows()
                    }
                    for index, row in cluster_means.iterrows():
                        cluster = int(row["trajectory_cluster"])
                        label = f"Cluster {cluster} (n={int(row['n_trajectories'])})"
                        axes[0].plot(
                            retained_bins,
                            [row[age_bin] for age_bin in retained_bins],
                            marker="o", linewidth=2, color=cluster_colors[str(cluster)],
                            label=label,
                        )
                    axes[0].axhline(0, color="#555555", linestyle="--", linewidth=0.8)
                    axes[0].set_xlabel("Age bin")
                    axes[0].set_ylabel("Mean standardized trajectory (SD units)")
                    axes[0].set_title("Mean trajectory by hierarchical cluster")
                    axes[0].tick_params(axis="x", rotation=35)
                    axes[0].legend(fontsize=8)
                    embedded = cluster_labels.dropna(subset=["umap_1", "umap_2"])
                    if embedded.empty:
                        axes[1].text(
                            0.5, 0.5,
                            f"At least {trajectory_settings.minimum_umap_trajectories} "
                            "trajectories are needed for UMAP",
                            ha="center", va="center", transform=axes[1].transAxes,
                        )
                    else:
                        plot_points = embedded.reset_index()
                        plot_points["trajectory_cluster"] = plot_points[
                            "trajectory_cluster"
                        ].astype(str)
                        sns.scatterplot(
                            data=plot_points, x="umap_1", y="umap_2",
                            hue="trajectory_cluster", palette=cluster_colors, s=24,
                            alpha=0.8, linewidth=0, ax=axes[1], legend="brief",
                        )
                        axes[1].set_title(
                            f"Trajectory UMAP ({trajectory_settings.distance_metric} "
                            "bin distance)"
                        )
                        axes[1].set_xlabel("UMAP 1")
                        axes[1].set_ylabel("UMAP 2")
                    plot_trajectory_merge_diagnostic(
                        axes[2], merge_diagnostics,
                        max_clusters=trajectory_settings.max_clusters,
                    )
                    figure.tight_layout()
                    display_collapsible_figure(
                        figure,
                        f"Show trajectory clusters, UMAP, and merge diagnostic ({len(profiles):,} genes, "
                        f"FDR ≤ {cluster_fdr:g})",
                        open_by_default=True,
                    )
                    display(Markdown(
                        "The merge-cost panel shows up to the final 20 hierarchical merges. "
                        "A sharp increase at **K → K−1** supports retaining K trajectory groups; "
                        "the red marker shows the configured cluster cap when that transition is available."
                    ))
                    display(Markdown(
                        "Download the [gene cluster assignments and standardized profiles]"
                        "(age_trajectory_clusters.csv) or [cluster mean trajectories]"
                        "(age_trajectory_cluster_means.csv)."
                    ))

        display(Markdown(
            "<a id=\"sample-residual-clustering\"></a>\n"
                "## Per-Sample Residual Clustering"
        ))
        display(Markdown(
            "These sample profiles come from the merged age-bin model. It replaces continuous "
            "age with configured decade bins and includes sex, plus study-site and "
            "log10(total counts) terms when estimable. For each sample and gene, the saved "
            "Pearson residual is the observed count minus the fitted mean, divided by the "
            "fitted standard deviation \\(\\sqrt{\\mu + \\alpha\\mu^2}\\); it is a scaled "
            "deviation from the model fit. PCA centers each gene's residuals across samples, "
            "then retains up to the first 100 principal components, capped by the available "
            "samples and genes. Those PC scores are the input to the two-dimensional UMAP "
            "shown below. Colors overlay sample covariates and residual-cluster labels on "
            "the same UMAP coordinates."
        ))
        display(Markdown("### Residual UMAP Covariate Overlays"))
        if de_dir and age_trajectory_fit:
            (
                residual_clusters, residual_pc_scores, residual_pcs, residual_matrix,
            ) = load_residual_sample_clusters(de_dir, age_trajectory_fit, trajectory_settings)
        else:
            residual_clusters, residual_pc_scores, residual_pcs, residual_matrix = (
                None, None, 0, None
            )
        if residual_clusters is None:
            residual_reason = age_trajectory_fit.get(
                "reason", "the age-bin model did not save Pearson residuals"
            )
            display(Markdown(f"_Residual clustering unavailable: {residual_reason}._"))
        else:
            residual_clusters.to_csv(output_dir / "residual_sample_clusters.csv", index=False)
            residual_pc_scores.to_csv(
                output_dir / "residual_sample_pc_scores.csv", index=False
            )
            display(Markdown(
                f"The embedding includes **{len(residual_clusters)} samples** and "
                f"**{residual_pc_scores.shape[1] - 1:,} genes**, with {residual_pcs} "
                "principal components retained. When enough samples are available, "
                "exploratory K-means cluster labels are selected by silhouette score on "
                "the two-dimensional UMAP."
            ))
            plot_residual_sample_umaps(
                residual_clusters,
                minimum_umap_samples=trajectory_settings.minimum_umap_trajectories,
            )
            plot_residual_cluster_markers(residual_matrix, residual_clusters)
            display(Markdown(
                "Download the [residual cluster assignments](residual_sample_clusters.csv) "
                "or [sample PCA scores](residual_sample_pc_scores.csv)."
            ))

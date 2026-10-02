# %%
import json
import os
import warnings
from html import escape
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

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 110})
from pbmc_pipeline.config import read_json
from pbmc_pipeline.covariates import (
    categorical_contrast_label,
    categorical_levels,
    normalize_cmv_status,
)
from pbmc_pipeline.reporting import (
    sample_cell_type_fractions,
    sample_cluster_fractions,
    signed_value_at_largest_absolute_magnitude,
)

primitive = sc.read_h5ad(Path(os.environ["CELL_TYPE_H5AD"]))
analysis_path = Path(os.environ["CELL_TYPE_ANALYSIS_H5AD"])
output_dir = analysis_path.parent
pipeline = read_json(Path(os.environ["CELL_TYPE_CONFIG"]))
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
summary_metrics_html = "".join(
    f"<div class=\"report-metric\"><dt>{escape(label)}</dt><dd>{escape(value)}</dd></div>"
    for label, value in summary_metrics
)
de_toc_link = (
    '<a href="#differential-expression">Pseudobulk differential expression</a>'
    if de_dir else ""
)
display(HTML(f"""
<style>
:root {{
  --report-ink: #182230;
  --report-muted: #5f6b7a;
  --report-line: #d8e0e8;
  --report-paper: #ffffff;
  --report-page: #f5f7fa;
  --report-accent: #176b87;
  --report-accent-soft: #e6f3f7;
}}
body.jp-Notebook {{
  background: var(--report-page);
  color: var(--report-ink);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  font-size: 17px;
}}
main {{ max-width: 1600px; margin: 0 auto; padding: 2rem clamp(1.25rem, 3vw, 2.5rem) 4rem; }}
.jp-Cell {{ margin: 1.4rem 0; }}
.jp-Cell-outputWrapper, .jp-OutputArea {{ width: 100%; }}
.jp-RenderedImage img {{ display: block; max-width: 100%; height: auto; margin: 0 auto; }}
.jp-RenderedMarkdown h2 {{
  margin-top: 3rem; padding-bottom: .55rem; border-bottom: 2px solid var(--report-line);
  color: var(--report-ink); font-size: clamp(1.45rem, 2.7vw, 2rem);
}}
.jp-RenderedMarkdown h3 {{ color: var(--report-accent); margin-top: 2.25rem; }}
.jp-RenderedMarkdown p, .jp-RenderedMarkdown li {{ color: #354152; font-size: 1.04rem; line-height: 1.65; }}
.jp-RenderedMarkdown blockquote {{
  margin: 1rem 0; padding: .75rem 1rem; border-left: 4px solid #bd6a29;
  background: #fff6eb; color: #6e3e18;
}}
.report-hero {{
  margin: 0 0 2.5rem; padding: clamp(1.4rem, 4vw, 2.6rem); border-radius: 18px;
  background: linear-gradient(135deg, #315c6d, #527990); color: #fff;
  box-shadow: 0 12px 30px rgba(49, 92, 109, .16);
}}
.report-eyebrow {{ margin: 0 0 .45rem; color: #bde7f0; font-size: .76rem; font-weight: 700; letter-spacing: .11em; text-transform: uppercase; }}
.report-hero h1 {{ margin: 0; color: #fff; font-size: clamp(2rem, 5vw, 3.4rem); line-height: 1.08; }}
.report-subtitle {{ max-width: 48rem; margin: .85rem 0 1.4rem; color: #e2f3f7; line-height: 1.55; }}
.report-metrics {{ display: flex; flex-wrap: wrap; align-items: flex-start; gap: .55rem; margin: 0; }}
.report-metric {{ display: flex; align-items: baseline; gap: .55rem; padding: .62rem .78rem; border: 1px solid rgba(255,255,255,.22); border-radius: 10px; background: rgba(255,255,255,.1); }}
.report-metric dt {{ color: #d6edf2; font-size: .9rem; font-weight: 600; }}
.report-metric dd {{ margin: 0; color: #fff; font-size: .9rem; font-weight: 650; text-align: right; }}
.report-toc {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(175px, 1fr)); gap: .75rem; margin-top: 1.5rem; }}
.report-toc-group {{ padding: .7rem .8rem; border: 1px solid rgba(255,255,255,.23); border-radius: 10px; background: rgba(16, 42, 53, .15); }}
.report-toc-title {{ margin: 0 0 .38rem; color: #cbeaf1; font-size: .72rem; font-weight: 750; letter-spacing: .07em; text-transform: uppercase; }}
.report-toc-links {{ display: flex; flex-direction: column; gap: .25rem; }}
.report-hero .report-toc a, .report-hero .report-toc a:visited {{
  color: #f7fcfd; font-size: .91rem; font-weight: 600; line-height: 1.35; text-decoration: none;
}}
.report-hero .report-toc a::before {{ content: "›"; display: inline-block; width: .8rem; color: #bde7f0; }}
.report-hero .report-toc a:hover {{ color: #fff; text-decoration: underline; text-underline-offset: .16em; }}
.report-provenance {{ padding: .75rem 1rem; border: 1px solid var(--report-line); border-radius: 10px; background: var(--report-paper); }}
.report-provenance summary {{ cursor: pointer; color: var(--report-accent); font-weight: 700; }}
.report-provenance table {{ margin-top: .8rem; }}
.report-evidence {{ margin: 2.5rem 0 1.25rem; padding: 1.25rem 1.4rem; border: 1px solid #c8dde5; border-radius: 14px; background: #f1f8fa; }}
.report-evidence h2 {{ margin: 0 0 .45rem; color: var(--report-ink); font-size: 1.45rem; }}
.report-evidence p {{ margin: 0 0 .85rem; color: #354152; line-height: 1.55; }}
.report-evidence-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: .65rem; }}
.report-evidence-item {{ padding: .7rem .8rem; border-radius: 9px; background: #fff; border: 1px solid #d7e7ec; }}
.report-evidence-item dt {{ color: var(--report-muted); font-size: .77rem; font-weight: 700; text-transform: uppercase; letter-spacing: .04em; }}
.report-evidence-item dd {{ margin: .22rem 0 0; color: var(--report-ink); font-weight: 650; line-height: 1.35; }}
.jp-RenderedHTMLCommon {{ max-width: 100%; overflow-x: auto; }}
.jp-RenderedHTMLCommon table {{ border-collapse: collapse; background: var(--report-paper); }}
.jp-RenderedHTMLCommon th {{ background: #edf3f6; color: #263544; font-weight: 700; }}
.jp-RenderedHTMLCommon th, .jp-RenderedHTMLCommon td {{ padding: .45rem .65rem; border: 1px solid var(--report-line); vertical-align: top; }}
@media print {{
  body.jp-Notebook {{ background: #fff; }} main {{ max-width: none; padding: 0; }}
  .report-hero {{ box-shadow: none; }}
}}
</style>
<header class=\"report-hero\">
  <p class=\"report-eyebrow\">PBMC ageing · per-cell-type report</p>
  <h1>{escape(cell_type_name)}</h1>
  <p class=\"report-subtitle\">Descriptive sample-level and single-cell diagnostics. Cell-level associations are not independent-sample inference.</p>
  <dl class=\"report-metrics\">{summary_metrics_html}</dl>
  <nav class=\"report-toc\" aria-label=\"Report table of contents\">
    <section class=\"report-toc-group\">
      <p class=\"report-toc-title\">Study support</p>
      <div class=\"report-toc-links\"><a href=\"#sample-coverage\">Sample coverage and protocol metadata</a></div>
    </section>
    <section class=\"report-toc-group\">
      <p class=\"report-toc-title\">Fraction analysis</p>
      <div class=\"report-toc-links\">
        <a href=\"#sample-fractions\">Sample fractions across age</a>
        <a href=\"#fraction-model-evidence\">Adjusted fraction-model evidence</a>
      </div>
    </section>
    <section class=\"report-toc-group\">
      <p class=\"report-toc-title\">Cell-state context</p>
      <div class=\"report-toc-links\">
        <a href=\"#global-location\">Global and local embeddings</a>
        <a href=\"#cluster-composition\">Cluster composition</a>
        <a href=\"#cluster-markers\">Cluster markers</a>
      </div>
    </section>
    <section class=\"report-toc-group\">
      <p class=\"report-toc-title\">Technical and gene-level diagnostics</p>
      <div class=\"report-toc-links\">
        <a href=\"#pc-study-technology\">PCA by study and protocol</a>
        <a href=\"#pc-age\">PCA loadings and PC–age</a>
        {de_toc_link}
      </div>
    </section>
  </nav>
</header>
"""))
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


def plot_umap(adata, *, color, title=None, label_clusters=False):
    """Render comparable low-resolution, equal-aspect UMAP panels."""
    plot_kwargs = {"color": color, "size": max(2, 100_000 / adata.n_obs), "show": False}
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
    plt.show()


def plot_combined_study_fits(fraction_specs):
    """Plot equal-study standardized predictions from study-adjusted linear models."""
    figure, axes = plt.subplots(
        1, len(fraction_specs), figsize=(7.8 * len(fraction_specs), 5.5), squeeze=False,
    )
    for axis, (data, fraction_column, label) in zip(axes.flat, fraction_specs):
        model_data = data[["study", "age", fraction_column]].copy()
        model_data["study"] = model_data["study"].astype(str)
        model_data["age"] = pd.to_numeric(model_data["age"], errors="coerce")
        model_data[fraction_column] = pd.to_numeric(model_data[fraction_column], errors="coerce")
        model_data = model_data.dropna().reset_index(drop=True)
        study_ranges = model_data.groupby("study", observed=True)["age"].agg(["min", "max"])
        variable_studies = model_data.groupby("study", observed=True)["age"].nunique()
        if model_data.empty or not (variable_studies >= 2).any():
            axis.text(
                0.5, 0.5, "Insufficient within-study age variation for a shared slope",
                transform=axis.transAxes, ha="center", va="center",
            )
        else:
            study_terms = pd.get_dummies(
                model_data["study"], prefix="study", drop_first=True, dtype=float,
            )
            design = pd.concat(
                [model_data[["age"]].astype(float), study_terms], axis=1,
            )
            design = sm.add_constant(design, has_constant="add")
            rank = np.linalg.matrix_rank(design.to_numpy(dtype=float))
            if len(model_data) <= rank:
                axis.text(
                    0.5, 0.5, "Insufficient residual degrees of freedom for interval",
                    transform=axis.transAxes, ha="center", va="center",
                )
            else:
                study_sample_weights = model_data.groupby(
                    "study", observed=True,
                )["study"].transform(lambda values: 1 / len(values))
                model = sm.WLS(
                    model_data[fraction_column].to_numpy(dtype=float), design,
                    weights=study_sample_weights.to_numpy(dtype=float),
                ).fit(cov_type="HC3")
                age_grid = np.linspace(model_data["age"].min(), model_data["age"].max(), 160)
                standardized_design = pd.DataFrame(
                    0.0, index=np.arange(len(age_grid)), columns=design.columns,
                )
                standardized_design["const"] = 1.0
                standardized_design["age"] = age_grid
                for column in study_terms.columns:
                    study = column.removeprefix("study_")
                    standardized_design[column] = np.mean(
                        model_data["study"].to_numpy() == study
                    )
                prediction = model.get_prediction(standardized_design).summary_frame()

                overlap_min = study_ranges["min"].max()
                overlap_max = study_ranges["max"].min()
                if overlap_min <= overlap_max:
                    axis.axvspan(
                        overlap_min, overlap_max, color="#176b87", alpha=0.08,
                        label="All studies observed",
                    )
                axis.fill_between(
                    age_grid,
                    prediction["mean_ci_lower"].to_numpy(),
                    prediction["mean_ci_upper"].to_numpy(),
                    color="#176b87", alpha=0.18, linewidth=0,
                    label="95% HC3 confidence interval",
                )
                axis.plot(
                    age_grid, prediction["mean"].to_numpy(),
                    color="#124e63", linewidth=3,
                    label="Equal-study adjusted estimate",
                )
                support_note = (
                    "Shading marks the age range observed in every study; outside it, "
                    "predictions extrapolate for some studies."
                    if overlap_min <= overlap_max
                    else "No age range is observed in every study; the standardized curve extrapolates."
                )
                axis.text(
                    0.02, 0.02, support_note, transform=axis.transAxes,
                    ha="left", va="bottom", fontsize=8, color="#4d5966",
                )
                axis.legend(loc="best", fontsize=8)
        axis.set_title(label)
        axis.set_xlabel("Age (years)")
        axis.set_ylabel(label)
    figure.suptitle("Study-adjusted linear fraction estimates", y=0.98)
    figure.text(
        0.5, 0.015,
        "Weighted least squares with study-specific intercepts and one shared age slope; "
        "each study has equal total model weight. Intervals use HC3 robust standard errors.",
        ha="center", fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.06, 1, 0.92))
    plt.show()


def plot_cluster_fractions(data, *, palette):
    """Keep study identity in points while showing one uncluttered pooled trend."""
    if data.empty:
        display(Markdown("_No samples meet the denominator cutoff for cluster fractions._"))
        return
    plot_data = data.assign(study=data["study"].astype(str))
    clusters = sorted(plot_data["cluster"].astype(str).unique())
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
    plt.show()


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


def plot_parent_fraction_coefficients(coefficients, residual_diagnostics, *, palette):
    """Plot adjusted coefficients plus observed and binomial residual variation."""
    if coefficients.empty:
        display(Markdown("_No estimable adult-only study models for this fraction._"))
        return
    covariate_order = [name for name in ["age", "sex", "bmi", "cmv"]
                       if name in coefficients["covariate"].unique()]
    studies = sorted(coefficients["study"].unique())
    residual_label = "residual_sd"
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
    figure, axis = plt.subplots(figsize=(10, max(4.8, 1.15 * (len(coefficient_rows) + 1))))
    for study_index, study in enumerate(studies):
        subset = coefficients.loc[coefficients["study"] == study]
        for _, row in subset.iterrows():
            offset = (study_index - (len(studies) - 1) / 2) * 0.055
            y = y_positions[(row["covariate"], row["contrast"])] + offset
            axis.errorbar(
                row["estimate"], y,
                xerr=[[row["estimate"] - row["ci_low"]], [row["ci_high"] - row["estimate"]]],
                fmt="o", color=palette[study], capsize=3, markersize=5,
            )
    for study_index, study in enumerate(studies):
        subset = residual_diagnostics.loc[residual_diagnostics["study"] == study]
        if subset.empty:
            continue
        row = subset.iloc[0]
        # Interpretation check: assess whether residual variation is too large
        # for the coefficient estimates above to support useful interpretation.
        offset = (study_index - (len(studies) - 1) / 2) * 0.055
        y = residual_y + offset
        axis.plot(row["residual_sd"], y, "o", color=palette[study], markersize=6)
        axis.plot(
            row["binomial_sampling_sd"], y, marker="D", linestyle="", markersize=5,
            markerfacecolor="none", markeredgecolor=palette[study], markeredgewidth=1.25,
        )
    axis.axvline(0, color="#333333", linewidth=1, linestyle="--")
    axis.set_yticks(
        list(range(residual_y + 1)),
        [*(label for _, _, label in coefficient_rows), "RESIDUAL SD (fraction)"],
    )
    axis.set_xlabel("OLS coefficient (95% CI)")
    axis.set_ylabel("Covariate / diagnostic")
    axis.invert_yaxis()
    study_legend = axis.legend(
        handles=[plt.Line2D([], [], marker="o", linestyle="", color=palette[study], label=study)
                 for study in studies],
        title="Study (AGE = per decade)", bbox_to_anchor=(1.02, 1), loc="upper left",
    )
    axis.add_artist(study_legend)
    axis.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="", color="#333333", label="Observed model residual SD"),
            plt.Line2D([], [], marker="D", linestyle="", color="#333333", markerfacecolor="none",
                       label="Conditional binomial sampling SD (bootstrap mean)"),
        ],
        title="Residual diagnostic", bbox_to_anchor=(1.02, 0.53), loc="upper left",
    )
    axis.set_title("Adult sample-level fraction model coefficients")
    figure.tight_layout()
    plt.show()


def directional_cluster_genes(markers, *, cluster, n_genes=8):
    """Return the strongest positive and negative genes for one cluster contrast."""
    subset = markers.loc[markers["group"].astype(str) == str(cluster)].copy()
    positive = subset.loc[subset["logfoldchanges"] > 0].nlargest(n_genes, "logfoldchanges")
    negative = subset.loc[subset["logfoldchanges"] < 0].nsmallest(n_genes, "logfoldchanges")
    result = pd.concat([negative, positive]).copy()
    result["direction"] = result["logfoldchanges"].ge(0).map({
        True: "Higher in cluster", False: "Lower in cluster",
    })
    return result.sort_values("logfoldchanges")


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
    plt.show()

# %% [markdown]
# <a id="sample-coverage"></a>
# ## Sample coverage among retained samples
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
display(study_coverage)
figure, axis = plt.subplots(figsize=(10, 4))
sns.stripplot(data=sample_covariates, x="study", y="age", hue="sex", dodge=True, ax=axis)
axis.set_title("Age coverage by study and sex")
axis.tick_params(axis="x", rotation=45)
figure.tight_layout()
plt.show()

# %% [markdown]
# <a id="sample-fractions"></a>
# ## Sample-level fraction across age
#
# Each point is one study × sample. The left panel uses all retained PBMCs as
# denominator; the right panel uses the manually configured AIFI-L1 parent
# derived from every cell's merged L2 call, never from a separate L1 classifier.
# Samples with fewer than the configured denominator-cell cutoff are excluded.
# Facets show descriptive within-study linear fits for all samples and for age
# ≥20. The combined view is a study-adjusted descriptive fit: one weighted
# least-squares model with study-specific intercepts and a shared age slope,
# equal total weight per study, and HC3 confidence intervals. Shading identifies
# the age interval observed in every study; predictions outside it extrapolate
# for at least one study. HC3 intervals treat sample rows as independent;
# repeat samples from the same subject are not clustered, which may understate
# uncertainty.

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
display(pd.DataFrame({
    "fraction": [spec[2] for spec in fraction_specs],
    "eligible_samples": [len(spec[0]) for spec in fraction_specs],
    "excluded_below_denominator_cutoff": [
        total - len(spec[0]) for total, spec in zip(all_fraction_samples, fraction_specs)
    ],
    "minimum_denominator_cells": min_fraction_denominator,
}))

studies = sorted(pd.concat([spec[0] for spec in fraction_specs])["study"].unique())
palette = dict(zip(studies, sns.color_palette("tab10", n_colors=len(studies))))
display(Markdown(f"### 1. Linear Age Models per Study ({cell_type_name})"))
for index, (fractions, fraction_column, title, ylabel) in enumerate(fraction_specs):
    plot_study_linear_trends(
        fractions, fraction_column=fraction_column, title=title, ylabel=ylabel, palette=palette,
        show_fit_legend=index == 0,
    )
display(Markdown(f"### 2. Study-adjusted Linear Fraction Estimate ({cell_type_name})"))
plot_combined_study_fits([
    (fractions, fraction_column, ylabel)
    for fractions, fraction_column, _, ylabel in fraction_specs
])

# %%
display(Markdown(
    f"### 3. Fraction of {cell_type_name} in {parent}: Linear Models, All Covariates"
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
# scale, those effects may not support a reasonable interpretation.

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
display(Markdown("#### Fraction-model residual and sampling diagnostics"))
display(fraction_residual_diagnostics)

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
    "<h2>Fraction-model evidence at a glance</h2>"
    "<p>These are study-specific adult-sample OLS fits, not a pooled effect estimate. "
    "Use them to assess direction, precision, and the residual diagnostic before interpreting individual coefficients.</p>"
    "<dl class=\"report-evidence-grid\">"
    f"<div class=\"report-evidence-item\"><dt>Estimable study models</dt><dd>{len(fraction_residual_diagnostics)}</dd></div>"
    f"<div class=\"report-evidence-item\"><dt>Age-effect direction</dt><dd>{escape(age_direction)}</dd></div>"
    f"<div class=\"report-evidence-item\"><dt>Residual check</dt><dd>{escape(residual_summary)}</dd></div>"
    f"<div class=\"report-evidence-item\"><dt>Covariates represented</dt><dd>{escape(', '.join(model_covariates) or 'None')}</dd></div>"
    "</dl></section>"
))
plot_parent_fraction_coefficients(
    fraction_coefficients, fraction_residual_diagnostics, palette=palette,
)
display(Markdown(
    "_How to read this plot: each colour is a study and intervals are its 95% OLS confidence "
    "intervals. The filled dot on the residual row is the observed post-model variation; the "
    "open diamond is the conditional cell-sampling reference. Compare studies for consistency "
    "and precision, rather than treating their collection as a pooled estimate._"
))

# %% [markdown]
# <a id="global-location"></a>
# ## Location on the global embedding
#
# These coordinates are inherited from the cross-study embedding and restricted
# to the current cell type. They provide context for where the type sits in the
# global PBMC landscape; this section does not fit a new embedding.

# %%
if "X_umap" in primitive.obsm:
    plot_umap(primitive, color="study", title="Former global UMAP coordinates, restricted to this cell type")
    display(Markdown(
        "_This is orientation context only: its coordinates come from the global embedding, and "
        "colours show study membership rather than a tested biological association._"
    ))

# %%
if report["status"] != "complete":
    display(Markdown("_Too few cells for a type-specific embedding and clustering._"))
else:
    plot_umap(adata, color="cluster", title="Local UMAP: cluster IDs", label_clusters=True)
    colors = ["study", "age"]
    if "aifi_l3_majority" in adata.obs:
        colors.append("aifi_l3_majority")
    for color in colors:
        plot_umap(adata, color=color)

# %% [markdown]
# <a id="cluster-composition"></a>
# ## Cluster composition
#
# The table shows local-cluster cell counts by study. The plots that follow use
# one study × sample fraction per point, so they describe sample-level cluster
# composition rather than treating cells as independent observations.

# %%
if report["status"] == "complete":
    composition = pd.read_csv(output_dir / "cluster_study_composition.tsv", sep="\t")
    display(composition.pivot(index="cluster", columns="study", values="n_cells").fillna(0).astype(int))
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
# <a id="cluster-markers"></a>
# ## Cluster markers
#
# Genes are ranked by a Wilcoxon comparison of each local cluster against all
# other local clusters. Signed log fold changes therefore show whether a gene
# is higher or lower in that cluster; these are descriptive contrasts.

# %%
if report["status"] == "complete":
    markers = pd.read_csv(output_dir / "markers.tsv", sep="\t")
    clusters = sorted(markers["group"].astype(str).unique())
    n_columns = min(3, len(clusters))
    n_rows = (len(clusters) + n_columns - 1) // n_columns
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5.5 * n_columns, 4.4 * n_rows), squeeze=False)
    for axis, cluster in zip(axes.flat, clusters):
        subset = directional_cluster_genes(markers, cluster=cluster)
        sns.barplot(
            data=subset, x="logfoldchanges", y="names", hue="direction", dodge=False,
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
    plt.show()

# %% [markdown]
# <a id="pc-study-technology"></a>
# ## PCA scores by study and protocol covariates
#
# These are the native, local PCA scores before Harmony adjustment. Colour uses
# the report's standard study palette and marker shape indicates the recorded
# single-cell technology in the top row. The lower row collapses technology to
# 3′ or 5′ marker shapes and maps intronic read inclusion to colour. The
# Harmony-adjusted representation is used only for the local neighbour graph,
# UMAP, and clustering.

# %%
if report["status"] == "complete":
    plot_pca_study_technology_and_introns(adata)
    display(Markdown(
        "_Read these panels as a technical-structure diagnostic: separation by study, technology, "
        "or intronic alignment can indicate design effects, but is not itself a biological test._"
    ))

# %% [markdown]
# <a id="pc-age"></a>
# ## PCA gene loadings, variance explained, and exploratory PC–age correlations
#
# Loadings show which genes contribute most to each native local PC, while the
# variance plot summarizes their cumulative contribution. PC–age correlations
# are descriptive cell-level summaries and are not subject-level inference.

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
    plt.show()

# %%
if report["status"] == "complete":
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
    cumulative_variance = np.cumsum(variance_ratio)
    figure, axis = plt.subplots(figsize=(9, 4.5))
    axis.plot(np.arange(1, len(cumulative_variance) + 1), cumulative_variance, marker="o", markersize=3)
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Cumulative variance explained")
    axis.set_ylim(0, 1.02)
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    plt.show()

# %%
if report["status"] == "complete":
    pc_age = pd.read_csv(output_dir / "pc_age_correlations.tsv", sep="\t")
    display(pc_age.reindex(pc_age["spearman_r"].abs().sort_values(ascending=False).index).head(10))
    plt.figure(figsize=(9, 4.5))
    sns.barplot(data=pc_age, x="pc", y="spearman_r", color="#4c72b0")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("Spearman correlation with age (cells; descriptive)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

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
    if merged_file is not None and merged_file.is_file():
        result = pd.read_csv(merged_file)
        if "covariate" not in result:
            result["covariate"] = "age"
        if "contrast" not in result:
            result["contrast"] = "per year"
        combined_results_by_covariate["age"] = result
    for result_path in combined_files:
        result = pd.read_csv(result_path)
        if "covariate" not in result:
            result["covariate"] = result_path.stem
        if "contrast" not in result:
            result["contrast"] = "per BMI unit" if result_path.stem == "bmi" else result_path.stem
        combined_results_by_covariate[result_path.stem] = result
    has_combined = bool(combined_results_by_covariate)
    if has_per_study or has_combined or de_run_metadata:
        display(Markdown("<a id=\"differential-expression\"></a>\n## Pseudobulk differential expression"))
        display(Markdown(
            "Pseudobulk models aggregate counts at the study × sample × cell-type level. "
            "They provide the sample-level complement to the descriptive single-cell views above."
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
    if has_per_study or has_combined:
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
        if not intersection_diagnostics.empty:
            display(Markdown(
                "For each covariate, each study's significant genes are compared with the "
                "intersection of genes tested by all available per-study fits for that "
                "covariate. The table reports the count and fraction outside the intersection."
            ))
            display(intersection_diagnostics.assign(
                **{"outside intersection (%)": intersection_diagnostics[
                    "outside intersection (%)"
                ].map(lambda value: f"{value:.1f}%")}
            ))

        combined_fit_records = [
            model for model in de_run_metadata.get("models", [])
            if model.get("model") in {"merged", "combined"}
        ]
        combined_covariate_summary = []
        for covariate, results in combined_results_by_covariate.items():
            padj = pd.to_numeric(results["padj"], errors="coerce")
            fit_record = next(
                (model for model in combined_fit_records
                 if covariate in model.get("results_by_covariate", {})),
                {},
            )
            combined_covariate_summary.append({
                "covariate": covariate,
                "design": fit_record.get("design", "—"),
                "FDR-significant genes": int((padj < de_alpha).sum()),
                "studies used": ", ".join(fit_record.get("studies", [])),
                "samples used": fit_record.get("n_samples", "—"),
            })
        if combined_covariate_summary:
            display(Markdown("### Combined covariate models"))
            display(pd.DataFrame(combined_covariate_summary))

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
                    "samples used": sum(sample_counts.values()) if sample_counts else "—",
                    "genes tested": subset["gene"].nunique(),
                    "FDR-significant genes": int(
                        (pd.to_numeric(subset["padj"], errors="coerce") < de_alpha).sum()
                    ),
                })
        if per_study_covariate_summary:
            display(Markdown("### Per-study covariate models"))
            display(pd.DataFrame(per_study_covariate_summary))

        display(Markdown(
            f"PyDESeq2 covariate effects; FDR threshold **{de_alpha:g}**. "
            "Age and sex share the study-adjusted model. BMI and CMV use separate "
            "models adjusted for age and sex, with study adjustment when multiple "
            "studies contribute. Missing covariate values are excluded from that "
            "covariate's fit. Categorical effects use female as the sex reference "
            "and no CMV (negative where that is the source label) as the CMV reference."
        ))
        if summary_rows:
            display(pd.DataFrame(summary_rows))

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
                top_genes = recurrence.head(20).sort_values("studies_associated")
                sns.barplot(
                    data=top_genes, x="studies_associated", y="gene", color="#4c72b0",
                    ax=axes[0],
                )
                axes[0].set_xlabel("Studies with FDR-significant age association")
                axes[0].set_ylabel("")
                axes[0].set_title("Genes recurring across studies (top 20)")
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
                display(recurrence.head(20))

        for covariate, combined_results in combined_results_by_covariate.items():
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
                volcano["association"] = np.where(
                    volcano["padj"] < de_alpha, "FDR significant", "Not significant"
                )
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
                n_samples_used = sum(sample_counts.values()) if sample_counts else "unavailable"
                display(Markdown(
                    f"**{covariate.upper()} contrast:** {contrast}. **Studies used:** "
                    f"{study_summary or 'unavailable in this legacy result set'}; "
                    f"**samples used:** {n_samples_used}."
                ))
                figure, axis = plt.subplots(figsize=(9, 6.5))
                sns.scatterplot(
                    data=volcano, x="log2FoldChange", y="minus_log10_padj",
                    hue="association", hue_order=["Not significant", "FDR significant"],
                    palette={"Not significant": "#b8b8b8", "FDR significant": "#c44e52"},
                    s=16, alpha=0.7, linewidth=0, ax=axis,
                )
                axis.axhline(-np.log10(de_alpha), color="#555555", linestyle="--", linewidth=1)
                axis.axvline(0, color="#555555", linewidth=0.8)
                if covariate == "age":
                    effect_label = "log2 fold change per year"
                elif covariate == "bmi":
                    effect_label = "log2 fold change per BMI unit"
                else:
                    effect_label = f"log2 fold change ({contrast})"
                axis.set_xlabel(f"{covariate.upper()} effect ({effect_label})")
                axis.set_ylabel("−log10(adjusted p-value)")
                if len(studies_used) > 1:
                    title = f"Combined study-adjusted {covariate.upper()} association"
                else:
                    study_name = studies_used[0] if studies_used else "single study"
                    title = f"{covariate.upper()} association in {study_name}"
                axis.set_title(title)
                labels = volcano.loc[volcano["padj"] < de_alpha].nsmallest(12, "padj")
                for _, row in labels.iterrows():
                    axis.annotate(
                        str(row["gene"]),
                        (row["log2FoldChange"], row["minus_log10_padj"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7,
                    )
                figure.tight_layout()
                plt.show()
        display(Markdown(
            "_DE interpretation: the recurrence panels emphasize agreement among available per-study "
            "models; each combined volcano summarizes all eligible samples with age and sex adjustment. "
            "When only one study contributes, study adjustment is omitted and the plot describes that "
            "study alone. Combined fits do not by themselves establish replication across studies._"
        ))

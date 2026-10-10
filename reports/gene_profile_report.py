# %%
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import HTML, Markdown, display
from matplotlib.lines import Line2D

from pbmc_pipeline.report_theme import render_report_header

sns.set_theme(style="whitegrid")
profile_dir = Path(os.environ["GENE_PROFILE_DIR"])
summary = json.loads((profile_dir / "fit_summary.json").read_text())
cell_type = summary["cell_type"]
gene = summary["gene"]
gene_result = pd.read_csv(profile_dir / "gene_age_result.csv")
gene_samples = pd.read_csv(profile_dir / "gene_sample_counts.csv")
model_fit_diagnostics = pd.read_csv(profile_dir / "model_fit_diagnostics.csv")
mean_variance = pd.read_csv(profile_dir / "mean_variance.csv")
residual_correlations = pd.read_csv(profile_dir / "gene_residual_correlations.csv")
top_gene_correlations = pd.read_csv(
    profile_dir / "top_gene_residual_correlations.csv",
)
top_correlated_residuals = pd.read_csv(
    profile_dir / "top_correlated_gene_residuals.csv",
)
age_bin_samples = pd.read_csv(profile_dir / "naive_age_bin_samples.csv")
age_bin_model_result = pd.read_csv(profile_dir / "naive_age_bin_model_result.csv")
age_bin_predictions = pd.read_csv(profile_dir / "naive_age_bin_predictions.csv")
study_slope_model_result = pd.read_csv(
    profile_dir / "naive_study_slope_model_result.csv",
)
study_slopes = pd.read_csv(profile_dir / "naive_study_slopes.csv")
study_slope_predictions = pd.read_csv(
    profile_dir / "naive_study_slope_predictions.csv",
)
mixed_model_fixed_effects = pd.read_csv(profile_dir / "naive_mixed_model_fixed_effects.csv")
mixed_model_study_slopes = pd.read_csv(profile_dir / "naive_mixed_model_study_slopes.csv")
mixed_model_predictions = pd.read_csv(profile_dir / "naive_mixed_model_predictions.csv")

display(HTML(render_report_header(
    title=f"{gene} · {cell_type}",
    eyebrow="PBMC ageing · focused gene profile",
    subtitle=(
        "One refit of the configured merged DESeq2 model, followed by focused "
        "expression and dispersion diagnostics for the selected gene."
    ),
    metrics=[
        ("Cell type", cell_type),
        ("Gene", gene),
        ("Samples", f"{summary['n_samples']:,}"),
        ("Studies", str(summary["n_studies"])),
        ("Fit genes", f"{summary['n_genes_fit']:,}"),
    ],
    title_tag="div",
)))


def display_fit_diagnostics(model_name):
    rows = model_fit_diagnostics.loc[model_fit_diagnostics["model"] == model_name]
    if rows.empty:
        raise ValueError(f"Fit diagnostics are missing for {model_name!r}")
    row = rows.iloc[0]
    residual_sd = (
        f"{row['residual_sd']:.4g} ({row['residual_sd_scale']})"
        if pd.notna(row["residual_sd"]) else "not available"
    )
    r_squared = (
        f"{row['r_squared']:.3f} ({row['r_squared']:.1%})"
        if pd.notna(row["r_squared"]) else "not available"
    )
    summary = (
        f"**Residual SD:** {residual_sd}  \n"
        f"**{row['r_squared_label']}:** {r_squared}"
    )
    if pd.notna(row["fixed_effects_r_squared"]):
        summary += (
            f"  \n**Fixed-effects-only R²:** "
            f"{row['fixed_effects_r_squared']:.3f} "
            f"({row['fixed_effects_r_squared']:.1%})"
        )
    display(Markdown(summary))

# %% [markdown]
# # PyDESeq2
#
# The sections below describe the configured merged DESeq2 fit for this cell
# type, with diagnostics focused on the selected gene.

# %% [markdown]
# ## Model and gene result
#
# The result below is the age coefficient from the configured merged model.
# Samples follow the production DE inclusion rules, and genes are restricted
# to those available in every study represented in this fit.

# %%
display(Markdown(f"**Design:** `{summary['design']}`"))
display_fit_diagnostics("PyDESeq2")
display(Markdown(
    "R² values are in-sample summaries, not cross-validated predictive scores. "
    "PyDESeq2 R² is calculated on raw counts; the naive models use "
    "`log(normalized count + 1)`. The mixed model reports fixed-effects-only "
    "and conditional R², with the latter including study-specific slope estimates."
))
display(gene_result)
display(Markdown(
    f"Included studies: {', '.join(summary['studies'])}. "
))
display(pd.Series(summary["metadata_exclusions"], name="pseudobulks").to_frame())

# %% [markdown]
# ## Mean–dispersion position
#
# Each point shows a gene's normalized mean and final dispersion estimate used
# by the fitted model. PyDESeq2 shrinks gene-wise estimates toward the fitted
# trend for most genes; flagged dispersion outliers retain their gene-wise
# estimate. For the selected gene, the open square marks the gene-wise estimate
# before shrinkage and the filled point marks its final dispersion. The orange
# curve is the fitted trend, not the selected gene's final dispersion value.

# %%
mean_column = "mean_normalized_count"
dispersion_column = "dispersions"
means = pd.to_numeric(mean_variance[mean_column], errors="coerce")
dispersions = pd.to_numeric(mean_variance[dispersion_column], errors="coerce")
valid = np.isfinite(means) & (means > 0) & np.isfinite(dispersions) & (dispersions > 0)
figure, axis = plt.subplots(figsize=(8.5, 5.7))
axis.scatter(
    means[valid], dispersions[valid], s=10, alpha=0.35,
    color="#777777", linewidths=0, label="Final dispersions used in model",
)
if "fitted_dispersions" in mean_variance:
    fitted = pd.to_numeric(mean_variance["fitted_dispersions"], errors="coerce")
    fitted_valid = valid & np.isfinite(fitted) & (fitted > 0)
    order = np.argsort(means[fitted_valid].to_numpy())
    axis.plot(
        means[fitted_valid].to_numpy()[order], fitted[fitted_valid].to_numpy()[order],
        color="#d55e00", linewidth=1.5, label="Fitted dispersion trend",
    )
target = mean_variance.loc[mean_variance["gene"].astype(str) == gene]
target_mean = float(pd.to_numeric(target[mean_column], errors="coerce").iloc[0])
target_dispersion = float(pd.to_numeric(target[dispersion_column], errors="coerce").iloc[0])
target_genewise_dispersion = float(
    pd.to_numeric(target["genewise_dispersions"], errors="coerce").iloc[0]
)
if np.isfinite(target_mean) and target_mean > 0 and np.isfinite(target_dispersion) and target_dispersion > 0:
    axis.scatter(
        [target_mean], [target_dispersion], s=75, color="#0072b2",
        edgecolor="white", zorder=4,
    )
    axis.annotate(gene, (target_mean, target_dispersion), xytext=(6, 5), textcoords="offset points")
    if np.isfinite(target_genewise_dispersion) and target_genewise_dispersion > 0:
        axis.plot(
            [target_mean, target_mean],
            [target_genewise_dispersion, target_dispersion],
            color="#0072b2", linewidth=1.2, linestyle=":", zorder=3,
        )
        axis.scatter(
            [target_mean], [target_genewise_dispersion],
            s=72, marker="s", facecolors="none", edgecolors="#d55e00",
            linewidths=1.7, zorder=5,
            label="Target gene-wise dispersion before shrinkage",
        )
axis.set_xscale("log")
axis.set_yscale("log")
axis.set_xlabel("Mean normalized count")
axis.set_ylabel("Estimated dispersion")
axis.grid(True, which="both", alpha=0.2)
axis.legend()
figure.tight_layout()
plt.show()

# %% [markdown]
# ## Pearson residuals
#
# Each point is the selected gene's Pearson residual from the fitted merged
# model. Residuals use the observed raw count, the model's fitted mean, and the
# final gene dispersion:
# `(count - fitted mean) / sqrt(fitted mean + dispersion * fitted mean²)`.
# Both plots split panels by study and color points by study site. The model
# already adjusts for the other terms in its design.

# %%
study_palette = dict(zip(
    summary["studies"],
    sns.color_palette("colorblind", n_colors=len(summary["studies"])),
    strict=True,
))

def make_residual_plot_config(metadata, residual_label):
    metadata = metadata[["study", "study_site"]].copy()
    studies = metadata["study"].astype(str)
    sites = metadata["study_site"].astype(str)
    study_order = [study for study in summary["studies"] if study in set(studies)]
    site_order = sorted(sites.unique())
    site_palette = dict(zip(
        site_order,
        sns.color_palette("colorblind", n_colors=len(site_order)),
        strict=True,
    ))
    return {
        "metadata": metadata,
        "facet_column": "study",
        "facet_order": study_order,
        "color_column": "study_site",
        "color_order": site_order,
        "palette": site_palette,
        "facet_wrap": int(np.ceil(np.sqrt(max(len(study_order), 1)))),
        "legend_title": "Study site",
        "residual_label": residual_label,
    }


# Keep each fitted model's metadata and residual plot choices together so the
# shared plot helper can be reused without misaligning samples across fits.
residual_plot_configs = {
    "PyDESeq2": make_residual_plot_config(
        gene_samples, "Pearson residual",
    ),
    "Naive linear age model": make_residual_plot_config(
        gene_samples, "OLS residual",
    ),
    "Naive age-bin model": make_residual_plot_config(
        age_bin_samples, "Age-bin OLS residual",
    ),
    "Naive age model with study-specific slopes": make_residual_plot_config(
        gene_samples, "Study-slope OLS residual",
    ),
    "Naive mixed age model": make_residual_plot_config(
        gene_samples, "Mixed-model residual",
    ),
}
active_residual_plot_config = residual_plot_configs["PyDESeq2"]


def plot_residuals(resids, x_values, x_label):
    metadata = active_residual_plot_config["metadata"]
    if len(resids) != len(metadata) or len(x_values) != len(metadata):
        raise ValueError("Residuals and x values must align with the model sample metadata")
    plot_data = metadata.copy()
    plot_data["residual"] = np.asarray(resids, dtype=float)
    plot_data["x_value"] = np.asarray(x_values, dtype=float)
    residual_grid = sns.relplot(
        data=plot_data,
        x="x_value",
        y="residual",
        hue=active_residual_plot_config["color_column"],
        hue_order=active_residual_plot_config["color_order"],
        palette=active_residual_plot_config["palette"],
        col=active_residual_plot_config["facet_column"],
        col_order=active_residual_plot_config["facet_order"],
        col_wrap=active_residual_plot_config["facet_wrap"],
        kind="scatter",
        legend=False,
        height=4.2,
        aspect=1.15,
        s=36,
        alpha=0.8,
        facet_kws={"sharex": True, "sharey": True},
    )
    residual_grid.set_axis_labels(x_label, active_residual_plot_config["residual_label"])
    residual_grid.set_titles("Study: {col_name}")
    for axis in residual_grid.axes.flat:
        axis.axhline(0, color="#555555", linewidth=1, linestyle="--")
        axis.grid(True, alpha=0.2)
    residual_grid.tight_layout(rect=(0, 0, 0.74, 1))
    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=active_residual_plot_config["palette"][site],
            markeredgecolor="none", label=site,
        )
        for site in active_residual_plot_config["color_order"]
    ]
    residual_grid.figure.legend(
        handles=legend_handles,
        labels=active_residual_plot_config["color_order"],
        title=active_residual_plot_config["legend_title"],
        bbox_to_anchor=(0.77, 0.5),
        loc="center left",
        frameon=False,
    )
    plt.show()

# %% [markdown]
# ### Against Age
#
# Age is included as a continuous model term in this fit.

# %%
plot_residuals(
    gene_samples["pearson_residual"], gene_samples["age"], "Age (years)",
)

# %% [markdown]
# ### Against Depth
#
# Depth is shown as `log10(total_counts)`, the library-depth covariate used in
# the fitted model.

# %%
plot_residuals(
    gene_samples["pearson_residual"],
    gene_samples["log10_total_counts"],
    "log10(total_counts)",
)

# %% [markdown]
# ### Correlated genes
#
# The table starts with the selected gene itself (r = 1), followed by its ten
# strongest residual correlates ranked by absolute Pearson correlation across
# included samples. It also shows each gene's mean normalized count and final
# dispersion. The scatter plots show the four strongest other genes; the
# signed correlation is shown in each title. These are descriptive residual
# associations from this fitted model.

# %%
plot_correlations = residual_correlations.head(4).copy()
if top_gene_correlations.empty:
    display(Markdown("_No other genes had finite residual correlations._"))
else:
    display(top_gene_correlations[[
        "gene", "pearson_r", "mean_normalized_count", "dispersions",
    ]].round({
        "pearson_r": 3,
        "mean_normalized_count": 2,
        "dispersions": 4,
    }))
    figure, axes = plt.subplots(
        2, 2, figsize=(12, 9), sharex=True, sharey=False,
    )
    axes = axes.ravel()
    for index, (axis, correlated) in enumerate(
        zip(axes, plot_correlations.itertuples(index=False), strict=False),
    ):
        rows = top_correlated_residuals.loc[
            top_correlated_residuals["correlated_gene"] == correlated.gene
        ]
        sns.scatterplot(
            data=rows,
            x="target_gene_pearson_residual",
            y="correlated_gene_pearson_residual",
            hue="study",
            hue_order=summary["studies"],
            palette=study_palette,
            s=36,
            alpha=0.8,
            legend="auto" if index == 0 else False,
            ax=axis,
        )
        axis.axhline(0, color="#777777", linewidth=0.8, linestyle="--")
        axis.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
        axis.set_title(f"{correlated.gene} (r={correlated.pearson_r:.2f})")
        axis.set_xlabel(f"{gene} Pearson residual")
        axis.set_ylabel("Correlated gene Pearson residual")
        axis.grid(True, alpha=0.2)
    for axis in axes[len(plot_correlations):]:
        axis.set_visible(False)
    figure.tight_layout()
    plt.show()

# %% [markdown]
# ## Sample-level normalized counts
#
# Points show DESeq2 size-factor normalized counts on a log scale, colored by
# study. This is a descriptive view of the selected gene across included
# pseudobulk samples.

# %%
figure, axis = plt.subplots(figsize=(8.5, 5.7))
for study, rows in gene_samples.groupby("study", observed=True, sort=True):
    axis.scatter(
        rows["age"], np.log10(rows["normalized_count"] + 1),
        s=36, alpha=0.3, label=str(study),
    )
axis.set_xlabel("Age (years)")
axis.set_ylabel("log10(normalized count + 1)")
axis.legend(title="Study", bbox_to_anchor=(1.02, 1), loc="upper left")
axis.grid(True, alpha=0.2)
figure.tight_layout()
plt.show()
display(Markdown(
    "Download [sample-level raw and normalized counts](gene_sample_counts.csv), "
    "including the per-sample VST value."
))

# %% [markdown]
# ## VST matrix for downstream modeling
#
# The compressed matrix contains all genes in the merged fit gene set, with one
# row per included pseudobulk. VST values are transformed data for downstream
# exploration or modeling; use raw integer counts for DESeq2 differential-
# expression inference.

# %%
display(Markdown(
    "Download the [full compressed VST matrix](vst_normalized_counts.csv.gz) "
    "and [included sample metadata](fit_sample_metadata.csv). The VST was fit "
    "with the selected cell type's DESeq2 model."
))

# %% [markdown]
# # Naive normalization
#
# All naive models use total-count size factors followed by `log(normalized +
# 1)`. The age-bin fit calculates its geometric-mean reference across all
# cell-type samples before applying the age-bin model's sample filters, so
# samples at lower ages contribute to its normalization scale.
# %% [markdown]
# ## Normalization procedure
#
# For each sample, divide raw counts by `total_counts / geometric_mean`. The
# continuous-age models use the geometric mean among their included samples;
# the age-bin model calculates its geometric mean across all samples in this
# cell type first, then applies model eligibility filters. Thus samples below
# age 20 affect the age-bin scale even if any samples are later excluded from
# the fit. The pseudocount is added after normalization, so `1` is one unit on
# the geometric-mean library scale.
# %%
display(Markdown(
    f"**Geometric mean for the continuous-age models:** "
    f"{summary['geometric_mean_total_counts']:.2e} counts  \n"
    f"**Geometric mean for the age-bin model:** "
    f"{summary['age_bin_model']['geometric_mean_total_counts']:.2e} counts "
    f"(across all {summary['age_bin_model']['normalization_n_samples']:,} "
    "cell-type samples before model filters; includes lower-age samples)"
))

# %% [markdown]
# ### Distribution of normalized counts
#
# The histogram shows target-gene counts after normalization and before the
# log transform and pseudocount are applied.
# %%
figure, axis = plt.subplots(figsize=(8.5, 5.0))
axis.hist(
    gene_samples["naive_normalized_count"], bins="auto",
    color="#4c78a8", edgecolor="white",
)
axis.set_xlabel("Naive normalized count")
axis.set_ylabel("Number of samples")
axis.set_title(f"{gene}: normalized counts")
axis.grid(True, axis="y", alpha=0.2)
figure.tight_layout()
plt.show()

# %% [markdown]
# ### Naive normalized expression versus VST
#
# Each point is the selected gene in one sample. The x axis shows naive
# normalized log expression and the y axis shows PyDESeq2 VST expression.
# Color indicates sample library depth.
# %%
figure, axis = plt.subplots(figsize=(8.5, 6.0))
scatter = axis.scatter(
    gene_samples["naive_log_expression"], gene_samples["vst_count"],
    c=gene_samples["log10_total_counts"], cmap="viridis", s=42, alpha=0.3,
    edgecolors="none",
)
axis.set_xlabel("Naive log(normalized count + 1)")
axis.set_ylabel("PyDESeq2 VST expression")
axis.grid(True, alpha=0.2)
colorbar = figure.colorbar(scatter, ax=axis)
colorbar.set_label("log10(total_counts)")
figure.tight_layout()
plt.show()

# %% [markdown]
# # Naive linear age model
#
# OLS of the selected gene's normalized log expression against the same
# covariates as the PyDESeq2 model.
# %% [markdown]
# ## Model and coefficient results
#
# This fit uses the continuous-age model's included samples and configured
# merged design.
# %%
naive_model_result = pd.read_csv(profile_dir / "naive_linear_model_result.csv")
display(Markdown(f"**Design:** `{summary['design']}`"))
display_fit_diagnostics("Naive linear age model")
display(naive_model_result.round({
    "coefficient": 4, "standard_error": 4, "t_value": 3,
    "p_value": 4, "ci_lower": 4, "ci_upper": 4,
}))

# %% [markdown]
# ## Residuals
#
# OLS residuals are observed minus fitted normalized log expression. Panels
# split by study and points are colored by study site.
# %%
active_residual_plot_config = residual_plot_configs["Naive linear age model"]

# %% [markdown]
# ### Against Age
# %%
plot_residuals(
    gene_samples["naive_residual"], gene_samples["age"], "Age (years)",
)

# %% [markdown]
# ### Against Depth
# %%
plot_residuals(
    gene_samples["naive_residual"],
    gene_samples["log10_total_counts"],
    "log10(total_counts)",
)

# %% [markdown]
# # Naive age-bin model
#
# OLS with configured-width age bins as categorical levels and the other
# varying terms from the merged DE design. This fit includes all valid ages,
# including samples under 20.
# %% [markdown]
# ## Model and gene result
# %%
display(Markdown(
    f"**Design:** `{summary['age_bin_model']['design']}`  \n"
    f"**Samples:** {summary['age_bin_model']['n_samples']:,}  \n"
    f"**Bin width:** {summary['age_bin_model']['bin_width_years']} years  \n"
    f"**Reference bin:** {summary['age_bin_model']['reference_bin']}"
))
display_fit_diagnostics("Naive age-bin model")
display(age_bin_model_result.round({
    "coefficient": 4, "standard_error": 4, "t_value": 3,
    "p_value": 4, "ci_lower": 4, "ci_upper": 4,
}))

# %% [markdown]
# ## Residuals
#
# Age-bin OLS residuals are observed minus fitted normalized log expression.
# Each panel is a study, with points colored by study site.
# %%
active_residual_plot_config = residual_plot_configs["Naive age-bin model"]

# %% [markdown]
# ### Against Age
# %%
plot_residuals(age_bin_samples["residual"], age_bin_samples["age"], "Age (years)")

# %% [markdown]
# ### Against Depth
# %%
plot_residuals(
    age_bin_samples["residual"], age_bin_samples["log10_total_counts"],
    "log10(total_counts)",
)

# %% [markdown]
# ## Sample-level expression and age-bin levels
#
# Points show sample expression, colored by study. Horizontal segments and
# confidence bars show age-bin fitted values standardized over the observed
# non-age covariate distribution.
# %%
figure, axis = plt.subplots(figsize=(10, 6.0))
for study, rows in age_bin_samples.groupby("study", observed=True, sort=True):
    axis.scatter(
        rows["age"], rows["naive_log_expression"],
        color=study_palette.get(str(study)), s=32, alpha=0.3, label=str(study),
    )
for prediction in age_bin_predictions.itertuples(index=False):
    midpoint = (prediction.age_bin_start + prediction.age_bin_end) / 2
    half_width = summary["age_bin_model"]["bin_width_years"] * 0.38
    axis.hlines(
        prediction.fitted_expression, midpoint - half_width, midpoint + half_width,
        color="#222222", linewidth=2.4, zorder=4,
    )
    axis.vlines(
        midpoint, prediction.ci_lower, prediction.ci_upper,
        color="#222222", linewidth=1.2, zorder=4,
    )
axis.set_xlabel("Age (years)")
axis.set_ylabel("log(naive normalized count + 1)")
axis.legend(title="Study", bbox_to_anchor=(1.02, 1), loc="upper left")
axis.grid(True, alpha=0.2)
figure.tight_layout()
plt.show()
display(age_bin_predictions.round({
    "fitted_expression": 3, "ci_lower": 3, "ci_upper": 3,
}))

# %% [markdown]
# # Naive age model with study-specific slopes
#
# This linear model allows a different age slope and intercept for each study.
# Study site is omitted because sites are nested within studies; sex and library
# depth remain as covariates.
# %% [markdown]
# ## Per-study age slopes
# %%
display(Markdown(f"**Design:** `{summary['study_slope_model']['design']}`"))
display_fit_diagnostics("Naive age model with study-specific slopes")
display(study_slopes.round({
    "slope": 4, "standard_error": 4, "t_value": 3,
    "p_value": 4, "ci_lower": 4, "ci_upper": 4,
}))

# %% [markdown]
# ## Residuals
#
# These are OLS residuals from the model with separate age slopes and
# intercepts for each study. Panels are split by study and colored by site.
# %%
active_residual_plot_config = residual_plot_configs[
    "Naive age model with study-specific slopes"
]

# %% [markdown]
# ### Against Age
# %%
plot_residuals(
    gene_samples["study_slope_residual"], gene_samples["age"], "Age (years)",
)

# %% [markdown]
# ### Against Depth
# %%
plot_residuals(
    gene_samples["study_slope_residual"],
    gene_samples["log10_total_counts"], "log10(total_counts)",
)

# %% [markdown]
# ## Sample-level expression and fitted slopes
#
# Points show each sample, colored by study. Lines show fitted expression over
# the observed age range for each study, evaluated at the reference sex and
# mean log library depth. Slopes are on the log(normalized count + 1) scale.
# %%
figure, axis = plt.subplots(figsize=(10, 6.0))
for study, rows in gene_samples.groupby("study", observed=True, sort=True):
    color = study_palette.get(str(study))
    axis.scatter(
        rows["age"], rows["naive_log_expression"],
        color=color, s=32, alpha=0.3, label=str(study),
    )
    predictions = study_slope_predictions.loc[
        study_slope_predictions["study"].astype(str) == str(study)
    ]
    axis.plot(
        predictions["age"], predictions["fitted_expression"],
        color=color, linewidth=2.2,
    )
axis.set_xlabel("Age (years)")
axis.set_ylabel("log(naive normalized count + 1)")
axis.legend(title="Study", bbox_to_anchor=(1.02, 1), loc="upper left")
axis.grid(True, alpha=0.2)
figure.tight_layout()
plt.show()
display(study_slope_model_result.round({
    "coefficient": 4, "standard_error": 4, "t_value": 3,
    "p_value": 4, "ci_lower": 4, "ci_upper": 4,
}))
display(Markdown(
    "Sample-level normalized counts, expression, fitted values, and residuals "
    "are available in [gene_sample_counts.csv](gene_sample_counts.csv). The "
    "age-bin sample data and adjusted levels are in "
    "[naive_age_bin_samples.csv](naive_age_bin_samples.csv) and "
    "[naive_age_bin_predictions.csv](naive_age_bin_predictions.csv)."
))

# %% [markdown]
# # Naive mixed age model with study fixed intercepts and random slopes
#
# This model keeps study intercepts as fixed effects and partially pools the
# study-specific age slopes around a shared mean slope. The random-slope
# variance describes between-study slope variation. With few studies, treat
# this variance and the conditional per-study slopes as exploratory.
# %% [markdown]
# ## Fixed effects and fit diagnostics
#
# The age coefficient is the across-study mean slope. Study fixed effects are
# intercept contrasts at the mean sample age. Conditional study slopes add the
# estimated random-slope deviation to the shared slope; they are partially
# pooled estimates, not independent study-by-study fits.
# %%
mixed_diagnostics = summary["mixed_model"]
display_fit_diagnostics("Naive mixed age model")
display(Markdown(
    f"**Design:** `{mixed_diagnostics['design']}`  \n"
    f"**Samples / studies:** {mixed_diagnostics['n_samples']:,} / "
    f"{mixed_diagnostics['n_studies']}  \n"
    f"**Converged:** {mixed_diagnostics['converged']}  \n"
    f"**Random-slope SD:** {mixed_diagnostics['random_slope_sd_per_year']:.4g} "
    "log-expression units per year  \n"
    f"**Residual SD:** {mixed_diagnostics['residual_sd']:.4g}  \n"
    f"**Fixed-effect uncertainty:** {mixed_diagnostics['fixed_effect_uncertainty']}"
))
display(mixed_model_fixed_effects.round({
    "coefficient": 4, "standard_error": 4, "test_statistic": 3,
    "p_value": 4, "ci_lower": 4, "ci_upper": 4,
}))
display(mixed_model_study_slopes.round({
    "fixed_age_slope": 4,
    "random_slope_deviation": 4,
    "conditional_age_slope": 4,
}))

# %% [markdown]
# ## Residuals
#
# Residuals use the fitted fixed effects plus each study's estimated random
# slope deviation. Panels are split by study and colored by study site.
# %%
active_residual_plot_config = residual_plot_configs["Naive mixed age model"]

# %% [markdown]
# ### Against Age
# %%
plot_residuals(
    gene_samples["mixed_model_residual"], gene_samples["age"], "Age (years)",
)

# %% [markdown]
# ### Against Depth
# %%
plot_residuals(
    gene_samples["mixed_model_residual"],
    gene_samples["log10_total_counts"], "log10(total_counts)",
)

# %% [markdown]
# ## Sample-level expression and partially pooled slopes
#
# Points show each sample. Lines show the fitted study intercept and
# conditional slope at reference categorical covariates and mean library depth.
# %%
figure, axis = plt.subplots(figsize=(10, 6.0))
for study, rows in gene_samples.groupby("study", observed=True, sort=True):
    color = study_palette.get(str(study))
    axis.scatter(
        rows["age"], rows["naive_log_expression"],
        color=color, s=32, alpha=0.3, label=str(study),
    )
    predictions = mixed_model_predictions.loc[
        mixed_model_predictions["study"].astype(str) == str(study)
    ]
    axis.plot(
        predictions["age"], predictions["fitted_expression"],
        color=color, linewidth=2.2,
    )
axis.set_xlabel("Age (years)")
axis.set_ylabel("log(naive normalized count + 1)")
axis.legend(title="Study", bbox_to_anchor=(1.02, 1), loc="upper left")
axis.grid(True, alpha=0.2)
figure.tight_layout()
plt.show()

# %% [markdown]
# # limma-voom
#
# TODO: implement R execution.

# %% [markdown]
# # DESeq2
#
# TODO: implement R execution.

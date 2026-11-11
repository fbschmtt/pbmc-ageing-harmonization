"""Importable sample-level cluster-fraction figures and summaries."""
from __future__ import annotations

import base64
import re
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from IPython.display import Markdown, display


def natural_sort_key(value):
    """Sort labels naturally so numbered clusters follow numeric order."""
    return tuple(
        (1, int(part)) if part.isdigit() else (0, part.casefold())
        for part in re.split(r"(\d+)", str(value))
    )


def has_age_span(data, *, min_samples=2):
    return len(data) >= min_samples and data["age"].nunique() >= 2


def draw_cluster_fraction_axis(axis, subset, studies, palette):
    sns.scatterplot(
        data=subset, x="age", y="fraction_within_cell_type", hue="study", palette=palette,
        s=20, alpha=0.5, edgecolor="none", ax=axis, legend=False,
    )
    # study-specific linear trends
    for study in studies:
        study_data = subset.loc[subset["study"] == study]
        if has_age_span(study_data):
            sns.regplot(
                data=study_data, x="age", y="fraction_within_cell_type", scatter=False,
                ci=None, color=palette[study],
                line_kws={"linewidth": 1.8, "linestyle": "--"}, ax=axis,
            )
    if has_age_span(subset):
        sns.regplot(
            data=subset, x="age", y="fraction_within_cell_type", scatter=False, ci=None,
            color="#1f1f1f", line_kws={"linewidth": 2.6}, ax=axis,
        )
    axis.set_ylabel("Fraction within cell type")


def cluster_fraction_fit_summary(data, study):
    """Summarize plotted samples and an unadjusted fraction-per-decade slope."""
    x = data["age"].to_numpy(dtype=float)
    y = data["fraction_within_cell_type"].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if not len(y):
        return None
    slope = float(np.polyfit(x, y, 1)[0] * 1000) if len(x) >= 2 and np.unique(x).size >= 2 else None
    return {"study": study, "samples": len(y), "meanFraction": float(np.mean(y)), "slopePpPerDecade": slope}


def plot_cluster_fractions(data, *, palette, display_figure=None):
    """Show descriptive pooled and study-specific linear trends."""
    if data.empty:
        display(Markdown("_No samples meet the denominator cutoff for cluster fractions._"))
        return []
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
        draw_cluster_fraction_axis(axis, subset, studies, palette)
        axis.set_title(f"Cluster {cluster}")
    for axis in axes.flat[len(clusters):]:
        axis.remove()
    handles = [
        plt.Line2D([], [], marker="o", linestyle="--", color=palette[study], label=study, markersize=5)
        for study in studies
    ]
    handles.append(plt.Line2D([], [], color="#1f1f1f", linewidth=2.6, label="Pooled"))
    figure.legend(handles=handles, title="Study / trend", loc="upper center",
                  ncol=min(4, len(handles)), frameon=False)
    figure.suptitle("Cluster fractions: linear trends (study-unadjusted)", y=0.98)
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    if display_figure is not None:
        display_figure(figure, "Show cluster-composition curves")
    else:
        plt.close(figure)
    details = []
    for cluster in clusters:
        subset = plot_data.loc[plot_data["cluster"].astype(str) == cluster]
        present_studies = sorted(subset["study"].unique())
        detail_figure, axis = plt.subplots(figsize=(7, 4.5))
        draw_cluster_fraction_axis(axis, subset, present_studies, palette)
        axis.set_title(f"Cluster {cluster}")
        detail_handles = [
            plt.Line2D([], [], marker="o", linestyle="--", color=palette[study], label=study, markersize=5)
            for study in present_studies
        ]
        detail_handles.append(plt.Line2D([], [], color="#1f1f1f", linewidth=2.6, label="Pooled"))
        axis.legend(handles=detail_handles, title="Study / trend", fontsize=8, frameon=False)
        detail_figure.tight_layout()
        image_buffer = BytesIO()
        detail_figure.savefig(
            image_buffer, format="png", dpi=150, bbox_inches="tight", transparent=True,
        )
        plt.close(detail_figure)
        fits = [cluster_fraction_fit_summary(subset, None)]
        fits.extend(cluster_fraction_fit_summary(subset.loc[subset["study"] == study], study)
                    for study in present_studies)
        details.append({"cluster": cluster, "image": base64.b64encode(image_buffer.getvalue()).decode("ascii"),
                        "fits": [fit for fit in fits if fit is not None]})
    return details


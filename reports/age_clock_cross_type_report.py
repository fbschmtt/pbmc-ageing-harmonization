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
report_dir = Path(os.environ["AGE_CLOCK_CROSS_TYPE_DIR"])
summary = json.loads((report_dir / "summary.json").read_text())
methods = summary.get("methods", [])
if not methods:
    raise ValueError(f"{report_dir}: no cross-cell-type clock methods are available")

display(HTML(render_report_header(
    title="Combining clocks across cell types",
    eyebrow="Transcriptomic age clocks",
    subtitle=(
        "For each donor, this analysis averages predictions from the same clock "
        "method fitted separately in multiple cell types. It evaluates matched donors "
        "only and does not combine different clock algorithms."
    ),
    metrics=[
        ("Cell types", f"{len(summary.get('cell_types', [])):,}"),
        ("Shared clock methods", f"{len(methods):,}"),
        ("Minimum per donor", f"{summary.get('minimum_cell_types_per_donor', 2)} cell types"),
    ],
)))

performance_rows = []
component_rows = []
for method in methods:
    for split in ("training", "held_out"):
        result = method.get(split)
        if not result:
            continue
        metrics = result["metrics"]
        performance_rows.append({
            "Clock method": method.get("name", method["id"]),
            "Evaluation": "Training · leave-one-study-out" if split == "training" else "Held-out studies",
            "Donors": metrics["n_donors"],
            "Studies": metrics["n_studies"],
            "Mean cell types per donor": metrics["mean_cell_types_per_donor"],
            "MAE (years)": metrics["mae_years"],
            "RMSE (years)": metrics["rmse_years"],
            "R²": metrics["r_squared"],
            "Bias (years)": metrics["bias_years"],
        })
    for split, records in method.get("component_metrics", {}).items():
        for record in records:
            component_rows.append({
                "Clock method": method.get("name", method["id"]),
                "Evaluation": "Training · leave-one-study-out" if split == "training" else "Held-out studies",
                "Cell type": record["cell_type"],
                "Matched donors": record["n_donors"],
                "Cell-type clock MAE (years)": record["clock_mae_years"],
                "Mean clock MAE on same donors (years)": record["mean_ensemble_mae_years"],
                "Change in MAE (years)": record["mean_ensemble_mae_years"] - record["clock_mae_years"],
            })

display(Markdown(
    "## Combined prediction performance\n\n"
    "The ensemble and its component clocks are compared on donors with predictions "
    "from at least two cell types. MAE improvement is descriptive and is not an "
    "independent validation beyond the listed study split."
))
display(pd.DataFrame(performance_rows).round(3))
if component_rows:
    display(Markdown("## Component clocks on the same donors"))
    display(pd.DataFrame(component_rows).round(3))

correlation_path = report_dir / summary.get("residual_correlations_file", "residual_correlations.csv")
correlations = pd.read_csv(correlation_path) if correlation_path.is_file() else pd.DataFrame()
if not correlations.empty:
    display(Markdown(
        "## Residual correlation between cell-type clocks\n\n"
        "Each correlation is calculated across donors shared by the two cell types, "
        "using prediction residuals (predicted minus observed age). Positive values "
        "indicate errors that tend to move together; negative values indicate errors "
        "that tend to move in opposite directions."
    ))
    display(correlations.sort_values(["split", "method", "pearson_r"], na_position="last").round(3))
    for (split, method), group in correlations.groupby(["split", "method"], sort=True):
        labels = sorted(set(group["cell_type_a"]) | set(group["cell_type_b"]))
        matrix = pd.DataFrame(index=labels, columns=labels, dtype=float)
        for row in group.itertuples(index=False):
            if pd.notna(row.pearson_r):
                matrix.loc[row.cell_type_a, row.cell_type_b] = row.pearson_r
                matrix.loc[row.cell_type_b, row.cell_type_a] = row.pearson_r
        for label in labels:
            matrix.loc[label, label] = 1.0
        figure, axis = plt.subplots(figsize=(max(6, len(labels) * 0.8), max(5, len(labels) * 0.65)))
        sns.heatmap(
            matrix, vmin=-1, vmax=1, center=0, cmap="vlag", annot=True,
            fmt=".2f", square=True, linewidths=0.5, ax=axis,
            cbar_kws={"label": "Pearson correlation of residuals"},
        )
        axis.set_title(f"{method} · {split.replace('_', ' ')} residuals")
        axis.set_xlabel("Cell type")
        axis.set_ylabel("Cell type")
        figure.tight_layout()
        display_accessible_figure(
            figure,
            alt_text=f"Heatmap of pairwise residual correlations for {method} clocks in the {split} split.",
        )
else:
    display(Markdown("Pairwise residual correlations were not available for this release."))

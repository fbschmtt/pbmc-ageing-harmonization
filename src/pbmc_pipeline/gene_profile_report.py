"""Execute and render the focused gene-profile notebook report."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .notebook import render_notebook


def format_model_results(table, *, digits: int = 4):
    """Round coefficient output and preserve small p-values in scientific notation."""
    formatted = table.copy()
    numeric = formatted.select_dtypes(include=["number"]).columns
    significance_columns = ["p_value", "pvalue", "padj"]
    formatted[numeric.difference(significance_columns)] = formatted[
        numeric.difference(significance_columns)
    ].round(digits)
    for column in significance_columns:
        if column in formatted:
            formatted[column] = pd.to_numeric(
                formatted[column], errors="coerce"
            ).map(lambda value: f"{value:.2e}" if pd.notna(value) else "—")
    return formatted


def generate_gene_profile_report(
    *, profile_dir: Path, template_path: Path, project_root: Path,
) -> Path:
    """Render fitted gene-profile artifacts and build their heading-based TOC."""
    import json

    summary_path = profile_dir / "fit_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = json.loads(summary_path.read_text())
    return render_notebook(
        template_path=template_path, output_dir=profile_dir, project_root=project_root,
        environment={"GENE_PROFILE_DIR": str(profile_dir.resolve())},
        title=f"{summary['gene']} · {summary['cell_type']} | PBMC ageing",
        log_path=profile_dir / "report_render.log", minimum_heading_level=1,
    )

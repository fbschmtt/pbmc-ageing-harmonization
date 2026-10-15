"""Render a standalone per-cell-type pseudobulk DE report."""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd


def _table(frame: pd.DataFrame, columns: list[str], *, limit: int = 25) -> str:
    available = [column for column in columns if column in frame.columns]
    if frame.empty or not available:
        return "<p>No result rows are available for this model.</p>"
    view = frame.loc[:, available].head(limit).copy()
    return view.to_html(index=False, escape=True, border=0, classes="results")


def _format_models(metadata: dict) -> str:
    rows = []
    for model in metadata.get("models", []):
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(model.get('model', '')))}</td>"
            f"<td>{html.escape(str(model.get('study') or 'All studies'))}</td>"
            f"<td>{html.escape(str(model.get('status', 'unknown')))}</td>"
            f"<td>{html.escape(str(model.get('n_samples', '')))}</td>"
            f"<td>{html.escape(str(model.get('n_genes', '')))}</td>"
            f"<td>{html.escape(str(model.get('n_significant_adjusted_p_lt_alpha', '')))}</td>"
            "</tr>"
        )
    if not rows:
        return "<p>No model metadata was recorded.</p>"
    return (
        "<table><thead><tr><th>Model</th><th>Study</th><th>Status</th>"
        "<th>Samples</th><th>Genes</th><th>Adjusted p &lt; alpha</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def render_differential_expression_report(result_dir: Path, output_dir: Path) -> Path:
    result_dir = result_dir.resolve()
    metadata_path = result_dir / "run_metadata.json"
    result_record_path = result_dir / "cell_type_result.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
    result_record = json.loads(result_record_path.read_text()) if result_record_path.is_file() else {}
    cell_type = str(metadata.get("cell_type") or result_record.get("cell_type") or result_dir.parent.name)
    warning = metadata.get("interpretation_warning") or ""
    warning_html = (
        f'<aside class="warning">{html.escape(str(warning))}</aside>' if warning else ""
    )
    result_sections = []
    result_files = sorted(
        path for path in result_dir.rglob("*.csv")
        if path.name not in {
            "age_model_diagnostics.csv",
            "age_bin_pearson_residuals_samples.csv",
        }
    )
    for path in result_files:
        frame = pd.read_csv(path)
        groups = (
            list(frame.groupby("covariate", sort=True, dropna=False))
            if path.name == "merged.csv" and "covariate" in frame
            else [(None, frame)]
        )
        for covariate, group in groups:
            sort_column = next(
                (column for column in ("padj", "omnibus_padj") if column in group), None
            )
            if sort_column:
                group = group.copy()
                group[sort_column] = pd.to_numeric(group[sort_column], errors="coerce")
                group = group.sort_values(sort_column, na_position="last")
            label = path.relative_to(result_dir).with_suffix("").as_posix().replace("_", " ")
            if covariate is not None:
                label = f"{label} — {covariate}"
            result_sections.append(
                f"<h3>{html.escape(label.title())}</h3>"
                + _table(group, [
                    "gene", "covariate", "contrast", "age_bin", "baseMean",
                    "trajectory_scaled", "omnibus_padj", "log2FoldChange", "lfcSE",
                    "pvalue", "padj",
                ])
            )
    results_html = "".join(result_sections) or "<p>No result tables were published.</p>"
    mode = html.escape(str(metadata.get("analysis_mode", "unknown")))
    samples = html.escape(str(result_record.get("n_pseudobulks", "unknown")))
    output_path = output_dir / "report.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ | PBMC ageing</title><style>
body{font:16px/1.5 system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#183044}
h1,h2{color:#075b78}table{border-collapse:collapse;width:100%;font-size:.9rem;display:block;overflow:auto}
th,td{padding:.45rem .6rem;border-bottom:1px solid #d8e1e8;text-align:left;white-space:nowrap}
thead{background:#eef5f8}.warning{padding:1rem;background:#fff4d6;border-left:4px solid #d69e00}
</style></head><body><h1>__TITLE__</h1>
<p>Analysis mode: <strong>__MODE__</strong>; pseudobulk samples: <strong>__SAMPLES__</strong>.</p>
__WARNING__<h2>Model summary</h2>__MODELS__
<h2>Differential-expression results</h2><p>Each table shows up to 25 rows, ordered by adjusted p-value where available.</p>__RESULT_TABLES__</body></html>"""
        .replace("__TITLE__", html.escape(cell_type))
        .replace("__MODE__", mode)
        .replace("__SAMPLES__", samples)
        .replace("__WARNING__", warning_html)
        .replace("__MODELS__", _format_models(metadata))
        .replace("__RESULT_TABLES__", results_html),
        encoding="utf-8",
    )
    return output_path


def render_per_type_reports(input_dir: Path, output_root: Path) -> list[Path]:
    rendered = []
    for result_dir in sorted(input_dir.glob("*/differential_expression")):
        if not (result_dir / "run_metadata.json").is_file():
            continue
        output_dir = output_root / result_dir.parent.name / "differential_expression"
        output_dir.mkdir(parents=True, exist_ok=True)
        rendered.append(render_differential_expression_report(result_dir, output_dir))
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    reports = render_per_type_reports(args.input_dir, args.output_root)
    print(f"Rendered {len(reports)} per-type DE reports")


if __name__ == "__main__":
    main()

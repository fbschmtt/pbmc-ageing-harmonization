"""Render a small HTML index for the published analysis reports."""
from __future__ import annotations

import argparse
import json
from html import escape
from html.parser import HTMLParser
from pathlib import Path


class _TitleParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_title = False
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "title":
            self.in_title = True

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "title":
            self.in_title = False

    def handle_data(self, data: str) -> None:
        if self.in_title:
            self.parts.append(data)


def _report_title(path: Path) -> str:
    parser = _TitleParser()
    try:
        with path.open(encoding="utf-8", errors="replace") as report:
            parser.feed(report.read(16_384))
    except OSError:
        pass
    title = "".join(parser.parts).strip()
    if title:
        return title.split("|", maxsplit=1)[0].strip()
    return path.parent.parent.name.replace("-", " ").title()


def render_report_index(outdir: Path) -> Path:
    outdir = outdir.resolve()
    analysis_dir = outdir / "analysis"
    harmonization_dir = outdir / "harmonization"
    qc_root = harmonization_dir / "harmonization_qc"
    merged_qc_report = qc_root / "merged" / "report.html"
    study_qc_reports = sorted(
        (path for path in qc_root.glob("*/report.html") if path.parent.name != "merged"),
        key=lambda path: path.parent.name,
    )
    benchmark_dir = harmonization_dir / "integration_benchmark"
    benchmark_report = benchmark_dir / "report.html"
    trajectory_dir = analysis_dir / "combined" / "trajectory"
    trajectory_report = trajectory_dir / "report.html"
    per_type_root = analysis_dir / "per_type"
    cell_type_manifest = per_type_root / "cell_type_manifest.json"
    if cell_type_manifest.is_file():
        document = json.loads(cell_type_manifest.read_text(encoding="utf-8"))
        manifest_slugs = [
            str(record["slug"])
            for record in document.get("cell_types", [])
            if isinstance(record, dict) and record.get("slug")
        ]
    else:
        manifest_slugs = []
    available_slugs = {
        path.parent.parent.name
        for pattern in ("*/analysis/report.html", "*/differential_expression/report.html")
        for path in per_type_root.glob(pattern)
    }
    slugs = sorted(set(manifest_slugs) | available_slugs)

    def link(path: Path, label: str) -> str:
        relative = path.relative_to(outdir).as_posix()
        return f'<a href="{escape(relative, quote=True)}">{escape(label)}</a>'

    harmonization_links = []
    if merged_qc_report.is_file():
        harmonization_links.append(
            f"<li>{link(merged_qc_report, 'Merged harmonization QC report')}</li>"
        )
    study_qc_table = ""
    if study_qc_reports:
        rows = []
        for report in study_qc_reports:
            study = _report_title(report)
            rows.append(
                "<tr>"
                f"<th scope=\"row\">{escape(study)}</th>"
                f"<td>{link(report, 'Open report')}</td>"
                f"<td><code>{escape(report.parent.name)}</code></td>"
                "</tr>"
            )
        study_qc_table = (
            "<h3>Per-study harmonization QC</h3>"
            "<table><thead><tr><th>Study</th><th>Report</th><th>Output folder</th>"
            "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
        )
    if harmonization_links or study_qc_table:
        harmonization_content = (
            ("<ul>" + "".join(harmonization_links) + "</ul>" if harmonization_links else "")
            + study_qc_table
        )
    else:
        harmonization_content = "<p>No harmonization QC reports are published yet.</p>"

    benchmark_links = []
    if benchmark_report.is_file():
        benchmark_links.append(link(benchmark_report, "Integration benchmark report (HTML)"))
    if benchmark_links:
        benchmark_content = "<ul>" + "".join(
            f"<li>{item}</li>" for item in benchmark_links
        ) + "</ul>"
    else:
        benchmark_content = "<p>No integration benchmark report is published yet.</p>"

    trajectory_links = []
    if trajectory_report.is_file():
        trajectory_links.append(link(trajectory_report, "Cross-cell-type analysis report (HTML)"))
    if not trajectory_links:
        trajectory_content = "<p>No cross-cell-type analysis report is published yet.</p>"
    else:
        trajectory_content = "<ul>" + "".join(
            f"<li>{item}</li>" for item in trajectory_links
        ) + "</ul>"

    if slugs:
        rows = []
        for slug in slugs:
            residual_report = per_type_root / slug / "analysis" / "report.html"
            de_report = per_type_root / slug / "differential_expression" / "report.html"
            report = residual_report if residual_report.is_file() else de_report
            cell_type = _report_title(report)
            residual_link = link(residual_report, "Open report") if residual_report.is_file() else "—"
            de_link = link(de_report, "Open report") if de_report.is_file() else "—"
            rows.append(
                "<tr>"
                f"<th scope=\"row\">{escape(cell_type)}</th>"
                f"<td>{residual_link}</td>"
                f"<td>{de_link}</td>"
                f"<td><code>{escape(slug)}</code></td>"
                "</tr>"
            )
        per_type_content = (
            "<table><thead><tr><th>Cell type</th><th>Residual single-cell report</th>"
            "<th>Differential-expression report</th><th>Output folder</th>"
            "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
        )
    else:
        per_type_content = "<p>No per-type reports are published yet.</p>"

    html = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>PBMC ageing analysis reports</title>
  <style>
    :root { color-scheme: light; font: 16px/1.5 system-ui, sans-serif; }
    body { max-width: 960px; margin: 3rem auto; padding: 0 1.25rem; color: #183044; }
    h1, h2 { color: #075b78; }
    section { margin: 2rem 0; }
    table { width: 100%; border-collapse: collapse; }
    th, td { padding: .65rem .8rem; border-bottom: 1px solid #d8e1e8; text-align: left; }
    thead { background: #eef5f8; }
    a { color: #075baf; }
    code { color: #475569; }
  </style>
</head>
<body>
  <h1>PBMC ageing analysis reports</h1>
  <p>Links are relative to this file and point to reports published under this output directory.</p>
  <section>
    <h2>Harmonization and merge</h2>
    __HARMONIZATION_REPORTS__
  </section>
  <section>
    <h2>Integration benchmark</h2>
    __BENCHMARK_REPORTS__
  </section>
  <section>
    <h2>Cross-cell-type analysis</h2>
    __TRAJECTORY_LINKS__
  </section>
  <section>
    <h2>Per-cell-type reports</h2>
    __CELL_TYPE_REPORTS__
  </section>
</body>
</html>
""".replace("__HARMONIZATION_REPORTS__", harmonization_content).replace(
        "__BENCHMARK_REPORTS__", benchmark_content
    ).replace("__TRAJECTORY_LINKS__", trajectory_content).replace(
        "__CELL_TYPE_REPORTS__", per_type_content
    )
    index_path = outdir / "index.html"
    outdir.mkdir(parents=True, exist_ok=True)
    index_path.write_text(html, encoding="utf-8")
    return index_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("output"))
    args = parser.parse_args()
    print(render_report_index(args.outdir.resolve()))


if __name__ == "__main__":
    main()

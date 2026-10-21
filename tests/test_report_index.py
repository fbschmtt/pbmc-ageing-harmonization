from html import unescape
from pathlib import Path

from scripts.render_report_index import render_report_index


def _write_report(path: Path, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"<html><head><title>{title}</title></head></html>", encoding="utf-8")


def test_report_index_links_published_reports_from_output_root(tmp_path: Path) -> None:
    _write_report(
        tmp_path / "harmonization/harmonization_qc/merged/report.html",
        "Merged QC | PBMC",
    )
    _write_report(
        tmp_path / "harmonization/harmonization_qc/study_a/report.html",
        "Study A QC | PBMC",
    )
    _write_report(
        tmp_path / "harmonization/integration_benchmark/report.html",
        "Integration benchmark | PBMC",
    )
    (tmp_path / "harmonization/integration_benchmark/executed.ipynb").touch()
    _write_report(tmp_path / "analysis/combined/trajectory/report.html", "Trajectory | PBMC")
    (tmp_path / "analysis/combined/trajectory/executed.ipynb").touch()
    _write_report(
        tmp_path / "analysis/per_type/cd14-monocyte/analysis/report.html",
        "CD14 monocyte | PBMC",
    )
    _write_report(
        tmp_path / "analysis/per_type/cd14-monocyte/differential_expression/report.html",
        "CD14 monocyte DE | PBMC",
    )
    _write_report(
        tmp_path / "analysis/per_type/naive-cd4-t-cell/differential_expression/report.html",
        "Naive CD4 T cell | PBMC",
    )
    _write_report(
        tmp_path / "analysis/per_type/cd14-monocyte/ageing_clock/report.html",
        "Age-clock comparison | PBMC",
    )
    _write_report(
        tmp_path / "analysis/per_type/cd14-monocyte/ageing_clock/elastic_net/report.html",
        "Elastic net age clock | PBMC",
    )
    _write_report(
        tmp_path / "analysis/per_type/cd14-monocyte/ageing_clock/equal_study_elastic_net/report.html",
        "Study-balanced age clock | PBMC",
    )
    (tmp_path / "analysis/per_type/cd14-monocyte/ageing_clock/age_clock.json").write_text(
        '{"cell_type": "CD14 monocyte", "methods": ['
        '{"method": "elastic_net", "method_name": "Elastic net"}, '
        '{"method": "equal_study_elastic_net", '
        '"method_name": "Study-balanced elastic net"}]}',
        encoding="utf-8",
    )
    (tmp_path / "analysis/per_type/age-only-type/ageing_clock").mkdir(parents=True)
    (tmp_path / "analysis/per_type/age-only-type/ageing_clock/age_clock.json").write_text(
        '{"cell_type": "Age-only type", "methods": []}', encoding="utf-8"
    )
    _write_report(
        tmp_path / "analysis/per_type/age-only-type/ageing_clock/report.html",
        "Age-only comparison | PBMC",
    )
    (tmp_path / "analysis/per_type/cell_type_manifest.json").write_text(
        '{"cell_types": [{"slug": "cd14-monocyte"}]}',
        encoding="utf-8",
    )

    index_path = render_report_index(tmp_path)
    html = index_path.read_text(encoding="utf-8")

    assert index_path == tmp_path / "index.html"
    assert 'href="harmonization/harmonization_qc/merged/report.html"' in html
    assert 'href="harmonization/harmonization_qc/study_a/report.html"' in html
    assert 'href="harmonization/integration_benchmark/report.html"' in html
    assert 'href="analysis/combined/trajectory/report.html"' in html
    assert 'href="analysis/per_type/cd14-monocyte/analysis/report.html"' in html
    assert 'href="analysis/per_type/cd14-monocyte/differential_expression/report.html"' in html
    assert 'href="analysis/per_type/naive-cd4-t-cell/differential_expression/report.html"' in html
    assert 'href="analysis/per_type/cd14-monocyte/ageing_clock/report.html"' in html
    assert 'href="analysis/per_type/age-only-type/ageing_clock/report.html"' in html
    assert 'href="analysis/per_type/cd14-monocyte/ageing_clock/elastic_net/report.html"' in html
    assert (
        'href="analysis/per_type/cd14-monocyte/ageing_clock/'
        'equal_study_elastic_net/report.html"' in html
    )
    assert "Per-cell-type analysis" in html
    assert "Generated " in html
    assert 'datetime="' in html
    assert "DE & abundance" in unescape(html)
    assert "Age clock" in html
    assert "Study-balanced elastic net" in unescape(html)
    assert "Age-only type" in html
    assert "<code>age-only-type</code>" in html
    assert (
        "<table><thead><tr><th>Cell type</th><th>Single-cell</th>"
        "<th>DE &amp; abundance</th>"
    ) in html
    assert "Integration benchmark report (HTML)" not in html
    assert "Open report" not in html
    assert ".ipynb" not in html
    assert "Study A QC" in html

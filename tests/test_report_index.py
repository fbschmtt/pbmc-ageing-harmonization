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
    assert ".ipynb" not in html
    assert "Study A QC" in html

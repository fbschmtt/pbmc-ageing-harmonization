import json
import sys

import pytest

from pbmc_pipeline.age_clock_report import main


@pytest.mark.parametrize("cell_type", [None, "cd14-monocyte"])
def test_root_render_refreshes_each_selected_cell_type_comparison(tmp_path, monkeypatch, cell_type):
    report_root = tmp_path / "per_type"
    for label, slug in [("CD14 monocyte", "cd14-monocyte"), ("CD8 T cell", "cd8-t-cell")]:
        for method in ("elastic_net", "rank_elastic_net"):
            directory = report_root / slug / "ageing_clock" / method
            directory.mkdir(parents=True)
            (directory / "age_clock.json").write_text(json.dumps({
                "cell_types": [{"cell_type": label}],
            }))
    calls = []

    def render(**kwargs):
        calls.append((kwargs["report_dir"], kwargs["template_path"].name))
        return kwargs["report_dir"] / "report.html"

    monkeypatch.setattr("pbmc_pipeline.age_clock_report.generate_age_clock_report", render)
    arguments = [
        "pbmc-age-clock-report", "--report-root", str(report_root),
        "--template", "reports/age_clock_report.py",
        "--comparison-template", "reports/age_clock_comparison_report.py",
    ]
    if cell_type:
        arguments += ["--cell-type", cell_type]
    monkeypatch.setattr(sys, "argv", arguments)
    main()
    expected_slugs = {"cd14-monocyte"} if cell_type else {"cd14-monocyte", "cd8-t-cell"}
    comparisons = {directory for directory, template in calls if template == "age_clock_comparison_report.py"}
    assert comparisons == {report_root / slug / "ageing_clock" for slug in expected_slugs}
    methods = {directory for directory, template in calls if template == "age_clock_report.py"}
    assert methods == {
        report_root / slug / "ageing_clock" / method
        for slug in expected_slugs for method in ("elastic_net", "rank_elastic_net")
    }

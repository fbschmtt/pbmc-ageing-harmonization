import hashlib
import json
from datetime import datetime, timezone

import pytest

from scripts.check_cell_type_reports import check_reports
from scripts.check_de_publication import check_publication


@pytest.fixture
def cell_type_output(tmp_path):
    root = tmp_path / "analysis/per_type"
    directory = root / "test-type/analysis"
    directory.mkdir(parents=True)
    for name in ("analysis.h5ad", "report.html", "executed.ipynb"):
        (directory / name).write_text("test output")
    (directory / "report.json").write_text(json.dumps({"cell_type": "Test type", "status": "complete"}))
    image = "site-assets/test-image.png"
    (directory / "site-assets").mkdir()
    (directory / image).write_bytes(b"test image bytes")
    (directory / "site-assets.json").write_text(json.dumps({
        "schemaVersion": 1, "cellType": "Test type", "image": image,
        "assets": {image: {"sha256": hashlib.sha256(b"test image bytes").hexdigest()}},
    }))
    manifest = {"status": "complete", "created_utc": datetime.now(timezone.utc).isoformat(),
                "failed_cell_types": [], "cell_types": [{"slug": "test-type", "cell_type": "Test type",
                "status": "complete", "report_html_sha256": hashlib.sha256(b"test output").hexdigest()}]}
    (root / "cell_type_manifest.json").write_text(json.dumps(manifest))
    return tmp_path, root / "cell_type_manifest.json", directory


def test_cell_type_check_rejects_ignored_failures_and_stale_success(cell_type_output):
    outdir, manifest_path, _directory = cell_type_output
    check_reports(outdir, cell_types="test-type")
    with pytest.raises(ValueError, match="selection differs"):
        check_reports(outdir, cell_types="missing-type")
    with pytest.raises(ValueError, match="stale"):
        check_reports(outdir, started_after=datetime.now(timezone.utc).timestamp() + 1)
    manifest = json.loads(manifest_path.read_text())
    manifest.update(status="partial", failed_cell_types=[{"slug": "another-type", "failed_stage": "report"}])
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="failed stages"):
        check_reports(outdir)


def test_cell_type_check_requires_site_assets_and_rejects_changed_reports(cell_type_output):
    outdir, _, directory = cell_type_output
    (directory / "report.html").write_text("stale HTML from another run")
    with pytest.raises(ValueError, match="differs from its manifest"):
        check_reports(outdir)
    (directory / "site-assets.json").unlink()
    with pytest.raises(ValueError, match="Missing cell-type deliverable"):
        check_reports(outdir)


def test_complete_cell_type_requires_the_named_umap_image(cell_type_output):
    outdir, _, directory = cell_type_output
    assets_path = directory / "site-assets.json"
    assets = json.loads(assets_path.read_text())
    assets["image"] = None
    assets_path.write_text(json.dumps(assets))
    with pytest.raises(ValueError, match="Missing named UMAP"):
        check_reports(outdir)


def test_de_publication_requires_all_three_stages_to_coexist(tmp_path):
    directory = tmp_path / "analysis/per_type/test-type/differential_expression"
    (directory / "residual_wgcna").mkdir(parents=True)
    for name in ("diagnostics.json", "fraction_model_coefficients.tsv", "residual_wgcna/metadata.json",
                 "report.html", "executed.ipynb", "report_render.log"):
        (directory / name).write_text("test output")
    (directory / "run_metadata.json").write_text(json.dumps({"cell_type": "Test type"}))
    (directory / "merged.csv").write_text("gene,log2FoldChange\nTEST_GENE,0.25\n")
    check_publication(tmp_path)
    (directory / "merged.csv").unlink()
    with pytest.raises(ValueError, match="lost required artifacts"):
        check_publication(tmp_path)

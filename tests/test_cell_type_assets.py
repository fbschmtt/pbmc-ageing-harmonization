import base64
import hashlib
import json

import pytest

from pbmc_pipeline.cell_type_assets import CellTypeAssets

PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII="
)


def test_explicit_assets_preserve_figures_fraction_semantics_and_checksums(tmp_path):
    assets = CellTypeAssets(tmp_path, "CD14 monocyte")
    image = assets.image("umap", PNG)
    assets.note("ageDeNote", "No age-DE results")
    fits = [{"study": None, "samples": 3, "meanFraction": 0.4, "slopePpPerDecade": 1.2}]
    assets.fractions([{"cluster": "0", "image": base64.b64encode(PNG).decode(), "fits": fits}],
                     note="One sample excluded")
    manifest = json.loads(assets.publish().read_text())
    assert manifest["image"] == image
    assert manifest["fractionDetails"] == [{"cluster": "0", "image": image, "fits": fits}]
    assert manifest["fractionIncludesZeros"] is True
    assert manifest["fractionHasStudyTrends"] is True
    assert manifest["assets"][image]["sha256"] == hashlib.sha256((tmp_path / image).read_bytes()).hexdigest()
    assert manifest["ageDeImage"] is None
    assert manifest["ageDeNote"] == "No age-DE results"


def test_asset_collection_rejects_invalid_image_bytes(tmp_path):
    with pytest.raises(ValueError, match="must be PNG"):
        CellTypeAssets(tmp_path, "T cell").image("umap", b"invalid")

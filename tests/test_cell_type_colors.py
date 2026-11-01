import anndata as ad
import pandas as pd

from pbmc_pipeline.cell_type_colors import aifi_l2_colors, set_aifi_l2_scanpy_palette


def test_optional_palette_is_loaded_when_present(tmp_path, monkeypatch):
    path = tmp_path / "colors.csv"
    path.write_text("AIFI_L2,AIFI_L2_color\nNaive CD4 T cell,#123456\n")
    monkeypatch.setenv("PBMC_AIFI_L2_COLORS", str(path))
    assert aifi_l2_colors() == {"Naive CD4 T cell": "#123456"}


def test_scanpy_palette_uses_aifi_colors_in_category_order(tmp_path, monkeypatch):
    path = tmp_path / "colors.csv"
    path.write_text(
        "AIFI_L2,AIFI_L2_color\nNaive CD4 T cell,#123456\nMemory CD4 T cell,#abcdef\n"
    )
    monkeypatch.setenv("PBMC_AIFI_L2_COLORS", str(path))
    adata = ad.AnnData(
        obs=pd.DataFrame(
            {"aifi_l2_majority": pd.Categorical(["Naive CD4 T cell", "Memory CD4 T cell"])}
        )
    )

    set_aifi_l2_scanpy_palette(adata)

    assert adata.uns["aifi_l2_majority_colors"] == [
        "#abcdef",
        "#123456",
    ]


def test_missing_palette_uses_stable_colors_across_category_subsets(tmp_path, monkeypatch):
    monkeypatch.setenv("PBMC_AIFI_L2_COLORS", str(tmp_path / "missing.csv"))
    assert aifi_l2_colors() == {}
    both = ad.AnnData(obs=pd.DataFrame({"aifi_l2_majority": ["Naive CD4 T cell", "Treg"]}))
    one = ad.AnnData(obs=pd.DataFrame({"aifi_l2_majority": ["Treg"]}))
    set_aifi_l2_scanpy_palette(both)
    set_aifi_l2_scanpy_palette(one)
    assert one.uns["aifi_l2_majority_colors"][0] == both.uns["aifi_l2_majority_colors"][1]
    assert both.uns["aifi_l2_majority_colors"][0] != both.uns["aifi_l2_majority_colors"][1]

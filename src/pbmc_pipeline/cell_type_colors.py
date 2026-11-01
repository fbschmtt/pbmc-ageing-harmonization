"""Shared AIFI cell-type colors for comparable report plots."""
import colorsys
import csv
import hashlib
import os
from pathlib import Path

_COLORSET_NAME = "AIFI_L2_imm_health_atlas_type_order_colors.csv"


def aifi_l2_colors() -> dict[str, str]:
    """Load the optional published AIFI L2 label-to-color mapping."""
    default = Path(__file__).resolve().parents[2] / "aifi_models" / _COLORSET_NAME
    path = Path(os.environ.get("PBMC_AIFI_L2_COLORS", str(default)))
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return {
                row["AIFI_L2"]: row["AIFI_L2_color"]
                for row in csv.DictReader(handle)
                if row.get("AIFI_L2") and row.get("AIFI_L2_color")
            }
    except FileNotFoundError:
        return {}


def _fallback_color(label: str) -> str:
    """Derive a stable color for a label, independent of category order."""
    hue = int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:4], "big") / 2**32
    rgb = colorsys.hsv_to_rgb(hue, 0.62, 0.78)
    return "#" + "".join(f"{round(component * 255):02x}" for component in rgb)


def set_aifi_l2_scanpy_palette(adata, column: str = "aifi_l2_majority") -> None:
    """Set Scanpy's categorical palette for an AIFI L2 observation column."""
    if column not in adata.obs:
        return
    values = adata.obs[column]
    categories = values.cat.categories if hasattr(values, "cat") else values.astype("category").cat.categories
    colors = aifi_l2_colors()
    adata.uns[f"{column}_colors"] = [
        colors.get(str(label)) or _fallback_color(str(label)) for label in categories
    ]

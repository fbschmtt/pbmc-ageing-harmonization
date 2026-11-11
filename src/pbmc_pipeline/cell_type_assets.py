"""Explicit, versioned figures and metadata consumed by the static site."""
from __future__ import annotations

import base64
import hashlib
import json
from io import BytesIO
from pathlib import Path


class CellTypeAssets:
    """Collect report figures without coupling export to notebook presentation."""

    def __init__(self, directory: Path, cell_type: str):
        self.directory = directory
        self.assets: dict[str, dict] = {}
        self.metadata = {
            "schemaVersion": 1, "cellType": cell_type, "image": None,
            "ageDeImage": None, "ageDeNote": None, "fractionNote": None,
            "fractionHasStudyTrends": False, "fractionIncludesZeros": False,
            "fractionDetails": None,
        }

    def image(self, name: str, png: bytes) -> str:
        if not png.startswith(b"\x89PNG\r\n\x1a\n"):
            raise ValueError(f"{name}: report image must be PNG")
        filename = f"site-assets/{hashlib.sha256(png).hexdigest()}.png"
        path = self.directory / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(png)
        self.assets[filename] = {"sha256": hashlib.sha256(png).hexdigest()}
        self.metadata["image" if name == "umap" else name] = filename
        return filename

    def figure(self, name: str, figure) -> str:
        buffer = BytesIO()
        figure.savefig(buffer, format="png", bbox_inches="tight")
        return self.image(name, buffer.getvalue())

    def note(self, name: str, text: str) -> None:
        self.metadata[name] = text

    def fractions(self, clusters: list[dict], *, note: str) -> None:
        details = []
        for cluster in clusters:
            filename = self.image(
                f"fraction-{cluster['cluster']}", base64.b64decode(cluster["image"], validate=True),
            )
            # Fraction image entries are referenced from their cluster records.
            self.metadata.pop(f"fraction-{cluster['cluster']}")
            details.append({**cluster, "image": filename})
        self.metadata.update({
            "fractionDetails": details, "fractionNote": note,
            "fractionHasStudyTrends": True, "fractionIncludesZeros": True,
        })

    def publish(self) -> Path:
        path = self.directory / "site-assets.json"
        path.write_text(json.dumps({**self.metadata, "assets": self.assets}, indent=2) + "\n")
        return path

"""Reject missing, stale or partially failed cell-type test workflow outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path


def check_reports(outdir: Path, *, cell_types: str = "", started_after: float = 0) -> None:
    root = outdir / "analysis/per_type"
    manifest = json.loads((root / "cell_type_manifest.json").read_text())
    created = datetime.fromisoformat(manifest["created_utc"]).timestamp()
    if created < started_after:
        raise ValueError("Cell-type manifest is stale; this workflow did not publish a new manifest")
    if manifest.get("status") != "complete" or manifest.get("failed_cell_types"):
        raise ValueError(f"Cell-type workflow has failed stages: {manifest.get('failed_cell_types')}")
    records = manifest.get("cell_types", [])
    found = {record["slug"] for record in records}
    selected = {value.strip() for value in cell_types.split(",") if value.strip()}
    if not records or len(found) != len(records) or (selected and selected != found):
        raise ValueError(f"Cell-type report selection differs: requested={sorted(selected)}, published={sorted(found)}")
    for record in records:
        directory = root / record["slug"] / "analysis"
        for name in ("analysis.h5ad", "report.json", "report.html", "executed.ipynb", "site-assets.json"):
            if not (directory / name).is_file():
                raise ValueError(f"Missing cell-type deliverable: {directory / name}")
        report = json.loads((directory / "report.json").read_text())
        assets = json.loads((directory / "site-assets.json").read_text())
        if (report.get("status") != record["status"] or report.get("cell_type") != record["cell_type"]
                or assets.get("schemaVersion") != 1 or assets.get("cellType") != record["cell_type"]):
            raise ValueError(f"Cell-type metadata disagrees: {record['slug']}")
        if report["status"] == "complete" and assets.get("image") not in assets["assets"]:
            raise ValueError(f"Missing named UMAP image: {record['slug']}")
        digest = hashlib.sha256((directory / "report.html").read_bytes()).hexdigest()
        if digest != record["report_html_sha256"]:
            raise ValueError(f"Cell-type report differs from its manifest: {record['slug']}")
        for name, metadata in assets["assets"].items():
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError(f"Unsafe cell-type image path: {name}")
            digest = hashlib.sha256((directory / name).read_bytes()).hexdigest()
            if digest != metadata["sha256"]:
                raise ValueError(f"Cell-type image checksum mismatch: {name}")
    print(f"Cell-type report check passed: {len(records)} selected types and their site assets.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cell-types", default="")
    parser.add_argument("--started-after", type=float, default=0)
    args = parser.parse_args()
    check_reports(args.outdir, cell_types=args.cell_types, started_after=args.started_after)


if __name__ == "__main__":
    main()

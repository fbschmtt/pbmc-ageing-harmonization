"""Write a compact, machine-readable manifest for published workflow artifacts."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Write provenance for a published pipeline run")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--artifact", type=Path, nargs="+", required=True)
    parser.add_argument("--report", type=Path, nargs="+", required=True)
    parser.add_argument("--requested-studies", required=True)
    parser.add_argument("--studies", required=True)
    parser.add_argument("--skipped-studies", default="")
    parser.add_argument("--pipeline-revision", required=True)
    parser.add_argument("--python-image", required=True)
    parser.add_argument("--r-image", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "requested_studies": [study for study in args.requested_studies.split(",") if study],
        "selected_studies": [study for study in args.studies.split(",") if study],
        "skipped_studies": [study for study in args.skipped_studies.split(",") if study],
        "pipeline_revision": args.pipeline_revision,
        "images": {"python": args.python_image, "r_conversion": args.r_image},
        "configuration": {"path": str(args.config), "sha256": _sha256(args.config)},
        "artifacts": [{"path": str(path), "sha256": _sha256(path)} for path in args.artifact],
        "reports": [{"path": str(path), "sha256": _sha256(path)} for path in args.report],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")

"""Best-effort, opt-in downloader for public pipeline inputs.

This convenience command intentionally does not authenticate, scrape portals, or
guarantee source availability. It records what happened and leaves manual-only
sources clearly visible to the caller.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    partial = destination.with_name(destination.name + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "pbmc-ageing-pipeline/0.1"})
    with urllib.request.urlopen(request) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    partial.replace(destination)


def extract(archive: Path, member: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if archive.suffix == ".zip":
        with zipfile.ZipFile(archive) as bundle, bundle.open(member) as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)
        return
    with tarfile.open(archive) as bundle:
        source = bundle.extractfile(member)
        if source is None:
            raise ValueError(f"Archive member is not a file: {member}")
        with source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)


def selected_artifacts(document: dict, studies: set[str]):
    return [item for item in document["artifacts"] if studies.intersection(item["studies"])]


def main() -> None:
    parser = argparse.ArgumentParser(description="Best-effort download of configured public inputs")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("config/input_sources.json"))
    parser.add_argument("--studies", default="all", help="Comma-separated study IDs or all")
    parser.add_argument("--force", action="store_true", help="Redownload files that already exist")
    parser.add_argument("--require-all", action="store_true", help="Fail if any selected artifact is manual-only")
    args = parser.parse_args()

    root = args.project_root.resolve()
    document = json.loads((root / args.manifest).read_text())
    all_studies = {study for item in document["artifacts"] for study in item["studies"]}
    studies = all_studies if args.studies == "all" else {value.strip() for value in args.studies.split(",") if value.strip()}
    unknown = studies - all_studies
    if unknown:
        parser.error(f"unknown studies: {', '.join(sorted(unknown))}")

    records = []
    for item in selected_artifacts(document, studies):
        destination = root / item["path"]
        record = {"id": item["id"], "path": item["path"], "timestamp_utc": datetime.now(timezone.utc).isoformat()}
        if destination.exists() and not args.force:
            observed_size = destination.stat().st_size
            record["observed_size_bytes"] = observed_size
            if item.get("expected_size_bytes") not in (None, observed_size):
                record["status"] = "size_mismatch"
                record["expected_size_bytes"] = item["expected_size_bytes"]
            else:
                record["status"] = "present"
        elif "url" not in item:
            record.update({"status": "manual", "manual_url": item.get("manual_url"), "note": item.get("note")})
        else:
            try:
                if item.get("archive_member"):
                    archive = destination.with_name(destination.name + ".download")
                    download(item["url"], archive)
                    extract(archive, item["archive_member"], destination)
                    archive.unlink()
                else:
                    download(item["url"], destination)
                record["status"] = "downloaded"
            except (OSError, ValueError, tarfile.TarError, urllib.error.URLError, zipfile.BadZipFile) as error:
                record.update({"status": "failed", "error": str(error), "source_page": item.get("source_page")})
        if destination.exists():
            record["observed_size_bytes"] = destination.stat().st_size
            record["sha256"] = sha256(destination)
            if item.get("sha256") and record["sha256"] != item["sha256"]:
                record["status"] = "checksum_mismatch"
            expected_size = item.get("expected_size_bytes")
            if (
                record["status"] != "failed"
                and expected_size is not None
                and record["observed_size_bytes"] != expected_size
            ):
                record["status"] = "size_mismatch"
                record["expected_size_bytes"] = expected_size
        records.append(record)
        print(f"{record['status']:>17}  {item['id']}  {destination}")
        if record["status"] == "manual":
            print(f"  obtain manually: {record.get('manual_url')}\n  {record.get('note', '')}")

    report_path = root / "input_data/download_manifest.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps({"records": records}, indent=2) + "\n")
    incomplete = {"manual", "failed", "checksum_mismatch", "size_mismatch"}
    if args.require_all and any(record["status"] in incomplete for record in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

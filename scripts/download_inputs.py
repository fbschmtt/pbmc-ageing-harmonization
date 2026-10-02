"""Best-effort, opt-in downloader for public pipeline inputs.

This convenience command intentionally does not authenticate, scrape portals, or
guarantee source availability. It records what happened and leaves manual-only
sources clearly visible to the caller.
"""
from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

_MAX_DOWNLOAD_ATTEMPTS = 3
_RETRY_DELAYS_SECONDS = (2, 5)


class _RetryableDownloadError(OSError):
    """A transient transfer failure for which a byte-range retry is useful."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _partial_metadata_path(partial: Path) -> Path:
    return partial.with_name(partial.name + ".json")


def _load_partial_metadata(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _discard_partial(partial: Path) -> None:
    partial.unlink(missing_ok=True)
    _partial_metadata_path(partial).unlink(missing_ok=True)


def _partial_matches_request(
    metadata: dict | None, *, url: str, expected_size: int | None,
    expected_sha256: str | None,
) -> bool:
    return metadata == {
        "url": url,
        "expected_size_bytes": expected_size,
        "expected_sha256": expected_sha256,
    } or (
        isinstance(metadata, dict)
        and metadata.get("url") == url
        and metadata.get("expected_size_bytes") == expected_size
        and metadata.get("expected_sha256") == expected_sha256
    )


def _write_partial_metadata(
    path: Path, *, url: str, expected_size: int | None, expected_sha256: str | None,
    response,
) -> None:
    path.write_text(json.dumps({
        "url": url,
        "expected_size_bytes": expected_size,
        "expected_sha256": expected_sha256,
        "etag": response.headers.get("ETag"),
        "last_modified": response.headers.get("Last-Modified"),
    }, indent=2) + "\n")


def _response_matches_partial(metadata: dict, response) -> bool:
    for key, header in (("etag", "ETag"), ("last_modified", "Last-Modified")):
        previous, current = metadata.get(key), response.headers.get(header)
        if previous is not None and current is not None and previous != current:
            return False
    return True


def _retryable_error(error: Exception) -> bool:
    if isinstance(error, urllib.error.HTTPError):
        return error.code in {408, 429} or error.code >= 500
    return isinstance(error, (http.client.HTTPException, OSError, urllib.error.URLError))


def _report_progress(
    *, copied: int, reported: int, resume_from: int, total: int | None, started: float,
    terminal: bool, step: int, force: bool,
) -> int:
    if not force and copied - reported < step:
        return reported
    rate = format_size(int((copied - resume_from) / max(time.monotonic() - started, 0.001))) + "/s"
    if total:
        percent = min(copied * 100 // total, 100)
        message = f"      {percent:3d}%  {format_size(copied)} / {format_size(total)}  {rate}"
    else:
        message = f"      {format_size(copied)} downloaded  {rate}"
    print(message, end="\r" if terminal else "\n", file=sys.stderr, flush=True)
    return copied


def download(
    url: str,
    destination: Path,
    expected_size: int | None = None,
    expected_sha256: str | None = None,
) -> dict:
    partial = destination.with_name(destination.name + ".part")
    partial_metadata_path = _partial_metadata_path(partial)
    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata = _load_partial_metadata(partial_metadata_path)
    if partial.exists() and not _partial_matches_request(
        metadata, url=url, expected_size=expected_size, expected_sha256=expected_sha256,
    ):
        print("      Discarding partial download with incompatible or missing metadata.", file=sys.stderr)
        _discard_partial(partial)
        metadata = None

    for attempt in range(1, _MAX_DOWNLOAD_ATTEMPTS + 1):
        resume_from = partial.stat().st_size if partial.exists() else 0
        if expected_size is not None and resume_from > expected_size:
            _discard_partial(partial)
            metadata = None
            resume_from = 0
        if expected_size is not None and resume_from == expected_size:
            if expected_sha256 is not None and sha256(partial) != expected_sha256:
                _discard_partial(partial)
                metadata = None
                continue
            partial.replace(destination)
            partial_metadata_path.unlink(missing_ok=True)
            return {"attempts": attempt - 1, "resumed_from_bytes": resume_from}
        headers = {"User-Agent": "pbmc-ageing-pipeline/0.1"}
        if resume_from:
            headers["Range"] = f"bytes={resume_from}-"
            print(
                f"      Resuming at {format_size(resume_from)} (attempt {attempt}/{_MAX_DOWNLOAD_ATTEMPTS}).",
                file=sys.stderr,
            )
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request) as response:
                status = response.getcode()
                content_range = response.headers.get("Content-Range", "")
                if resume_from and (
                    status != 206 or not content_range.startswith(f"bytes {resume_from}-")
                ):
                    print(
                        "      Server did not accept a compatible byte-range request; restarting download.",
                        file=sys.stderr,
                    )
                    _discard_partial(partial)
                    return download(url, destination, expected_size, expected_sha256)
                if resume_from and metadata is not None and not _response_matches_partial(metadata, response):
                    print(
                        "      Source validators changed since the partial download; restarting download.",
                        file=sys.stderr,
                    )
                    _discard_partial(partial)
                    return download(url, destination, expected_size, expected_sha256)
                content_length_header = response.headers.get("Content-Length")
                content_length = (
                    int(content_length_header)
                    if content_length_header and content_length_header.isdigit()
                    else None
                )
                total = resume_from + content_length if content_length is not None else expected_size
                if expected_size is not None and total is not None and total != expected_size:
                    raise ValueError(
                        f"server advertised {total} bytes; expected size is {expected_size}"
                    )
                _write_partial_metadata(
                    partial_metadata_path,
                    url=url,
                    expected_size=expected_size,
                    expected_sha256=expected_sha256,
                    response=response,
                )
                copied = resume_from
                reported = copied
                started = time.monotonic()
                terminal = sys.stderr.isatty()
                step = max((total or 0) // 20, 10 * 1024 * 1024) if total else 100 * 1024 * 1024

                with partial.open("ab" if resume_from else "wb") as handle:
                    while block := response.read(1024 * 1024):
                        handle.write(block)
                        copied += len(block)
                        reported = _report_progress(
                            copied=copied, reported=reported, resume_from=resume_from,
                            total=total, started=started, terminal=terminal, step=step, force=False,
                        )
                _report_progress(
                    copied=copied, reported=reported, resume_from=resume_from,
                    total=total, started=started, terminal=terminal, step=step, force=True,
                )
                if terminal:
                    print(file=sys.stderr, flush=True)
            if content_length is not None and copied != total:
                raise _RetryableDownloadError(
                    f"download ended at {copied} bytes; server advertised {total}"
                )
            if expected_size is not None and copied != expected_size:
                raise _RetryableDownloadError(
                    f"download ended at {copied} bytes; expected {expected_size}"
                )
            if expected_sha256 is not None:
                observed_sha256 = sha256(partial)
                if observed_sha256 != expected_sha256:
                    _discard_partial(partial)
                    raise ValueError(
                        f"download checksum {observed_sha256} differs from expected {expected_sha256}"
                    )
            partial.replace(destination)
            partial_metadata_path.unlink(missing_ok=True)
            return {"attempts": attempt, "resumed_from_bytes": resume_from}
        except Exception as error:  # Keep only transient transport failures retryable.
            if not _retryable_error(error):
                raise
            preserved = partial.stat().st_size if partial.exists() else 0
            metadata = _load_partial_metadata(partial_metadata_path)
            if attempt == _MAX_DOWNLOAD_ATTEMPTS:
                print(
                    f"      Download failed after {attempt} attempts; preserving {format_size(preserved)} "
                    "for a future resume.",
                    file=sys.stderr,
                )
                raise
            delay = _RETRY_DELAYS_SECONDS[attempt - 1]
            print(
                f"      Attempt {attempt}/{_MAX_DOWNLOAD_ATTEMPTS} failed after "
                f"{format_size(preserved)}: {error}. Retrying in {delay}s; the next attempt "
                "will resume if the server supports HTTP Range.",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise AssertionError("unreachable")


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


def format_size(size_bytes: int) -> str:
    if size_bytes < 1024:
        return f"{size_bytes} bytes"
    for unit in ("KiB", "MiB", "GiB", "TiB"):
        size_bytes /= 1024
        if size_bytes < 1024 or unit == "TiB":
            return f"{size_bytes:.1f} {unit}"
    raise AssertionError("unreachable")


def warn_about_pending_downloads(artifacts: list[dict], root: Path, force: bool) -> None:
    """Make the deliberately separate, potentially large acquisition visible."""
    pending = [
        item for item in artifacts
        if item.get("url") and (force or not (root / item["path"]).exists())
    ]
    if not pending:
        return
    known_size = sum(item.get("expected_size_bytes") or 0 for item in pending)
    unknown_sizes = sum(item.get("expected_size_bytes") is None for item in pending)
    size_summary = (
        f"known sizes total at least {format_size(known_size)}"
        if known_size
        else "no expected sizes recorded"
    )
    if unknown_sizes:
        unknown_ids = [item["id"] for item in pending if item.get("expected_size_bytes") is None]
        size_summary += f"; size unavailable for: {', '.join(unknown_ids)}"
    print(
        f"Input downloads: {len(pending)} pending; {size_summary}. "
        "Downloads run only through explicit Make targets.",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Best-effort download of configured public inputs")
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, default=Path("config/input_sources.json"))
    parser.add_argument("--studies", default="all", help="Comma-separated study IDs or all")
    parser.add_argument("--force", action="store_true", help="Redownload files that already exist")
    parser.add_argument(
        "--require-all",
        action="store_true",
        help="Fail unless every selected artifact is present and passes its recorded checks",
    )
    args = parser.parse_args()

    root = args.project_root.resolve()
    document = json.loads((root / args.manifest).read_text())
    all_studies = {study for item in document["artifacts"] for study in item["studies"]}
    studies = all_studies if args.studies == "all" else {value.strip() for value in args.studies.split(",") if value.strip()}
    unknown = studies - all_studies
    if unknown:
        parser.error(f"unknown studies: {', '.join(sorted(unknown))}")

    artifacts = selected_artifacts(document, studies)
    warn_about_pending_downloads(artifacts, root, args.force)
    records = []
    study_order = list(dict.fromkeys(
        study for item in artifacts for study in item["studies"] if study in studies
    ))
    for overall_index, item in enumerate(artifacts, start=1):
        associated_studies = [study for study in study_order if study in item["studies"]]
        print(
            f"\n[{overall_index}/{len(artifacts)}] {item['id']} "
            f"(study: {', '.join(associated_studies)})",
            flush=True,
        )
        destination = root / item["path"]
        record = {"id": item["id"], "path": item["path"], "timestamp_utc": datetime.now(timezone.utc).isoformat()}
        if destination.exists() and (not args.force or "url" not in item):
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
                    transfer = download(
                        item["url"], archive, item.get("expected_size_bytes"), item.get("sha256")
                    )
                    extract(archive, item["archive_member"], destination)
                    archive.unlink()
                else:
                    transfer = download(
                        item["url"], destination, item.get("expected_size_bytes"),
                        item.get("sha256"),
                    )
                record.update({"status": "downloaded", "transfer": transfer})
            except (
                OSError, ValueError, http.client.HTTPException, tarfile.TarError,
                urllib.error.URLError, zipfile.BadZipFile,
            ) as error:
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
        observed_size = record.get("observed_size_bytes")
        expected_size = item.get("expected_size_bytes")
        size_parts = []
        if observed_size is not None:
            size_parts.append(f"actual {format_size(observed_size)}")
        if expected_size is not None:
            size_parts.append(f"expected {format_size(expected_size)}")
        size_text = ", ".join(size_parts) if size_parts else "size unavailable"
        print(
            f"{record['status']:>17}  {item['id']}  {size_text}  {destination}",
            flush=True,
        )
        if record["status"] == "manual":
            print(f"  source: {record.get('manual_url')}\n  reason: {record.get('note', '')}")

    report_path = root / "input_data/download_manifest.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps({"records": records}, indent=2) + "\n")
    invalid_downloads = {"failed", "checksum_mismatch", "size_mismatch"}
    if any(record["status"] in invalid_downloads for record in records):
        raise SystemExit(1)
    if args.require_all and any(record["status"] == "manual" for record in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

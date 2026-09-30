"""Explicit local contract check for the opt-in input downloader."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOWNLOADER = ROOT / "scripts" / "download_inputs.py"


def run_downloader(project_root: Path, manifest: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable, str(DOWNLOADER), "--project-root", str(project_root),
            "--manifest", str(manifest), *arguments,
        ],
        check=False, capture_output=True, text=True,
    )


def write_manifest(path: Path, artifacts: list[dict]) -> None:
    path.write_text(json.dumps({"schema_version": 1, "artifacts": artifacts}))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="pbmc-download-check-") as temporary_dir:
        root = Path(temporary_dir)
        project_root = root / "project"
        project_root.mkdir()
        source = root / "source.bin"
        source.write_bytes(b"verified fixture")
        manifest = root / "sources.json"
        write_manifest(manifest, [{
            "id": "fixture", "studies": ["demo"], "path": "input_data/demo/fixture.bin",
            "url": source.as_uri(),
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "expected_size_bytes": source.stat().st_size,
        }])
        direct = run_downloader(project_root, manifest, "--studies", "demo", "--require-all")
        destination = project_root / "input_data/demo/fixture.bin"
        if (
            direct.returncode
            or destination.read_bytes() != source.read_bytes()
            or "input acquisition is strictly opt-in" not in direct.stderr
        ):
            raise SystemExit(f"direct-download check failed:\n{direct.stdout}{direct.stderr}")

        write_manifest(manifest, [{
            "id": "manual_fixture", "studies": ["demo"], "path": "input_data/demo/manual.bin",
            "manual_url": "https://example.test/manual", "note": "Download after accepting terms.",
        }])
        manual = run_downloader(project_root, manifest, "--studies", "demo", "--require-all")
        if manual.returncode != 1:
            raise SystemExit(f"strict manual-input check failed:\n{manual.stdout}{manual.stderr}")
    print("Downloader local contract checks passed.")


if __name__ == "__main__":
    main()

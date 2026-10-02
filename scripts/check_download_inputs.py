"""Explicit local contract check for the opt-in input downloader."""
from __future__ import annotations

import hashlib
import http.server
import json
import subprocess
import sys
import tempfile
import threading
from pathlib import Path
from typing import ClassVar

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


def start_range_server(payload: bytes):
    """Serve one deliberately interrupted range response, then resume normally."""
    class Handler(http.server.BaseHTTPRequestHandler):
        requests: ClassVar[list[str | None]] = []
        interrupted: ClassVar[bool] = False

        def do_GET(self):
            requested_range = self.headers.get("Range")
            type(self).requests.append(requested_range)
            start = int(requested_range.removeprefix("bytes=").removesuffix("-") or 0) \
                if requested_range else 0
            body = payload[start:]
            self.send_response(206 if requested_range else 200)
            if requested_range:
                self.send_header("Content-Range", f"bytes {start}-{len(payload) - 1}/{len(payload)}")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("ETag", '"fixture-v1"')
            self.end_headers()
            if requested_range and not type(self).interrupted:
                type(self).interrupted = True
                self.wfile.write(body[:3])
                self.wfile.flush()
                self.connection.shutdown(2)
                return
            self.wfile.write(body)

        def log_message(self, _format, *_args):
            pass

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, Handler


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
        ):
            raise SystemExit(f"direct-download check failed:\n{direct.stdout}{direct.stderr}")

        payload = b"range-resume-fixture"
        server, handler = start_range_server(payload)
        try:
            resume_path = project_root / "input_data/demo/resume.bin"
            partial = resume_path.with_name(resume_path.name + ".part")
            partial.parent.mkdir(parents=True, exist_ok=True)
            partial.write_bytes(payload[:4])
            partial.with_name(partial.name + ".json").write_text(json.dumps({
                "url": f"http://127.0.0.1:{server.server_port}/fixture.bin",
                "expected_size_bytes": len(payload),
                "expected_sha256": hashlib.sha256(payload).hexdigest(),
                "etag": '"fixture-v1"',
                "last_modified": None,
            }))
            write_manifest(manifest, [{
                "id": "resume_fixture", "studies": ["demo"],
                "path": "input_data/demo/resume.bin",
                "url": f"http://127.0.0.1:{server.server_port}/fixture.bin",
                "sha256": hashlib.sha256(payload).hexdigest(),
                "expected_size_bytes": len(payload),
            }])
            resumed = run_downloader(project_root, manifest, "--studies", "demo", "--require-all")
            if (
                resumed.returncode
                or resume_path.read_bytes() != payload
                or handler.requests != ["bytes=4-", "bytes=7-"]
                or "Retrying in 2s" not in resumed.stderr
                or "Resuming at 7 bytes" not in resumed.stderr
            ):
                raise SystemExit(
                    f"resume-and-retry check failed:\n{resumed.stdout}{resumed.stderr}\n"
                    f"range requests: {handler.requests}"
                )
        finally:
            server.shutdown()
            server.server_close()

        write_manifest(manifest, [{
            "id": "manual_fixture", "studies": ["demo"], "path": "input_data/demo/manual.bin",
            "manual_url": "https://example.test/manual", "note": "Download after accepting terms.",
        }])
        manual = run_downloader(project_root, manifest, "--studies", "demo", "--require-all")
        if manual.returncode != 1:
            raise SystemExit(f"strict manual-input check failed:\n{manual.stdout}{manual.stderr}")
    print("Downloader local contract checks passed, including interrupted range resume.")


if __name__ == "__main__":
    main()

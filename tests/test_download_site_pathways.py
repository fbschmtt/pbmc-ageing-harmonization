import hashlib
import json

import pytest

from scripts import download_site_pathways as downloader


def test_extracts_examples_and_preserves_existing_files_and_labels(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "COLLECTIONS", {
        "example": ("digest", {"SET_A": "A", "SET_B": "B"}),
    })
    existing = tmp_path / f"SET_A.v{downloader.VERSION}.gmt"
    existing.write_text("custom\n")
    (tmp_path / "pathway-labels.json").write_text('{"SET_A": "Custom"}')
    monkeypatch.setattr(
        downloader, "download_collection", lambda *_: "SET_A\turl\tOLD\nSET_B\turl\tGENE1\tGENE2\n"
    )

    downloader.acquire(tmp_path)

    assert existing.read_text() == "custom\n"
    assert (tmp_path / f"SET_B.v{downloader.VERSION}.gmt").read_text() == "SET_B\turl\tGENE1\tGENE2\n"
    assert json.loads((tmp_path / "pathway-labels.json").read_text()) == {"SET_A": "Custom", "SET_B": "B"}


def test_missing_set_does_not_publish_partial_downloads(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "COLLECTIONS", {
        "first": ("digest", {"SET_A": "A"}),
        "second": ("digest", {"SET_B": "B"}),
    })
    monkeypatch.setattr(downloader, "download_collection", lambda *_: "SET_A\turl\tGENE\n")

    with pytest.raises(ValueError, match="Missing gene sets"):
        downloader.acquire(tmp_path)
    assert not list(tmp_path.iterdir())


def test_collection_checksum_is_verified(monkeypatch):
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_):
            pass

        def read(self):
            return b"SET\turl\tGENE\n"

    monkeypatch.setattr(downloader.urllib.request, "urlopen", lambda *_args, **_kwargs: Response())
    with pytest.raises(ValueError, match="Checksum mismatch"):
        downloader.download_collection("test", "incorrect")
    expected = hashlib.sha256(b"SET\turl\tGENE\n").hexdigest()
    assert downloader.download_collection("test", expected) == "SET\turl\tGENE\n"

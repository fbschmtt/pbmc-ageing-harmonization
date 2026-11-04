"""Download the pinned MSigDB collections and extract the atlas example sets."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

VERSION = "2026.1.Hs"
BASE_URL = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{VERSION}"
COLLECTIONS = {
    "h.all": ("eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596", {
        "HALLMARK_DNA_REPAIR": "DNA repair",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    }),
    "c5.go.bp": ("9be09dd06d6652566eb52eed530d62e6dfecc4365c1e81afd6f0b7f2e86dd4f9", {
        "GOBP_REGULATION_OF_TELOMERE_MAINTENANCE": "Telomere maintenance",
        "GOBP_RESPONSE_TO_NUTRIENT": "Nutrient response",
    }),
    "c2.cp.reactome": ("5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341", {
        "REACTOME_CELLULAR_SENESCENCE": "Cellular senescence",
        "REACTOME_MITOPHAGY": "Mitophagy",
    }),
    "c2.cp.wikipathways": ("42c7a6dedea6a7abbc37f5d77d1ab84cc8fdaad44d7f664026047d23cc78f52c", {
        "WP_AUTOPHAGY": "Autophagy",
        "WP_DNA_DAMAGE_RESPONSE": "DNA damage response",
        "WP_OXIDATIVE_STRESS_RESPONSE": "Oxidative stress response",
    }),
}


def download_collection(collection: str, expected_sha256: str) -> str:
    url = f"{BASE_URL}/{collection}.v{VERSION}.symbols.gmt"
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read()
    observed = hashlib.sha256(data).hexdigest()
    if observed != expected_sha256:
        raise ValueError(f"Checksum mismatch for {collection}: {observed}")
    return data.decode("utf-8")


def acquire(output: Path, *, force: bool = False) -> None:
    labels_path = output / "pathway-labels.json"
    labels = json.loads(labels_path.read_text(encoding="utf-8")) if labels_path.exists() else {}
    if not isinstance(labels, dict) or not all(isinstance(v, str) for v in labels.values()):
        raise ValueError(f"Invalid pathway labels: {labels_path}")
    original_labels = labels.copy()
    pending: dict[Path, str] = {}
    for collection, (digest, examples) in COLLECTIONS.items():
        needed = {
            name: output / f"{name}.v{VERSION}.gmt"
            for name in examples
            if force or not (output / f"{name}.v{VERSION}.gmt").exists()
        }
        for name, label in examples.items():
            labels.setdefault(name, label)
        if not needed:
            continue
        found: dict[str, str] = {}
        for line in download_collection(collection, digest).splitlines():
            name = line.split("\t", 1)[0]
            if name in needed:
                if name in found or len(line.split("\t")) < 3:
                    raise ValueError(f"Invalid or duplicate gene set {name} in {collection}")
                found[name] = line + "\n"
        missing = needed.keys() - found.keys()
        if missing:
            raise ValueError(f"Missing gene sets in {collection}: {', '.join(sorted(missing))}")
        pending.update({path: found[name] for name, path in needed.items()})
    output.mkdir(parents=True, exist_ok=True)
    for path, content in pending.items():
        path.write_text(content, encoding="utf-8")
        print(f"Wrote {path}")
    if labels != original_labels or not labels_path.exists():
        labels_path.write_text(json.dumps(labels, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Ready: {len(labels)} pathway labels; {len(pending)} gene sets downloaded")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("site/public/pathways"))
    parser.add_argument("--force", action="store_true", help="Replace existing example GMT files")
    args = parser.parse_args()
    try:
        acquire(args.output_dir, force=args.force)
    except (OSError, UnicodeError, ValueError, urllib.error.URLError) as error:
        print(f"Pathway download failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

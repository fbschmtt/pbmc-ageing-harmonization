"""Print the content-derived Docker tag for one repository build context."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = {
    "python": [".dockerignore", "docker/python.Dockerfile", "config", "src",
               "pyproject.toml", "requirements.lock", "aifi_models/README"],
    "r": ["docker/r-conversion.Dockerfile", "scripts/convert_rds.R"],
    "site": [".dockerignore", "docker/site.Dockerfile", "config/studies.json",
             "scripts/render_report_index.py", "site"],
}
SITE_EXCLUDED = (
    "site/node_modules/", "site/.npm-cache/", "site/.playwright-browsers/",
    "site/dist/", "site/public/data/", "site/data-inputs/", "site/public/pathways/",
)
SITE_INCLUDED = {"site/data-inputs/cell-type-hierarchy.json", "site/public/pathways/.gitkeep"}


def context_files(kind: str) -> list[str]:
    def tracked_and_untracked(*args):
        return subprocess.check_output(
            ["git", "ls-files", "-z", *args], cwd=ROOT,
        ).decode().split("\0")[:-1]

    paths = tracked_and_untracked("-co", "--exclude-standard", *CONTEXTS[kind])
    if kind == "python":
        color_table = "aifi_models/AIFI_L2_imm_health_atlas_type_order_colors.csv"
        if (ROOT / color_table).is_file():
            paths.append(color_table)
    if kind == "site":
        paths += tracked_and_untracked("-oi", "--exclude-standard", "site")
        paths = [path for path in paths if path in SITE_INCLUDED or (
            not path.startswith(SITE_EXCLUDED)
            and not {"__pycache__", ".ipynb_checkpoints"}.intersection(Path(path).parts)
            and not path.endswith(".log")
        )]
    return sorted(set(paths))


def context_hash(kind: str) -> str:
    digest = hashlib.sha256()
    for name in context_files(kind):
        path = ROOT / name
        if path.is_file():
            file_hash = hashlib.sha256()
            with path.open("rb") as source:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    file_hash.update(block)
            digest.update(f"{file_hash.hexdigest()}  {name}\n".encode())
    return digest.hexdigest()[:12]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("context", choices=CONTEXTS)
    print(context_hash(parser.parse_args().context))

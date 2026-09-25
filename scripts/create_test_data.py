#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pbmc_pipeline.config import load_configuration


def parse_args():
    parser = argparse.ArgumentParser(description="Create deterministic cell-downsampled test H5ADs")
    parser.add_argument("--cells", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--studies", default="all", help="Comma-separated study IDs or all")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    import scanpy as sc

    from pbmc_pipeline.source import read_study_input

    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    pipeline, document, _ = load_configuration(root, Path("config/pipeline.json"))
    requested = set(document["studies"]) if args.studies == "all" else {
        value.strip() for value in args.studies.split(",") if value.strip()
    }
    unknown = requested - set(document["studies"])
    if unknown:
        raise ValueError(f"Unknown studies: {', '.join(sorted(unknown))}")
    destination = root / pipeline["test_input_root"]
    destination.mkdir(parents=True, exist_ok=True)
    for study_id, study in document["studies"].items():
        if study_id not in requested:
            continue
        if study.get("conversion"):
            print(f"SKIP {study_id}: provide its downsampled RDS test fixture separately")
            continue
        source = root / study["input"]
        target = destination / study["test_input"]
        if not source.exists():
            print(f"SKIP {study_id}: missing {source}")
            continue
        if target.exists() and not args.overwrite:
            print(f"SKIP {study_id}: {target} exists")
            continue
        adata = read_study_input(root, source, study)
        if adata.n_obs > args.cells:
            sc.pp.sample(adata, n=args.cells, rng=args.seed, replace=False)
        adata.write_h5ad(target, compression="gzip")
        print(f"WROTE {study_id}: {target} ({adata.n_obs} cells)")


if __name__ == "__main__":
    main()

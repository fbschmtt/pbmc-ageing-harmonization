#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from pbmc_pipeline.config import load_configuration


def parse_args():
    parser = argparse.ArgumentParser(description="Create deterministic cell-downsampled test H5ADs")
    parser.add_argument("--cells", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    import scanpy as sc

    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    pipeline, document, _ = load_configuration(root, Path("config/pipeline.json"))
    destination = root / pipeline["test_input_root"]
    destination.mkdir(parents=True, exist_ok=True)
    for study_id, study in document["studies"].items():
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
        adata = sc.read_h5ad(source)
        if adata.n_obs > args.cells:
            sc.pp.sample(adata, n=args.cells, rng=args.seed, replace=False)
        adata.write_h5ad(target, compression="gzip")
        print(f"WROTE {study_id}: {target} ({adata.n_obs} cells)")


if __name__ == "__main__":
    main()

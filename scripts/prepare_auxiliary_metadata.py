#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create compact, pipeline-specific lookups from supplementary metadata"
    )
    parser.add_argument(
        "--terekhova-source",
        type=Path,
        default=Path("input_data/terekhova23/all_pbmcs_metadata.csv"),
    )
    parser.add_argument(
        "--terekhova-output",
        type=Path,
        default=Path("input_data/terekhova23/cell_to_tube.csv.gz"),
    )
    args = parser.parse_args()
    mapping = pd.read_csv(args.terekhova_source, usecols=["Unnamed: 0", "Tube_id"])
    mapping = mapping.rename(columns={"Unnamed: 0": "cell_id"})
    if mapping["cell_id"].duplicated().any():
        raise ValueError("Terekhova metadata contains duplicate cell IDs")
    if mapping["Tube_id"].isna().any():
        raise ValueError("Terekhova metadata contains cells without Tube_id")
    args.terekhova_output.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(args.terekhova_output, index=False, compression="gzip")
    print(f"Wrote {len(mapping):,} cell-to-tube mappings to {args.terekhova_output}")


if __name__ == "__main__":
    main()

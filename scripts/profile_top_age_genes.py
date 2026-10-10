"""Render focused gene profiles for top merged age-DE genes in each cell type."""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd


def _top_genes_by_cell_type(results_dir: Path, top_n: int) -> list[tuple[str, str]]:
    selected: list[tuple[str, str]] = []
    result_paths = sorted(results_dir.glob("*/merged.csv"))
    if not result_paths:
        raise FileNotFoundError(f"No merged DE result files found under {results_dir}")

    for result_path in result_paths:
        results = pd.read_csv(result_path)
        required = {"cell_type", "gene", "covariate", "padj"}
        missing = required - set(results.columns)
        if missing:
            raise ValueError(f"{result_path} is missing required columns: {sorted(missing)}")
        age = results.loc[results["covariate"].astype(str).str.lower() == "age"].copy()
        if age.empty:
            continue
        age = age.dropna(subset=["cell_type", "gene"])
        age["padj"] = pd.to_numeric(age["padj"], errors="coerce")
        age = age.loc[age["padj"].map(lambda value: math.isfinite(value) and value >= 0)]
        age = age.sort_values(["padj", "gene"], kind="stable").drop_duplicates("gene")
        cell_types = age["cell_type"].dropna().astype(str).unique()
        if len(cell_types) > 1:
            raise ValueError(f"{result_path} contains multiple cell types: {sorted(cell_types)}")
        if len(cell_types) == 0:
            continue
        cell_type = cell_types[0]
        selected.extend((cell_type, str(gene)) for gene in age["gene"].head(top_n))

    if not selected:
        raise ValueError(f"No finite merged age-DE adjusted p-values found under {results_dir}")
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Render detailed gene profiles for the top merged age-DE genes "
            "within each cell type."
        ),
    )
    parser.add_argument("--input", type=Path, required=True, help="Merged pseudobulk H5AD")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--differential-expression-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n must be at least 1")
    if not args.input.is_file():
        parser.error(f"Merged pseudobulk input does not exist: {args.input}")
    if not args.differential_expression_dir.is_dir():
        parser.error(
            "Differential-expression directory does not exist: "
            f"{args.differential_expression_dir}"
        )

    scripts_dir = Path(__file__).resolve().parent
    profile_script = scripts_dir / "profile_gene.py"
    selected = _top_genes_by_cell_type(args.differential_expression_dir, args.top_n)
    print(
        f"Rendering {len(selected)} gene profiles across "
        f"{len({cell_type for cell_type, _ in selected})} cell types."
    )
    for cell_type, gene in selected:
        print(f"\n=== {cell_type}: {gene} ===", flush=True)
        subprocess.run(
            [
                sys.executable,
                str(profile_script),
                "--input", str(args.input),
                "--config", str(args.config),
                "--cell-type", cell_type,
                "--gene", gene,
                "--output-dir", str(args.output_dir),
                "--cpus", str(args.cpus),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()

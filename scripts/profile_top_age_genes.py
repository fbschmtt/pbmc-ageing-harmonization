"""Render focused gene profiles for top merged age-DE genes in each cell type."""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


def _top_genes_by_cell_type(results_dir: Path, top_n: int, *, result_paths: list[Path] | None = None) -> list[tuple[str, str]]:
    selected: list[tuple[str, str]] = []
    if result_paths is None:
        result_paths = sorted(results_dir.glob("*/differential_expression/merged.csv"))
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

"""Reusable, offline preranked gene-set enrichment utilities."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


@dataclass(frozen=True)
class GeneSetEnrichmentConfig:
    """Source and statistical settings for a preranked GSEA run."""

    library_name: str
    organism: str
    gene_set_resource: str
    ontology_resource: str | None
    gene_set_source_url: str | None = None
    ontology_source_url: str | None = None
    minimum_gene_set_size: int = 10
    maximum_gene_set_size: int = 500
    minimum_ranked_genes: int = 10
    permutations: int = 1000
    seed: int = 123
    native_fdr_threshold: float = 0.05
    leading_edge_jaccard_threshold: float = 0.5


DEFAULT_GO_BP_CONFIG = GeneSetEnrichmentConfig(
    library_name="GO_Biological_Process_2023",
    organism="Human",
    gene_set_resource="go_biological_process_2023.gmt",
    ontology_resource="go-basic.obo.gz",
    gene_set_source_url=(
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&"
        "libraryName=GO_Biological_Process_2023"
    ),
    ontology_source_url="https://purl.obolibrary.org/obo/go/go-basic.obo",
)
GO_ENRICHMENT_DATA_DIR_ENV = "PBMC_GO_ENRICHMENT_DATA_DIR"


@dataclass
class GeneSetEnrichmentRun:
    """Results, provenance, and isolated failures from one enrichment batch."""

    results: pd.DataFrame
    gene_sets: dict[str, tuple[str, ...]]
    provenance: dict[str, object]
    failures: dict[str, str]


def _resource_path(name: str) -> Path:
    """Find a downloaded GO input, staged by Nextflow for report rendering."""
    root = Path(os.environ.get(GO_ENRICHMENT_DATA_DIR_ENV, "input_data/go"))
    path = root / name
    if not path.is_file():
        raise FileNotFoundError(
            f"GO enrichment input is missing: {path}. "
            "Run make download-inputs before rendering DE reports."
        )
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@cache
def load_gene_set_library(resource_name: str) -> dict[str, tuple[str, ...]]:
    """Load a downloaded GMT library without contacting Enrichr at report time."""
    path = _resource_path(resource_name)
    gene_sets = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            genes = tuple(gene.strip() for gene in fields[2:] if gene.strip())
            if genes:
                gene_sets[fields[0]] = genes
    if not gene_sets:
        raise ValueError(f"Gene-set library {path} is empty")
    return gene_sets


@cache
def load_go_ancestors(resource_name: str) -> tuple[dict[str, frozenset[str]], dict[str, str]]:
    """Load a downloaded GO DAG and its provenance for redundancy reduction."""
    path = _resource_path(resource_name)
    with gzip.open(path, mode="rt", encoding="utf-8") as handle:
        text = handle.read()
    parents_by_term: dict[str, tuple[str, ...]] = {}
    current_id: str | None = None
    parents: list[str] = []
    obsolete = False

    def save_term() -> None:
        if current_id and not obsolete:
            parents_by_term[current_id] = tuple(parents)

    for line in text.splitlines():
        if line.startswith("["):
            save_term()
            current_id, parents, obsolete = None, [], False
        elif line.startswith("id: GO:"):
            current_id = line.removeprefix("id: ").strip()
        elif line.startswith("is_a: GO:"):
            parents.append(line.split()[1])
        elif line.startswith("relationship:"):
            fields = line.split()
            if len(fields) >= 3 and fields[1] == "part_of" and fields[2].startswith("GO:"):
                parents.append(fields[2])
        elif line == "is_obsolete: true":
            obsolete = True
    save_term()

    resolved: dict[str, frozenset[str]] = {}

    def visit(term: str, active: frozenset[str]) -> frozenset[str]:
        if term in resolved:
            return resolved[term]
        if term in active:
            return frozenset()
        found = set()
        for parent in parents_by_term.get(term, ()):
            found.add(parent)
            found.update(visit(parent, active | {term}))
        resolved[term] = frozenset(found)
        return resolved[term]

    for term in parents_by_term:
        visit(term, frozenset())
    version = next(
        (line.removeprefix("data-version: ") for line in text.splitlines()
         if line.startswith("data-version: ")),
        "unknown",
    )
    return resolved, {
        "ontology_resource": resource_name,
        "ontology_path": str(path),
        "ontology_sha256": _sha256(path),
        "ontology_data_version": version,
    }


def _gene_set_signature(gene_sets: Mapping[str, tuple[str, ...]]) -> str:
    content = "\n".join(
        "\t".join((str(term), *sorted(map(str, genes))))
        for term, genes in sorted(gene_sets.items())
    )
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def enrichment_provenance(
    config: GeneSetEnrichmentConfig, gene_sets: Mapping[str, tuple[str, ...]] | None = None,
) -> dict[str, object]:
    """Return the exact local resources and settings used by an analysis."""
    downloaded_gene_sets = gene_sets is None
    resolved_gene_sets = gene_sets or load_gene_set_library(config.gene_set_resource)
    gene_set_path = _resource_path(config.gene_set_resource) if downloaded_gene_sets else None
    provenance: dict[str, object] = {
        "gene_set_library": config.library_name,
        "organism": config.organism,
        "gene_set_resource": config.gene_set_resource if downloaded_gene_sets else "caller-provided",
        "gene_set_path": str(gene_set_path) if gene_set_path else None,
        "gene_set_source_url": config.gene_set_source_url if downloaded_gene_sets else None,
        "gene_set_library_sha256": (
            _sha256(gene_set_path) if gene_set_path else _gene_set_signature(resolved_gene_sets)
        ),
        "gsea_parameters": asdict(config),
    }
    if config.ontology_resource:
        _, ontology_metadata = load_go_ancestors(config.ontology_resource)
        provenance.update(ontology_metadata)
        provenance["ontology_source_url"] = config.ontology_source_url
    return provenance


def _normalise_ranking(ranking: pd.Series, config: GeneSetEnrichmentConfig) -> pd.Series:
    values = pd.Series(ranking, dtype=float).copy()
    values.index = values.index.astype(str).str.strip()
    values = values.loc[values.index != ""]
    values = values.replace([np.inf, -np.inf], np.nan).dropna()
    if values.index.has_duplicates:
        duplicates = values.index[values.index.duplicated()].unique().tolist()
        raise ValueError(f"ranking contains duplicate gene identifiers: {duplicates[:5]}")
    if len(values) < config.minimum_ranked_genes:
        raise ValueError(
            f"ranking has {len(values)} finite genes; need at least {config.minimum_ranked_genes}"
        )
    return values


def run_preranked_gsea(
    rankings: Mapping[str, pd.Series], *, config: GeneSetEnrichmentConfig = DEFAULT_GO_BP_CONFIG,
    gene_sets: Mapping[str, tuple[str, ...]] | None = None,
    calculate_family_bh: bool = False,
) -> GeneSetEnrichmentRun:
    """Run one independent preranked GSEA per named signed ranking.

    Significance is the GSEA implementation's FDR q-value within each ranking.
    ``family_bh`` is optional descriptive metadata across all ranking-term pairs;
    it is never used by the specificity summary.
    """
    resolved_gene_sets = dict(gene_sets) if gene_sets is not None else load_gene_set_library(
        config.gene_set_resource
    )
    provenance = enrichment_provenance(
        config, resolved_gene_sets if gene_sets is not None else None,
    )
    rows = []
    failures = {}
    try:
        import gseapy
    except ImportError as exc:
        return GeneSetEnrichmentRun(
            pd.DataFrame(), resolved_gene_sets, provenance, {"gseapy": str(exc)},
        )

    for analysis_set, raw_ranking in rankings.items():
        label = str(analysis_set)
        try:
            ranking = _normalise_ranking(raw_ranking, config)
            ranked_genes = ranking.rename_axis("gene").reset_index(name="score")
            prerank = gseapy.prerank(
                rnk=ranked_genes,
                gene_sets=resolved_gene_sets,
                min_size=config.minimum_gene_set_size,
                max_size=config.maximum_gene_set_size,
                permutation_num=config.permutations,
                seed=config.seed,
                threads=1,
                outdir=None,
                no_plot=True,
                verbose=False,
            )
            if prerank.res2d is None or prerank.res2d.empty:
                failures[label] = "No gene sets passed the ranking-specific size and overlap filters"
                continue
            result = prerank.res2d.rename(columns={
                "Term": "gene_set",
                "ES": "enrichment_score",
                "NES": "normalized_enrichment_score",
                "NOM p-val": "nominal_pvalue",
                "FDR q-val": "analysis_set_fdr",
                "Lead_genes": "leading_edge_genes",
            }).copy()
            result.insert(0, "analysis_set", label)
            result.insert(1, "rank_metric", str(ranking.name or "signed_score"))
            result.insert(2, "n_ranked_genes", len(ranking))
            rows.append(result)
        except Exception as exc:  # noqa: BLE001 - GSEA failures are isolated per ranking.
            failures[label] = f"{type(exc).__name__}: {exc}"

    if not rows:
        return GeneSetEnrichmentRun(
            pd.DataFrame(), resolved_gene_sets, provenance, failures,
        )
    results = pd.concat(rows, ignore_index=True)
    for column in ("nominal_pvalue", "analysis_set_fdr", "normalized_enrichment_score"):
        results[column] = pd.to_numeric(results[column], errors="coerce")
    if calculate_family_bh:
        results["family_bh_nominal_pvalue"] = np.nan
        valid = results["nominal_pvalue"].notna()
        if valid.any():
            results.loc[valid, "family_bh_nominal_pvalue"] = multipletests(
                results.loc[valid, "nominal_pvalue"], method="fdr_bh",
            )[1]
    results = results.sort_values(
        ["analysis_set", "analysis_set_fdr", "nominal_pvalue", "gene_set"],
        na_position="last",
    ).reset_index(drop=True)
    return GeneSetEnrichmentRun(
        results, resolved_gene_sets, provenance, failures,
    )


def signed_de_gene_ranking(results: pd.DataFrame) -> pd.Series:
    """Build a signed DE ranking, preferring a complete model Wald statistic."""
    data = results.loc[results["gene"].notna()].copy()
    data["gene"] = data["gene"].astype(str)
    if "stat" in data:
        statistic = pd.to_numeric(data["stat"], errors="coerce")
        if statistic.notna().all() and not statistic.empty:
            return pd.Series(statistic.to_numpy(), index=data["gene"], name="wald_statistic")

    if not {"log2FoldChange", "pvalue"}.issubset(data.columns):
        return pd.Series(dtype=float)
    effect = pd.to_numeric(data["log2FoldChange"], errors="coerce")
    pvalue = pd.to_numeric(data["pvalue"], errors="coerce")
    signed_log_p = np.sign(effect) * -np.log10(pvalue.clip(
        lower=np.finfo(float).tiny, upper=1,
    ))
    return pd.Series(signed_log_p.to_numpy(), index=data["gene"], name="signed_minus_log10_p")


def _go_id(value: object) -> str | None:
    match = re.search(r"(GO:\d{7})", str(value))
    return match.group(1) if match else None


def leading_edge_gene_preview(value: object, *, max_genes: int = 12) -> tuple[int, str]:
    """Return the complete leading-edge count and a compact ranked gene preview."""
    genes = [gene.strip() for gene in re.split(r"[;,]", str(value)) if gene.strip()]
    preview = ", ".join(genes[:max_genes])
    if len(genes) > max_genes:
        preview += ", …"
    return len(genes), preview


def go_specific_gsea_summary(
    gsea_results: pd.DataFrame,
    go_ancestors: Mapping[str, frozenset[str]],
    *,
    group_column: str = "analysis_set",
    fdr_column: str = "analysis_set_fdr",
    max_terms: int = 10,
    fdr_threshold: float = 0.05,
    leading_edge_jaccard_threshold: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarise significant GO terms without hiding a displayed representative."""
    audit = gsea_results.copy()
    required = {group_column, fdr_column, "gene_set", "leading_edge_genes", "normalized_enrichment_score"}
    missing = required.difference(audit.columns)
    if missing:
        raise ValueError(f"GSEA results are missing columns required for summary: {sorted(missing)}")
    audit[fdr_column] = pd.to_numeric(audit[fdr_column], errors="coerce")
    audit["go_id"] = audit["gene_set"].map(_go_id)
    audit["specific_summary_status"] = "not_significant"
    audit["represented_by"] = ""
    audit["leading_edge_jaccard"] = np.nan
    significant = audit.loc[audit[fdr_column].le(fdr_threshold)].copy()
    audit.loc[significant.index, "specific_summary_status"] = "selected"

    visible_indices: set[int] = set()
    collapsed_indices: dict[int, tuple[int, float]] = {}
    for _, group_rows in significant.groupby(group_column, observed=True, sort=False):
        indices = group_rows.index.tolist()
        edges = {
            index: {gene.strip() for gene in re.split(r"[;,]", str(row["leading_edge_genes"])) if gene.strip()}
            for index, row in group_rows.iterrows()
        }
        candidates: dict[int, list[tuple[int, float]]] = {}
        for parent_index, parent in group_rows.iterrows():
            parent_go = parent["go_id"]
            parent_edges = edges[parent_index]
            parent_sign = np.sign(float(parent["normalized_enrichment_score"]))
            if not parent_go or not parent_edges or not parent_sign:
                continue
            supported = []
            for child_index in indices:
                child = group_rows.loc[child_index]
                child_go = child["go_id"]
                child_edges = edges[child_index]
                child_sign = np.sign(float(child["normalized_enrichment_score"]))
                if (
                    child_index == parent_index or not child_go or parent_go not in go_ancestors.get(child_go, ())
                    or child_sign != parent_sign or not child_edges
                ):
                    continue
                overlap = len(parent_edges & child_edges) / len(parent_edges | child_edges)
                if overlap >= leading_edge_jaccard_threshold:
                    supported.append((child_index, overlap))
            if supported:
                candidates[parent_index] = sorted(
                    supported,
                    key=lambda item: (
                        len(go_ancestors.get(str(group_rows.loc[item[0], "go_id"]), ())),
                        -float(group_rows.loc[item[0], fdr_column]),
                    ),
                    reverse=True,
                )

        def representative(
            index: int, active: frozenset[int] = frozenset(), candidates=candidates,
        ) -> tuple[int, float]:
            if index in active or index not in candidates:
                return index, 1.0
            child_index, overlap = candidates[index][0]
            terminal, _ = representative(child_index, active | {index})
            return terminal, overlap

        groups: dict[int, list[int]] = {}
        for index in indices:
            terminal, overlap = representative(index)
            groups.setdefault(terminal, []).append(index)
            if terminal != index:
                collapsed_indices[index] = (terminal, overlap)
        priority = {
            terminal: min(float(group_rows.loc[index, fdr_column]) for index in members)
            for terminal, members in groups.items()
        }
        visible_indices.update(sorted(
            groups, key=lambda index: (priority[index], float(group_rows.loc[index, fdr_column]))
        )[:max_terms])

    for parent_index, (representative_index, overlap) in collapsed_indices.items():
        if representative_index in visible_indices:
            audit.loc[parent_index, "specific_summary_status"] = "collapsed_parent"
            audit.loc[parent_index, "represented_by"] = str(audit.loc[representative_index, "gene_set"])
            audit.loc[parent_index, "leading_edge_jaccard"] = overlap
    for index in significant.index:
        if index not in visible_indices and audit.loc[index, "specific_summary_status"] == "selected":
            audit.loc[index, "specific_summary_status"] = f"outside_top_{max_terms}"
    summary = audit.loc[sorted(visible_indices)].copy()
    summary = summary.sort_values(
        [group_column, fdr_column, "gene_set"], na_position="last"
    ).reset_index(drop=True)
    return summary, audit


def save_go_specificity_outputs(
    run: GeneSetEnrichmentRun,
    *,
    output_dir: Path,
    prefix: str,
    group_column: str = "analysis_set",
    fdr_column: str = "analysis_set_fdr",
    config: GeneSetEnrichmentConfig = DEFAULT_GO_BP_CONFIG,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Save a GO-DAG-aware summary, audit, and complete reproducibility manifest."""
    if not config.ontology_resource:
        raise ValueError("GO-specific summaries require an ontology resource")
    go_ancestors, _ = load_go_ancestors(config.ontology_resource)
    summary, audit = go_specific_gsea_summary(
        run.results,
        go_ancestors,
        group_column=group_column,
        fdr_column=fdr_column,
        fdr_threshold=config.native_fdr_threshold,
        leading_edge_jaccard_threshold=config.leading_edge_jaccard_threshold,
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / f"{prefix}_specific_summary.csv", index=False)
    audit.to_csv(output_dir / f"{prefix}_specificity_audit.csv", index=False)
    manifest = {**run.provenance, "failures": run.failures}
    (output_dir / f"{prefix}_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary, audit

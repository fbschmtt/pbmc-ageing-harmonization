import gzip
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pbmc_pipeline import go_enrichment
from pbmc_pipeline.go_enrichment import (
    DEFAULT_GO_BP_CONFIG,
    GeneSetEnrichmentRun,
    go_specific_gsea_summary,
    leading_edge_gene_preview,
    run_preranked_gsea,
    save_go_specificity_outputs,
    signed_de_gene_ranking,
)


@pytest.fixture(autouse=True)
def downloaded_go_inputs(tmp_path, monkeypatch) -> None:
    """Provide the two downloaded GO inputs required by default enrichment."""
    data_dir = tmp_path / "go"
    data_dir.mkdir()
    (data_dir / "go_biological_process_2023.gmt").write_text(
        "Fixture process (GO:0000001)\t\tG0\tG1\tG2\n"
    )
    with gzip.open(data_dir / "go-basic.obo.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "format-version: 1.2\n"
            "data-version: fixture-1\n\n"
            "[Term]\n"
            "id: GO:0000001\n"
            "name: Fixture process\n"
        )
    monkeypatch.setenv(go_enrichment.GO_ENRICHMENT_DATA_DIR_ENV, str(data_dir))
    go_enrichment.load_gene_set_library.cache_clear()
    go_enrichment.load_go_ancestors.cache_clear()


def _gsea_rows() -> pd.DataFrame:
    return pd.DataFrame({
        "Term": ["Child process (GO:0000002)"],
        "ES": [0.4],
        "NES": [1.8],
        "NOM p-val": [0.001],
        "FDR q-val": [0.02],
        "Lead_genes": ["A;B;C"],
    })


def test_preranked_gsea_keeps_successful_rankings_when_one_ranking_fails(monkeypatch) -> None:
    def prerank(*, rnk, **_kwargs):
        if rnk["gene"].iloc[0] == "BAD":
            raise ValueError("bad ranking")
        return SimpleNamespace(res2d=_gsea_rows())

    monkeypatch.setitem(sys.modules, "gseapy", SimpleNamespace(prerank=prerank))
    ranking = pd.Series(np.arange(10, dtype=float), index=[f"G{i}" for i in range(10)])
    failed = ranking.rename(index={"G0": "BAD"})

    run = run_preranked_gsea({"good": ranking, "bad": failed})

    assert run.results["analysis_set"].tolist() == ["good"]
    assert "bad" in run.failures
    assert run.results["analysis_set_fdr"].tolist() == [0.02]
    assert "global_fdr" not in run.results
    assert run.provenance["gene_set_resource"] == "go_biological_process_2023.gmt"


def test_preranked_gsea_accepts_a_caller_provided_gene_set_collection(monkeypatch) -> None:
    captured = {}

    def prerank(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(res2d=_gsea_rows())

    monkeypatch.setitem(sys.modules, "gseapy", SimpleNamespace(prerank=prerank))
    ranking = pd.Series(np.arange(10, dtype=float), index=[f"G{i}" for i in range(10)])
    supplied_sets = {"A supplied term (GO:0000001)": ("G0", "G1", "G2")}

    run = run_preranked_gsea({"custom": ranking}, gene_sets=supplied_sets)

    assert captured["gene_sets"] == supplied_sets
    assert run.provenance["gene_set_resource"] == "caller-provided"


def test_preranked_gsea_rejects_duplicate_gene_identifiers_per_ranking(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "gseapy", SimpleNamespace(prerank=lambda **_: None))
    ranking = pd.Series(np.arange(10, dtype=float), index=["DUP", "DUP", *[f"G{i}" for i in range(8)]])

    run = run_preranked_gsea({"duplicate": ranking})

    assert run.results.empty
    assert "duplicate gene identifiers" in run.failures["duplicate"]


def test_de_ranking_uses_one_consistent_scale_when_wald_statistics_are_incomplete() -> None:
    results = pd.DataFrame({
        "gene": ["A", "B"],
        "stat": [3.0, np.nan],
        "log2FoldChange": [1.0, -1.0],
        "pvalue": [0.01, 0.001],
    })

    ranking = signed_de_gene_ranking(results)

    assert ranking.name == "signed_minus_log10_p"
    assert ranking.index.tolist() == ["A", "B"]
    assert ranking.tolist() == [2.0, -3.0]


def test_leading_edge_preview_keeps_count_and_marks_truncation() -> None:
    count, preview = leading_edge_gene_preview("A; B, C; D", max_genes=3)

    assert count == 4
    assert preview == "A, B, C, …"


def test_specific_summary_uses_visible_same_direction_representative_and_custom_fdr_column() -> None:
    results = pd.DataFrame({
        "module": ["blue", "blue", "blue"],
        "module_fdr": [0.001, 0.04, 0.02],
        "gene_set": [
            "Broad process (GO:0000001)",
            "Specific process (GO:0000002)",
            "Independent process (GO:0000003)",
        ],
        "leading_edge_genes": ["A;B;C", "A;B;C", "Q;R;S"],
        "normalized_enrichment_score": [2.0, 1.5, 1.7],
    })
    ancestors = {"GO:0000002": frozenset({"GO:0000001"})}

    summary, audit = go_specific_gsea_summary(
        results, ancestors, group_column="module", fdr_column="module_fdr", max_terms=1,
    )

    assert summary["gene_set"].tolist() == ["Specific process (GO:0000002)"]
    broad = audit.loc[audit["gene_set"].str.startswith("Broad")].iloc[0]
    assert broad["specific_summary_status"] == "collapsed_parent"
    assert broad["represented_by"] == "Specific process (GO:0000002)"


def test_specific_summary_keeps_opposite_direction_or_low_overlap_parent() -> None:
    results = pd.DataFrame({
        "analysis_set": ["age", "age", "age"],
        "analysis_set_fdr": [0.01, 0.02, 0.03],
        "gene_set": [
            "Opposing parent (GO:0000001)",
            "Opposing child (GO:0000002)",
            "Low-overlap child (GO:0000003)",
        ],
        "leading_edge_genes": ["A;B;C", "A;B;C", "X;Y;Z"],
        "normalized_enrichment_score": [2.0, -2.0, 2.0],
    })
    ancestors = {
        "GO:0000002": frozenset({"GO:0000001"}),
        "GO:0000003": frozenset({"GO:0000001"}),
    }

    summary, audit = go_specific_gsea_summary(results, ancestors)

    assert len(summary) == 3
    assert set(audit["specific_summary_status"]) == {"selected"}


def test_specificity_outputs_record_downloaded_resource_provenance(tmp_path) -> None:
    run = GeneSetEnrichmentRun(
        results=pd.DataFrame({
            "analysis_set": ["age"],
            "analysis_set_fdr": [0.01],
            "gene_set": ["A process (GO:0000001)"],
            "leading_edge_genes": ["A;B"],
            "normalized_enrichment_score": [2.0],
        }),
        gene_sets={},
        provenance={"gene_set_resource": DEFAULT_GO_BP_CONFIG.gene_set_resource},
        failures={"other": "no overlap"},
    )

    save_go_specificity_outputs(run, output_dir=tmp_path, prefix="age")

    provenance = (tmp_path / "age_provenance.json").read_text()
    assert "go_biological_process_2023.gmt" in provenance
    assert '"other": "no overlap"' in provenance

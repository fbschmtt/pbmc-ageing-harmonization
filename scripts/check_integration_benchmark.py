"""Check the integration-benchmark artifact published by the test workflow."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate integration benchmark test outputs")
    parser.add_argument("--outdir", type=Path, default=Path("output/test"))
    args = parser.parse_args()
    directory = args.outdir / "integration_benchmark"
    artifact_path = directory / "integration_benchmark.h5ad"
    report_path = directory / "integration_benchmark.json"
    html_path = directory / "report.html"
    notebook_path = directory / "executed.ipynb"
    if not all(path.is_file() for path in (artifact_path, report_path, html_path, notebook_path)):
        raise SystemExit(f"Missing integration benchmark artifacts under {directory}")

    artifact = ad.read_h5ad(artifact_path)
    expected_methods = {
        "harmony", "unintegrated_pca", "scanorama", "bbknn", "scvi", "scvi_unintegrated",
    }
    report = json.loads(report_path.read_text())
    method_records = report.get("methods", [])
    methods = {method["name"] for method in method_records}
    if methods != expected_methods:
        raise SystemExit(f"Integration benchmark has unexpected methods: {sorted(methods)}")
    required_columns = {"aifi_l2_study_majority"} | {
        method["label_column"] for method in method_records
    }
    missing = required_columns - set(artifact.obs)
    if artifact.n_obs == 0 or artifact.n_vars != 0 or missing:
        raise SystemExit(
            "Invalid thin integration benchmark artifact: "
            f"shape={artifact.shape}, missing_columns={sorted(missing)}"
        )
    expected_embeddings = {method["embedding_key"] for method in method_records}
    expected_embeddings |= {
        method["representation_key"]
        for method in method_records
        if method.get("representation_key")
    }
    if set(artifact.obsm) != expected_embeddings or any(
        artifact.obsm[method["embedding_key"]].shape != (artifact.n_obs, 2)
        for method in method_records
    ):
        raise SystemExit("Integration benchmark artifact is missing its method embeddings")
    expected_graphs = {
        key for method in method_records for key in method["graph_keys"].values()
    }
    if set(artifact.obsp) != expected_graphs:
        raise SystemExit(f"Integration benchmark has unexpected graph keys: {sorted(artifact.obsp)}")
    metadata = artifact.uns.get("integration_benchmark", {})
    if metadata.get("source_matrix_retained") is not False:
        raise SystemExit("Integration benchmark artifact unexpectedly retains its source matrix")
    if set(map(str, metadata.get("method_names", []))) != expected_methods:
        raise SystemExit("Integration benchmark artifact has unexpected method names")

    if report.get("status") != "ok" or not report.get("method_graphs_reused_for_celltypist_and_umap"):
        raise SystemExit("Integration benchmark JSON does not record graph reuse for labels and UMAPs")
    html = html_path.read_text()
    if "<title>Integration benchmark | PBMC ageing</title>" not in html:
        raise SystemExit("Integration benchmark report has an uninformative browser title")
    if "Embedding integration" not in html or "Type calling" not in html:
        raise SystemExit("Integration benchmark report lacks embedding or type-calling sections")
    print(f"Validated integration benchmark: {artifact.n_obs:,} cells, {len(methods)} methods")


if __name__ == "__main__":
    main()

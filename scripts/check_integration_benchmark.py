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
    required_columns = {
        "aifi_l2_study_majority",
        "benchmark_harmony_aifi_l2_majority",
        "benchmark_unintegrated_pca_aifi_l2_majority",
    }
    missing = required_columns - set(artifact.obs)
    if artifact.n_obs == 0 or artifact.n_vars != 0 or missing:
        raise SystemExit(
            "Invalid thin integration benchmark artifact: "
            f"shape={artifact.shape}, missing_columns={sorted(missing)}"
        )
    if set(artifact.obsm) != {"X_umap"} or artifact.obsm["X_umap"].shape != (artifact.n_obs, 2):
        raise SystemExit("Integration benchmark must retain exactly one two-dimensional UMAP")
    if artifact.obsp:
        raise SystemExit("Integration benchmark must not retain neighbor graphs")
    metadata = artifact.uns.get("integration_benchmark", {})
    if metadata.get("source_matrix_retained") is not False:
        raise SystemExit("Integration benchmark artifact unexpectedly retains its source matrix")
    if set(map(str, metadata.get("method_names", []))) != {"harmony", "unintegrated_pca"}:
        raise SystemExit("Integration benchmark artifact has unexpected method names")

    report = json.loads(report_path.read_text())
    methods = {method["name"] for method in report.get("methods", [])}
    if report.get("status") != "ok" or not report.get("harmony_graph_reused_for_umap"):
        raise SystemExit("Integration benchmark JSON does not record a successful reused Harmony graph")
    if methods != {"harmony", "unintegrated_pca"}:
        raise SystemExit(f"Integration benchmark JSON has unexpected methods: {sorted(methods)}")
    if "AIFI L2 label concordance" not in html_path.read_text():
        raise SystemExit("Integration benchmark report lacks its concordance section")
    print(f"Validated integration benchmark: {artifact.n_obs:,} cells, {len(methods)} methods")


if __name__ == "__main__":
    main()

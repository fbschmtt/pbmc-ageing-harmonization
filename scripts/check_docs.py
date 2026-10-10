"""Check local Markdown links and documented Make targets without extra tooling."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = (
    ROOT / "AGENTS.md",
    ROOT / "README.md",
    ROOT / "IMPLEMENTATION.md",
    ROOT / "INPUT_FILES.md",
    ROOT / "PLAN.md",
    ROOT / "input_data/README.md",
    ROOT / "aifi_models/README",
)
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")
MAKE_COMMAND = re.compile(r"^\s*make\s+([A-Za-z0-9_-]+)\b", re.MULTILINE)
MAKE_TARGET = re.compile(r"^([A-Za-z0-9_-]+):", re.MULTILINE)
AGE_SETTINGS_START = "<!-- generated-age-trajectory-settings:start -->"
AGE_SETTINGS_END = "<!-- generated-age-trajectory-settings:end -->"
AGE_SETTING_MEANINGS = {
    "bin_width_years": "Width of the half-open age bins, such as `[20,30)`.",
    "minimum_samples_per_bin": "Minimum eligible sample-by-cell-type pseudobulks for a bin to be retained.",
    "minimum_bins": "Minimum retained bins for a cell type to enter trajectory displays and cross-type analyses.",
    "reference_bin_start_age": "Start age of the zero-reference bin; this is `[20,30)` with the current width.",
    "strict_age_cutoff_exclusive": "Samples with age ≥90 are excluded before trajectory binning and residual fitting.",
    "cross_cell_type_residual_minimum_sample_coverage": "Target fraction of the union of fitted samples present in every cell type retained for the combined residual embedding.",
    "cluster_fdr_threshold": "Omnibus adjusted-p-value cutoff for trajectories entering clustering.",
    "de_fdr_threshold": "Omnibus adjusted-p-value cutoff for cross-type significance and recurrence summaries.",
    "max_clusters": "Upper bound on trajectory hierarchical clusters and candidate residual K-means groups.",
    "linkage_method": "Hierarchical clustering linkage. Supported values are `single`, `complete`, `average`, `weighted`, `centroid`, `median`, and `ward`.",
    "distance_metric": "Profile distance. Supported values are `euclidean`, `cityblock`, `cosine`, and `correlation`; Ward, centroid, and median require Euclidean distance.",
    "umap_neighbors": "UMAP neighborhood size.",
    "minimum_umap_trajectories": "Minimum number of profiles required to compute UMAP coordinates.",
    "umap_min_dist": "UMAP minimum-distance parameter.",
    "report_top_n_genes": "Maximum gene rows shown in a per-cell-type trajectory table.",
    "cross_report_top_n_genes": "Maximum gene rows shown in cross-cell-type trajectory tables.",
    "random_state": "Seed for reproducible UMAP coordinates.",
    "minimum_shared_bins": "Minimum common age bins required for the cross-cell-type trajectory comparison.",
    "residual_pca_components": "Maximum number of centered PCA components used before residual UMAP. The count is capped by the available samples and genes.",
    "residual_umap_neighbors": "Neighbourhood size for residual-sample UMAPs; capped at one fewer than the available samples.",
    "residual_hdbscan_min_cluster_size": "Smallest sample group retained by the HDBSCAN alternatives; points outside selected clusters are labeled noise.",
    "residual_hdbscan_min_samples": "Local-density conservativeness, used for both HDBSCAN comparison embeddings.",
    "residual_hdbscan_cluster_selection_method": "HDBSCAN selection strategy, used for both UMAP and PC-space comparison labels. No HDBSCAN parameter sweep is shown.",
    "residual_pc_hdbscan_components": "Maximum number of leading residual PCs used for the HDBSCAN-in-PC-space alternative and primary Leiden graph.",
    "residual_leiden_neighbors": "Number of PC-space nearest neighbors used to construct the Leiden graph; capped at one fewer than the available samples.",
    "residual_leiden_resolutions": "Candidate Leiden resolutions tested on the same PC-space graph.",
    "residual_leiden_min_clusters": "Target used to select the lowest candidate resolution yielding at least this many clusters; if none do, use the candidate yielding the most clusters.",
    "residual_marker_max_clusters": "Maximum number of largest clusters included in marker and feature-origin contrasts.",
}


def _format_setting_default(value: object) -> str:
    if isinstance(value, (str, list)):
        rendered = json.dumps(value, ensure_ascii=False, separators=(", ", ": "))
        return f"`{rendered}`"
    if isinstance(value, bool):
        return f"`{str(value).lower()}`"
    return str(value)


def generated_age_trajectory_table() -> str:
    pipeline = json.loads((ROOT / "config/pipeline.json").read_text(encoding="utf-8"))
    settings = pipeline["differential_expression"]["age_trajectory"]
    missing = set(settings) - set(AGE_SETTING_MEANINGS)
    stale = set(AGE_SETTING_MEANINGS) - set(settings)
    if missing or stale:
        raise ValueError(
            "Age-trajectory documentation descriptions and config keys differ "
            f"(missing descriptions={sorted(missing)}, stale descriptions={sorted(stale)})"
        )
    lines = ["| Setting | Default | Meaning |", "| --- | ---: | --- |"]
    for name, value in settings.items():
        lines.append(
            f"| `{name}` | {_format_setting_default(value)} | {AGE_SETTING_MEANINGS[name]} |"
        )
    return "\n".join(lines)


def sync_age_trajectory_table(*, write: bool) -> list[str]:
    document = ROOT / "IMPLEMENTATION.md"
    contents = document.read_text(encoding="utf-8")
    if AGE_SETTINGS_START not in contents or AGE_SETTINGS_END not in contents:
        return ["IMPLEMENTATION.md: generated age-trajectory settings markers are missing"]
    start = contents.index(AGE_SETTINGS_START) + len(AGE_SETTINGS_START)
    end = contents.index(AGE_SETTINGS_END)
    rendered = "\n\n" + generated_age_trajectory_table() + "\n"
    if contents[start:end] == rendered:
        return []
    if write:
        document.write_text(contents[:start] + rendered + contents[end:], encoding="utf-8")
        return []
    return [
        "IMPLEMENTATION.md: age-trajectory settings table is stale; run `make docs-update`"
    ]


def local_link_errors(document: Path) -> list[str]:
    errors: list[str] = []
    for destination in MARKDOWN_LINK.findall(document.read_text(encoding="utf-8")):
        destination = destination.strip().strip("<>")
        path, _, _fragment = destination.partition("#")
        if not path or "://" in path or path.startswith(("mailto:", "/")):
            continue
        if not (document.parent / path).exists():
            errors.append(f"{document.relative_to(ROOT)}: missing link target {path}")
    return errors


def documented_make_target_errors() -> list[str]:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    targets = set(MAKE_TARGET.findall(makefile))
    errors: list[str] = []
    for document in DOCUMENTS:
        for target in MAKE_COMMAND.findall(document.read_text(encoding="utf-8")):
            if target not in targets:
                errors.append(
                    f"{document.relative_to(ROOT)}: documented Make target does not exist: {target}"
                )
    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true",
        help="update generated documentation from config/pipeline.json",
    )
    args = parser.parse_args()
    try:
        generated_errors = sync_age_trajectory_table(write=args.write)
    except (KeyError, ValueError) as error:
        print(f"Documentation generation failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    errors = [error for document in DOCUMENTS for error in local_link_errors(document)]
    errors.extend(documented_make_target_errors())
    errors.extend(generated_errors)
    if errors:
        print("Documentation checks failed:", *errors, sep="\n", file=sys.stderr)
        raise SystemExit(1)
    print(f"Documentation checks passed for {len(DOCUMENTS)} files.")


if __name__ == "__main__":
    main()

"""Check local Markdown links and documented Make targets without extra tooling."""
from __future__ import annotations

import argparse
import json
import re
import shlex
import sys
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")
MAKE_TARGET = re.compile(r"^([A-Za-z0-9_-]+(?:[ \t]+[A-Za-z0-9_-]+)*):", re.MULTILINE)
SETTINGS_DOCUMENT = Path("docs/reference/analysis-settings.md")
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
    "residual_leiden_min_clusters": "Target cluster count; select the candidate resolution whose cluster count is closest to this value, with ties favoring fewer clusters.",
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
    document = ROOT / SETTINGS_DOCUMENT
    contents = document.read_text(encoding="utf-8")
    if AGE_SETTINGS_START not in contents or AGE_SETTINGS_END not in contents:
        return [f"{SETTINGS_DOCUMENT}: generated age-trajectory settings markers are missing"]
    start = contents.index(AGE_SETTINGS_START) + len(AGE_SETTINGS_START)
    end = contents.index(AGE_SETTINGS_END)
    rendered = "\n\n" + generated_age_trajectory_table() + "\n"
    if contents[start:end] == rendered:
        return []
    if write:
        document.write_text(contents[:start] + rendered + contents[end:], encoding="utf-8")
        return []
    return [
        f"{SETTINGS_DOCUMENT}: age-trajectory settings table is stale; run `make docs-update`"
    ]


def documentation_files(root: Path = ROOT) -> list[Path]:
    """Discover guides, reference pages and repository entry points."""
    paths = {*root.glob("*.md"), *(root / "docs").rglob("*.md")}
    paths.update((
        root / "input_data/README.md", root / "site/DATA_CONTRACT.md", root / "aifi_models/README",
    ))
    return sorted(path for path in paths if path.is_file())


def prose_lines(contents: str):
    """Exclude fenced examples when interpreting headings and local links."""
    fence = None
    for line in contents.splitlines():
        if match := re.match(r"^\s*(`{3,}|~{3,})", line):
            marker = match.group(1)[0]
            fence = None if fence == marker else marker if fence is None else fence
            continue
        if fence is None:
            yield line


def heading_anchors(contents: str) -> set[str]:
    anchors = set()
    counts: dict[str, int] = {}
    for line in prose_lines(contents):
        if heading := re.match(r"^#{1,6}\s+(.+?)\s*#*$", line):
            text = re.sub(r"\[([^]]+)\]\([^)]*\)", r"\1", heading.group(1)).replace("`", "")
            base = re.sub(r"[^\w\- ]", "", text.lower()).replace(" ", "-")
            count = counts.get(base, 0)
            anchors.add(base + (f"-{count}" if count else ""))
            counts[base] = count + 1
        anchors.update(re.findall(r'<a\s+[^>]*(?:id|name)=["\x27]([^"\x27]+)', line))
    return anchors


def local_link_errors(document: Path, root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    contents = "\n".join(prose_lines(document.read_text(encoding="utf-8")))
    for destination in MARKDOWN_LINK.findall(contents):
        destination = destination.strip().strip("<>")
        path, _, fragment = destination.partition("#")
        if "://" in path or path.startswith(("mailto:", "/")):
            continue
        target = document.parent / unquote(path) if path else document
        if not target.exists():
            errors.append(f"{document.relative_to(root)}: missing link target {path}")
        elif (fragment and target.is_file() and (target.suffix == ".md" or target.name == "README")
              and unquote(fragment) not in heading_anchors(target.read_text(encoding="utf-8"))):
            errors.append(f"{document.relative_to(root)}: missing heading {destination}")
    return errors


def documented_make_target_errors(root: Path = ROOT) -> list[str]:
    makefile = "\n".join(
        path.read_text(encoding="utf-8")
        for path in [root / "Makefile", *sorted((root / "make").glob("*.mk"))]
    )
    targets = {target for group in MAKE_TARGET.findall(makefile) for target in group.split()}
    errors: list[str] = []
    for document in documentation_files(root):
        contents = document.read_text(encoding="utf-8")
        commands = re.findall(r"^\s*(make\s+[^\n]+)", contents, re.MULTILINE)
        commands += re.findall(r"`(make\s+[^`\n]+)`", contents)
        for command in commands:
            try:
                words = shlex.split(command)
            except ValueError:
                continue  # A continued shell example is checked one line at a time.
            for target in words[1:]:
                if target.startswith("#"):
                    break
                if "=" in target or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", target):
                    continue
                if target in targets:
                    continue
                errors.append(
                    f"{document.relative_to(root)}: documented Make target does not exist: {target}"
                )
    return sorted(set(errors))


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
    documents = documentation_files()
    errors = [error for document in documents for error in local_link_errors(document)]
    errors.extend(documented_make_target_errors())
    errors.extend(generated_errors)
    if errors:
        print("Documentation checks failed:", *errors, sep="\n", file=sys.stderr)
        raise SystemExit(1)
    print(f"Documentation checks passed for {len(documents)} files.")


if __name__ == "__main__":
    main()

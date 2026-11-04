# Static atlas data contract

`make site-export` reads existing pipeline outputs; it does not run or resume
harmonization or any analysis. It generates schema version 1 assets in a fresh
staging directory, verifies them, and then replaces `site/public/data`. Failure
before promotion leaves the previous data directory in place. Promotion uses a
rename with rollback, not an atomic update of an already hosted site. Publish a
whole validated build as an immutable version and switch the host to it as one
operation. Configure cache invalidation for `data/release.json` and the site
entry point. The browser checks the release ID of each asset against that
manifest; a mixed deployment is rejected rather than silently combined.

`release.json` records the export ID and time, exporter Git revision, hashes of
observed analysis summary inputs, output paths, sizes and SHA-256 hashes,
availability, per-type model status, and input/exported row counts. When the
combined DE manifest exists, selected analysis metadata is recorded under
`analysis`, including the configuration digest, method and inclusion rules. The
Git revision is the exporter checkout revision, not proof of the analysis code
revision. The header shows the recorded analysis mode and export date; the export
ID and date identify the site export, not an analysis version. Input hashes establish what was read but do not prove that files came
from one run or that unversioned sources were authorised. When model-level
records are available in the combined DE manifest, export compares their cell
type, mode, model metadata and result paths against each per-type run metadata
file. Contradictions fail export; unavailable or incomplete comparisons are
recorded in `provenance`. A consistent comparison does not prove table contents
or report lineage. Before publication,
compare the manifest against the intended analysis run and document cohorts,
source citations and licences, annotation decisions, age support and any
omissions. Preserve the release manifest and the originating workflow outputs.

The exporter inventories every cell type with DE metadata, including skipped
models and reasons when supplied. If a combined DE manifest is available, each
listed type must have metadata. Missing paired trajectory inputs, invalid
required columns, missing completed-model result files, invalid generated asset
schemas and a nonempty trajectory input yielding no valid records fail export.
Optional views can be absent; `availability` records their presence. A row can
be excluded when its required values are missing or invalid, or its cell type
has no exported DE model. `counts` records input and retained totals and groups excluded rows by reason;
export prints a warning when these differ. Review these counts and model statuses
before release; the counts do not identify individual excluded rows.

All published JSON data assets have `schemaVersion: 1` and the same `releaseId`.
The shared executable structural and range checks are in
[`src/data/contracts.mjs`](src/data/contracts.mjs). The release verifier checks
those contracts and hashes; browser loaders check the contract and release ID.
Unsupported versions are rejected. The internal `age-previews/manifest.json` is
a local cleanup list and is not a browser data asset.

| Asset | Meaning and key fields |
|---|---|
| `de-index.json` and `de-models/*.json` | The index lists models and points to a default chunk containing adjusted p < 0.05 plus gene lookup chunks for all other records. A complete gene symbol retrieves its lookup chunk. `effect` is log2 fold change; `padj` and `pvalue` are probabilities, `se` is effect standard error, and `baseMean` is PyDESeq2’s mean size-factor-normalized count. `studySupport`, when available, counts per-study results for the same covariate and contrast with adjusted p < 0.05 and records the number with a finite adjusted p-value. The count ignores effect direction; studies with missing results or adjusted p-values are not counted as tested. It is descriptive, and per-study models can have different covariate adjustments. Null numeric values mean no finite estimate was supplied. |
| `age-previews/<slug>.json` and gene chunks | Age-bin effects from the separate combined age-bin model are partitioned by gene. `effect` is log2 fold change relative to `referenceBin`; null denotes a missing estimate. |
| `pipeline-trajectory-index.json` | Only gene–cell-type trajectories selected for clustering. Each record has UMAP `x`/`y`, cluster, omnibus adjusted p-value, seven scaled age-bin values, and the signed maximum absolute **raw** log2 fold change relative to 20–30. When available, reference expression is study-balanced reference-bin CPM with contributing sample and study counts. Cell count is unavailable. This is not the universe of tested genes. |
| `pipeline-composition-index.json` | Study-level `records` describe mean adult fractions within an AIFI L1 parent and model diagnostics. `coefficients` are per-study effects with confidence bounds; age is per decade. `ageEstimates` are adjusted fraction changes relative to `referenceBin` for `l1_parent` or `retained_pbmc`, with confidence bounds and per-bin sample counts. Reference mean and adjusted fractions can be null. `ageStatus` records model completion or skip reasons. Fractions measure relative representation, not absolute abundance. |
| `pathway-index.json` | Optional curated sets with identifier, label, description, uppercase genes and source filename. Source-file hashes are recorded in the release manifest; permissions and citations require separate review. |
| `single-cell-index.json` and `single-cell/*.json` | Optional per-type summaries from completed single-cell analyses. When a cell-type manifest exists, only types marked complete in it are included; otherwise available completed reports are used. Each type contains the report's local cluster UMAP, up to eight positive and negative markers per cluster ranked by log₂ fold change, and cluster-fraction figures when available. Newly rendered reports provide a per-cluster figure and pooled/per-study summaries: sample counts, mean observed fraction, and unadjusted linear slopes in percentage points per decade. Figures show sample-level age and fraction points without sample identifiers. Newly rendered reports include eligible samples with zero cells in a cluster; older reports omit them. Trends are descriptive and the pooled trend is study-unadjusted. The original per-cell data and identifiers are not exported. Input report, notebook and marker-table hashes are included in the source inventory; hashes do not prove common run lineage. |
| `gene-metadata.json` | Optional PANTHER annotations for retained trajectory genes, including source version, URL and source SHA-256. Rebuilt only when the pinned source and trajectory index are available. |
| `reports/<release-id>/analysis/…/report.html` | Available per-type DE and combined trajectory HTML reports, copied unchanged and linked from the modelling page. These reports may contain sample-level identifiers and diagnostics. The export does not include their linked CSV or other companion downloads; those links can be unavailable on the site. |
| `reports/<release-id>/harmonization/harmonization_qc/<study>/report.html` | Available study QC reports, discovered from the harmonization output and linked on the datasets page. These may contain sample-level information; review them before publication. |
| `images/<release-id>/studies/<study>-umap.png` | The standalone AIFI L2 UMAP from the study QC notebook, with a transparent background. Older QC notebooks without this figure do not produce a UMAP asset; rerender the study report to include it. |
| `images/<release-id>/studies/<study>-heatmap.png` | The study-label versus AIFI L2 comparison extracted from the study QC notebook, with a transparent background, when available. |

DE model chunks are fetched only for the currently selected cell type and
model. Gene-set prefetching is limited to that selected model. Older exports
with all model records embedded in `de-index.json` are filtered in the view,
but necessarily download the entire index; regenerate the export to enable
on-demand model loading.

The manifest includes study identifiers observed in exported models and selected
study harmonization reports. The dataset page indicates which studies have each
kind of result. If a workflow run manifest is available, only its selected studies
are used for harmonization reports, avoiding stale study outputs. Available
DOIs, annotation and technology descriptions come from the current checkout of
`config/studies.json`, whose fingerprint is recorded; this is not a verified
analysis-time snapshot. Optional `site_comments` strings in individual study
entries are displayed alongside that metadata. Data licences are recorded as unknown until reviewed.
Age-bin support is drawn from the trajectory analysis metadata when available.
Regenerating a study QC report updates its notebook and HTML; the site continues
to show the previous export until `make site-export` is run again. See the
[site setup and refresh instructions](../README.md#static-presentation-site).

HTML reports are included by default when present. The workflow report-index
renderer generates `reports/<release-id>/index.html` using only the reports
packaged with the site; it does not link other workflow outputs. The manifest and modelling
page list relative companion links that are unavailable in the exported site. Review their contents and sharing
permissions for the intended audience before deploying the site. The export hashes
the report files but cannot establish that they were rendered from the same run
as the structured summaries. HTML reports are separate files and load when opened.

JSON assets with worthwhile compression have gzip companions with a `.gz.bin`
suffix and their sizes and hashes recorded in the manifest. Browsers supporting
`DecompressionStream` fetch and decompress these bytes directly; other browsers
fetch the original JSON. The host must serve the companions unchanged and
include them in deployment. The manifest also has a gzip companion. HTTP
compression and caching for HTML, JavaScript, CSS and fallback JSON still
depend on hosting configuration.

The trajectory map starts with retained trajectory records having adjusted
omnibus p < 0.05. Explicit gene lookup also loads matching records above that
threshold from `trajectory-extras/`. The source contains only trajectories
selected for clustering; the lookup cannot supply genes absent from that source.
The selected pathway or entered gene set is prefetched in the background, with two workers per view;
large pathways can still fetch many chunks. A bounded in-memory cache shares
in-flight and recently used assets for the current page session. It is not a
persistent cache or a substitute for host-side HTTP caching.

PANTHER metadata and pathway indexes are rebuilt in staging when their source
inputs are available. They are not copied from a previous release. Check the
manifest's availability and source hashes; optional data can be absent on a
fresh checkout. `make site-build` verifies the manifest and every listed asset
before producing `site/dist`. The exporter retains warnings for test-only and
synthetic data; selecting inputs appropriate for the intended audience remains a
release decision.

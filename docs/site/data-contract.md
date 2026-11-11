# Static atlas data contract

[Documentation](../README.md)

`make site-export` reads existing pipeline outputs; it does not run or resume
harmonization or any analysis. It generates schema version 1 assets in a fresh
staging directory, verifies them, and then replaces `site/public/data`. Failure
before promotion leaves the previous data directory in place. Promotion uses a
rename with rollback, not an atomic update of an already hosted site. Publish a
whole validated build as an immutable version and switch the host to it as one
operation. Configure cache invalidation for `data/release.json` and the site
entry point. JSON loaders check each asset's release ID against the manifest.
Binary images use release-versioned URLs; CSV preparation checks the manifest
release ID and that the requested download is listed. These browser checks do
not replace publishing a consistent release: binary files are not checked
against their manifest SHA-256 hashes at runtime. `make site-build` verifies
all listed files and their hashes before building.

`release.json` records the export ID and time, exporter Git revision, hashes of
observed analysis summary inputs, output paths, sizes and SHA-256 hashes,
availability, per-type model status, and input/exported row counts. When the
combined DE manifest exists, selected analysis metadata is recorded under
`analysis`, including the configuration digest, method and inclusion rules. The
Git revision is the exporter checkout revision, not proof of the analysis code
revision. The header shows the recorded analysis mode and export date; the export
ID and date identify the site export, not an analysis version. Input hashes
establish what was read but do not prove that files came from one run or that
unversioned sources were authorised. When model-level records are available in the combined DE manifest, export compares their cell
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
required columns, missing completed-model result files needed by exported views
or study-support counts, invalid generated asset schemas and a nonempty
trajectory input yielding no valid records fail export. Optional views can be
absent; `availability` records selected view families, while the asset list and
per-view indexes identify the published data. A row can be excluded when its required values are missing or invalid, or its cell type
has no exported DE model. `counts` records input and retained totals for
trajectory and composition exports and groups excluded rows by reason; export
prints a warning when these differ. Review these counts and model statuses
before release; the counts do not identify individual excluded rows.

Single-cell views consume `analysis/site-assets.json` (schema version 1) and
its indexed PNG files. This manifest records the cell-type label, named UMAP
and optional age-DE images, fraction notes, explicit trend/zero-inclusion
flags, and per-cluster images and fit summaries. Every image has a SHA-256
checksum checked during export. Notebook headings, cell positions and source
comments are presentation details and are not used by this exporter. Older
complete reports without the manifest must be rerendered with
`make render-cell-type CELL_TYPE=<slug>` (or `make render-cell-type-test` for
test outputs) before export. `make site-export` reports this migration need
without replacing the previous release.

All published JSON data assets have `schemaVersion: 1` and the same `releaseId`.
The shared executable structural and range checks are in
[`src/data/contracts.mjs`](../../site/src/data/contracts.mjs). The release verifier checks
those contracts and hashes; browser loaders check the contract and release ID.
Unsupported versions are rejected. The internal `age-previews/manifest.json` is
a local cleanup list and is not a browser data asset.

| Asset | Meaning and key fields |
|---|---|
| `de-index.json` and `de-models/*.json` | The index lists combined models plus the developmental `juvenile_age` and `age_25_vs_5` results and points to a default chunk containing adjusted p < 0.05 plus gene lookup chunks for all other records. A complete gene symbol retrieves its lookup chunk. `effect` is log2 fold change; `padj` and `pvalue` are probabilities, `se` is effect standard error, and `baseMean` is PyDESeq2’s mean size-factor-normalized count. `studySupport`, when available, counts per-study results for the same covariate and contrast with adjusted p < 0.05 and records the number with a finite adjusted p-value. The count ignores effect direction; studies with missing results or adjusted p-values are not counted as tested. It is descriptive, and per-study models can have different covariate adjustments. Null numeric values mean no finite estimate was supplied. |
| `clock-data.json` and `clock-models/<cell-type>/<method>.json` | The index lists every available clock method for each cell type in baseline-to-more-advanced order, includes performance metrics and a shared age-axis domain, and points to one prediction asset per method. The browser fetches prediction assets only when a method is selected. Training predictions are leave-one-study-out; held-out predictions come from studies excluded from fitting and tuning. Cross-type mean methods include donors with predictions from at least two types and consistent study/age metadata. Published points contain study, age and predicted age (plus the number of cell types for means); donor identifiers and model coefficients are not exported. |
| `reference-expression/<cell-type>.json` | Study-balanced CPM for genes with reference-bin pseudobulk support, with the reference bin, method, and contributing sample and study counts. CPM is normalized expression, not an absolute molecule count. |
| `de-downloads/<cell-type>/<model>-<covariate>.csv.gz` | Gzip-compressed combined-model gene-level differential-expression results, with no uncompressed duplicate. The browser extracts the archive on demand and saves a plain CSV; browsers without `DecompressionStream` receive a compressed download to extract manually. Available statistic columns are `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. `#` metadata lines precede the CSV header with cell type, model, scope, analysis mode, covariate, contrast, design, sample counts, and any interpretation warning. Study and site sample-count maps identify the included studies and sites, so separate included-name lists are omitted. CSV does not standardize comments; programmatic readers must skip lines beginning with `#`. |
| `de-previews/<cell-type>.svg` | Static preview of the merged age-model volcano plot, shown while the interactive DE records load. |
| `age-previews/<slug>.json` and gene chunks | Age-bin effects from the separate combined age-bin model are partitioned by gene. `effect` is log2 fold change relative to `referenceBin`; null denotes a missing estimate. |
| `pipeline-trajectory-map.json` | Compact significant-point coordinates and labels for the global map. Full seven-bin profiles are fetched from `trajectory-details/<gene-bucket>.json` when a point or gene needs detail. |
| `pipeline-trajectory-preview.svg` | Static preview shown while the interactive global map loads. |
| `trajectory-details/<gene-bucket>.json` | The single published copy of full gene–cell-type trajectories selected for clustering. Each record has UMAP `x`/`y`, cluster, omnibus adjusted p-value, seven scaled age-bin values, and the signed maximum absolute **raw** log2 fold change relative to 20–30. When available, reference expression is study-balanced reference-bin CPM with contributing sample and study counts. Cell count is unavailable. This is not the universe of tested genes. |
| `pipeline-composition-index.json` | Study-level `records` describe mean adult fractions within an AIFI L1 parent and model diagnostics. `coefficients` are per-study effects with confidence bounds; age is per decade. `ageEstimates` are adjusted fraction changes relative to `referenceBin` for `l1_parent` or `retained_pbmc`, with confidence bounds and per-bin sample counts. Reference mean and adjusted fractions can be null. `ageStatus` records model completion or skip reasons. Fractions measure relative representation, not absolute abundance. |
| `pathway-index.json` | Optional curated sets with identifier, label, description, uppercase genes and source filename. Source-file hashes are recorded in the release manifest; permissions and citations require separate review. |
| `single-cell-index.json`, `single-cell/*.json` and `single-cell-images/<sha256>.png` | Optional per-type summaries from completed single-cell analyses. When a cell-type manifest exists, only types marked complete in it are included; otherwise available completed reports are used. Each type references the named local cluster UMAP in `analysis/site-assets.json`, an optional image of significant merged age-DE gene expression on that UMAP, up to ten positive and ten negative markers per cluster ranked by log₂ fold change, and cluster-fraction figures when available. The age-DE figure selects up to five positive and four negative age effects by adjusted p-value and exports only the rendered PNG, not per-cell expression values. Newly rendered reports provide a per-cluster figure and pooled/per-study summaries: sample counts, mean observed fraction, and unadjusted linear slopes in percentage points per decade. Figures show sample-level age and fraction points without sample identifiers. The manifest explicitly records whether eligible zero-count samples are included. Trends are descriptive and the pooled trend is study-unadjusted. Figures are separate PNG files referenced by the JSON, including per-cluster fraction figures. Identical optimized PNG bytes share one content-addressed file; the filename is the SHA-256 hash of those bytes. The original per-cell data and identifiers are not exported. Input report, notebook, marker-table, asset-manifest and PNG hashes are included in the source inventory; hashes do not prove common run lineage. |
| `gene-metadata.json` | Optional PANTHER annotations for retained trajectory genes, including source version, URL and source SHA-256. Built during export when the pinned source and trajectory inputs are available. Annotations currently cover retained trajectory genes, so a gene with a linear DE result can still have no function annotation. |
| `reports/<release-id>/analysis/…/report.html` | Available per-type DE and combined trajectory HTML reports, copied unchanged and linked from the modelling page. These reports may contain sample-level identifiers and diagnostics. The export does not include their linked CSV or other companion downloads; those links can be unavailable on the site. |
| `reports/<release-id>/harmonization/harmonization_qc/<study>/report.html` | Available study QC reports, discovered from the harmonization output and linked on the datasets page. These may contain sample-level information; review them before publication. |
| `images/<release-id>/studies/<study>-umap.png` | The standalone AIFI L2 UMAP from the study QC notebook, with a transparent background. Older QC notebooks without this figure do not produce a UMAP asset; rerender the study report to include it. |
| `images/<release-id>/studies/<study>-heatmap.png` | The study-label versus AIFI L2 comparison extracted from the study QC notebook, with a transparent background, when available. |

## Browser display semantics

Per-study DE chunks are used in staging to compute study-support counts, then removed
before publication. The developmental `adult_age` result is omitted from the site;
it remains in the analysis outputs. It comes from a separate developmental fit
and is not asserted to be numerically identical to the combined age result.
The release inventory still records source model statuses, including models that
are not published as interactive assets. These publication choices do not remove
or refit any analysis outputs.

The per-gene volcano uses the merged linear-age model's effect per decade on
x and its negative log10 adjusted p-value on y. The adjusted p-value comes from
the same model as the effect. Points require finite effect and adjusted p-value
estimates; requested genes are not restricted to adjusted p < 0.05. Point colors
identify cell types. The pathway heatmap uses the same linear-age effects per
decade and includes finite effects even when an adjusted p-value is unavailable.
It does not perform a pathway enrichment test.

Gene and pathway views obtain their cell types and results from the DE
catalogue, including types and genes without clustered trajectories. They do
not load the UMAP trajectory data. A selected gene's age-bin preview loads
`age-previews/` from the separate per-type age-bin model on demand. Those
previews are distinct from the scaled seven-bin profiles used to cluster and
position UMAP points. The trajectory UMAP detail card shows the signed maximum
raw age-bin log2 fold change relative to 20–30 immediately below the combined
linear age effect per decade. The maximum retains the sign of the bin with the
largest absolute change; it is not a per-decade slope or a maximum of the scaled
profile. The card's omnibus adjusted p-value belongs to the age-bin model.
Maximum age-bin change is not used by the linear-effect volcanoes or pathway
heatmap. Both volcano views use the same square sizing, capped
at the available width or 68% of the small viewport height, whichever is smaller.

The per-cell-type differential-expression view colors points by
`studySupport.significant`; zero is gray, and the color scale's upper bound is
the largest `tested` count among the plotted genes. Missing support remains distinct from zero. The
reference-expression switch uses the study-balanced CPM in
`reference-expression/<cell-type>.json`, transformed as log₂(CPM + 1), and
scales to the range across available genes in that asset, independent of the
current model's plotted subset. These color choices do not alter the published
values. Shared estimate formatting uses two decimal places for absolute values
of at least one and two significant digits below one; plot axes can use their
own tick formatting. Downloaded CSV statistics keep their source precision.

`clock-data.json` carries one age domain shared by the available cell types and
methods. Each method's summary metrics are in that index, while its predictions
are in a separate asset fetched when selected. This keeps the plot axes fixed
and avoids loading unselected prediction sets.

On startup, the browser loads the selected cell type's significant combined
age-model records and prefetches the corresponding assets for up to five other
types, ranked by contributing studies and sample count. A fresh session prefers
Naive CD8 T cell when available, otherwise the first available type. Other cell
types, model selections and exact-gene lookup chunks load on demand. The compact
trajectory map loads when its view opens; full profiles load when a point or
specific gene needs them. Older exports without the compact map fall back to
`pipeline-trajectory-index.json`.

The manifest includes study identifiers observed in exported models and selected
study harmonization reports. The dataset page indicates which studies have each
kind of result. If a workflow run manifest is available, only its selected studies
are used for harmonization reports, avoiding stale study outputs. Available
publication titles, formal citations, DOIs, annotation and technology descriptions
come from the current checkout of `config/studies.json`, whose fingerprint is
recorded; this is not a verified analysis-time snapshot. Donor and cell counts
are derived from the available merged single-cell and study metadata. Data
licences are recorded as unknown until reviewed.
Age-bin support is drawn from the trajectory analysis metadata when available.
Regenerating a study QC report updates its notebook and HTML; the site continues
to show the previous export until `make site-export` is run again. See the
[site setup and refresh instructions](README.md#static-presentation-site).

HTML reports are included by default when present. The workflow report-index
renderer generates `reports/<release-id>/index.html` using only the reports
packaged with the site; it does not link other workflow outputs. The manifest
and modelling page list relative companion links that are unavailable in the
exported site. Review their contents and sharing permissions for the intended
audience before deploying the site. The export hashes the report files but cannot establish that they were rendered from the same run
as the structured summaries. HTML reports are separate files and load when opened.

Standalone single-cell and study PNGs are automatically recompressed at zlib
level 9 during every export. The decompressed filtered scanlines and all
non-image PNG chunks are preserved, retaining the pixels, transparency,
dimensions and metadata.
The exporter compares one level-9 recompression with the original and publishes
the smaller file automatically, including in production exports. It does not
search multiple settings, resize figures or reduce color precision. The source
notebook and its images are not rewritten; animated PNGs are kept unchanged.
PNGs have no extra gzip copy. HTML reports remain unchanged, including their
embedded figures. Older releases with base64 images remain readable.

Combined-model CSV downloads are prepared from source result tables during
export and stored as `.csv.gz` files. They are not generated from model JSON
at download time. Downloads retain source statistics at their original scale
and precision; linear-age effects remain per year even when the plot shows
per-decade effects. Developmental results do not have CSV downloads.

CSV downloads are fetched only when requested and decompressed in the browser.
The download shows a preparation state and allows retry after a failed request
or invalid gzip stream; no file is saved until decompression finishes. Changing
the selected download or leaving the view cancels an in-progress request.
Download assets must be served unchanged, without HTTP `Content-Encoding`.
Legacy releases with plain CSV downloads remain supported.

JSON assets and SVG previews with worthwhile compression have raw gzip
companions with a `.gz.bin` suffix; their sizes and hashes are recorded in the
manifest. JSON loaders supporting `DecompressionStream` fetch and decompress
the compressed JSON directly; other browsers fetch the original JSON. SVG image
links use the original SVG files. The manifest also has a gzip companion.
The host must serve raw gzip companions unchanged, without HTTP
`Content-Encoding`, and include them in deployment. HTTP compression and caching
for HTML, JavaScript, CSS, SVG and fallback JSON depend on hosting configuration.

The compact map starts with retained trajectory records having adjusted
omnibus p < 0.05. Selecting a point or focusing a gene loads its full records
from `trajectory-details/`; these chunks also contain records above that
threshold. New exports omit the former duplicate `pipeline-trajectory-index.json`
and `trajectory-extras/` layout. The UMAP loader still accepts older releases.
The source contains only trajectories selected for clustering; lookup cannot
supply genes absent from that source. JSON assets are kept in a bounded in-memory cache and fetched with a release-versioned URL
through the browser's HTTP cache, so a repeat visit can reuse downloaded data.

PANTHER metadata and pathway indexes are rebuilt in staging when their source
inputs are available. They are not copied from a previous release. Check the
manifest's availability and source hashes; optional data can be absent on a
fresh checkout. `make site-build` verifies the manifest and every listed asset
before producing `site/dist`. The exporter retains warnings for test-only and
synthetic data; selecting inputs appropriate for the intended audience remains a
release decision.

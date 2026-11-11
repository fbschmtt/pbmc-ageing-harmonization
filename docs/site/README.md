# Static presentation site

[Documentation](../README.md)

The React/Vite app in `site/` reads a static release created from existing
workflow outputs. It does not run analysis or require an API.

The site includes differential-expression, composition, Cell States, age-bin
trajectory, transcriptomic-clock, pathway, dataset, modelling, overview, and
about views. Cell-type and gene selections are shared across views when matching
results exist, and cell-type selectors use the shared type colors. A shared
gene set can be selected on the Pathways page or edited from a plot; it
highlights matching genes in the per-cell-type differential-expression plot and
the age-bin trajectory map. In the per-cell-type volcano plot, gene-set points
have a distinct outline and are drawn above ordinary points; the selected point is drawn last.

The per-gene volcano and pathway heatmap use combined linear age effects per
decade. The per-gene volcano plots the matching combined-model adjusted
p-values on its y axis, colors points by cell type, and uses the same square
plot sizing as the per-cell-type volcano. Exact-gene lookups include
non-significant results; the pathway heatmap does not run an enrichment test.
These views do not require clustered trajectory results. Full trajectory profiles are stored
once for the UMAP; the separate per-type age-bin preview remains available.
Standalone figures are exported as PNG files, deduplicated for single-cell
views and automatically recompressed losslessly when that reduces their size.
HTML reports and their embedded figures are copied unchanged.

Per-cell-type volcano-plot color shows the number of studies with a significant
result for the same contrast, regardless of direction. Zero is gray; missing
support is shown as unavailable. The viridis scale runs from one to the largest tested-
study count among the plotted genes. A switch beside reference-bin expression
colors points by study-balanced expression (log₂ CPM + 1) instead. This mode is
selected by default when study-support data are unavailable. Its color-bar
limits span available genes in the selected cell type's reference-expression
asset. Shared estimate formatting uses two decimals for magnitudes of at least
one and two significant digits below one; plot axes can use their own tick
formatting. Downloaded CSV statistics retain their source precision.

Transcriptomic clocks have cell-type and model selectors. Available models are
listed from baseline to more advanced, the age axes stay fixed across selections,
and prediction points load only for the selected model. The across-cell-types
option shows donor-matched mean predictions when available. New or updated
clock outputs must be exported before they appear in the site.

### Run and refresh the site

The Docker-backed Make targets provide the locked Node environment; a local
Node installation is optional. Export the available pipeline results, then
start the live-reload development server:

```bash
make site-export
make serve-frontend
```

The server opens at `http://localhost:8000/`. Set `SITE_PORT` to change the
port. To serve it to other machines on the network, set `SITE_BIND=0.0.0.0`;
use a static build and a suitable web host for production. The server does not
run analysis or export results.

`make site-export` reads `output/analysis/per_type` by default. Set
`SITE_INPUT` to use another per-type directory, for example:

```bash
make site-export SITE_INPUT=output/test/analysis/per_type
```

Exporting replaces `site/public/data` after staging and validating a complete
release. It does not rerun analysis. Changes to analysis outputs, newly rendered
QC figures, or newly available clocks require another export. `make site-build`
validates the current export and creates the static site in `site/dist/`; run
both targets after changing source results before deploying a static build.
The interactive DE catalogue publishes combined models and the developmental
juvenile-age and age-25-versus-5 results. Per-study DE chunks are removed after
computing study-support counts. The developmental adult-age result is omitted
from the site but remains in the analysis outputs. Combined DE downloads are
prepared from source result tables during export and stored only as `.csv.gz`
archives; the browser extracts them when downloading and saves a plain
CSV with the original source values (per year for linear-age effects). Browsers
without decompression support receive the gzip archive instead.
These storage reductions take effect after re-exporting and rebuilding the site:

```bash
make site-export
make site-build
```

Deploy the complete resulting `site/dist/`, including the PNG assets and gzip
files. Download `.csv.gz` files and JSON/SVG `.gz.bin` companions must be served
as raw compressed bytes without HTTP `Content-Encoding`; the frontend handles
their decompression. Keep the original JSON files for browser fallback support.
The build expects to be served from the host root. See
[the data contract](data-contract.md) for asset contents, provenance,
release validation, and publication considerations.

For test data, run the relevant workflow first and point `SITE_INPUT` at its
per-type output. `make site-test-workflow` exercises existing `output/test`
artifacts without replacing the local export. Browser tests build an isolated
synthetic release and do not change the export. They check downloaded CSV
filenames and exact contents, decompression failures and fallback, matching
volcano sizes on desktop and mobile, model values, and PNG image loading.
Useful frontend checks are:

```bash
make site-test
make site-test-browser
make site-test-workflow
```

Use `make site-test-browser-navigation` for navigation and visualization
journeys or `make site-test-browser-release` for release and loading behavior.
`make verify` includes frontend checks and the workflow-backed site check.

Generated exports, dependencies, and `site/dist/` are ignored by Git. Optional
gene annotations and curated pathway sets are rebuilt during export when their
local source inputs are available; see [INPUT_FILES.md](../inputs/files.md) for
pathway provenance and setup. Review the release manifest and any linked HTML
reports before publication: structured assets omit raw counts and per-cell
identifiers, but reports can contain sample-level details. The report links do
not include their companion downloads.

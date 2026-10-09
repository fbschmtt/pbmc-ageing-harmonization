# PBMC Ageing Harmonization

Configuration-driven processing of PBMC ageing scRNA-seq studies. It produces
one harmonized H5AD per study, per-study sample × AIFI-L2 pseudobulks, and a
cross-study pseudobulk matrix. An explicitly enabled pathway also
creates a raw-count cross-study single-cell merge. A separate, opt-in
integration benchmark creates global embeddings and diagnostic annotations.

Study-specific inputs and adapter dependencies live in `config/studies.json`;
the tracked `config/studies.schema.json` and
`config/input_sources.schema.json` enforce their document structures before a
workflow starts. Explicit adapters under `src/pbmc_pipeline/studies/` own
source-specific joins and reshaping. The output metadata contract lives in
`config/harmonized_obs_schema.json`. See `PLAN.md` for scope, progress, and
unresolved scientific questions. See `IMPLEMENTATION.md` for the high-level
architecture.

See `INPUT_FILES.md` for the complete list of external expression files,
supplementary metadata, CellTypist models, expected paths, and checksums.

## Quick start

The supported workflow environment requires [Nextflow](https://www.nextflow.io/)
version 25.04 or newer and Docker on the host. To fetch the publicly available
inputs, run:

```bash
make download-inputs STUDIES=all
```

> **Large downloads.** Some source files are multi-gigabyte. The downloader
> handles configured public URLs only; sources requiring authentication or a
> browser session must be acquired separately. Workflows skip studies with
> unavailable inputs.

The three [AIFI CellTypist models](aifi_models/README) must be downloaded
separately. Their links require interactive Allen authentication, which the
project downloader does not support. AIDA25, Terekhova23, and Wang25 also have
manual-only inputs; `STUDIES=all` skips studies with incomplete inputs and
records the decision in `output/run_manifest.json`.

After the required AIFI models and the inputs for the studies you want to
include are in their configured paths, run:

```bash
make run-all STUDIES=all
```

This builds the workflow images when needed, runs the core harmonization and
single-cell merge, fits pseudobulk differential-expression models, and renders
all cell-type reports. It uses the Docker profile by default; no local Python
or R installation is required.

## Test quick start

Run this before a production workflow: it confirms that Docker, Nextflow,
models, and the available inputs work before committing to the much longer full
analysis. It creates missing test fixtures for locally available inputs, then
runs the complete ordered Docker test workflow:

```bash
make test-data
make run-all-test
```

By default this requests every configured study and skips unavailable fixtures.
It runs the core pipeline on available downsampled real-data fixtures, then the
integration benchmark and report checker, then a generated pseudobulk fixture
for positive DE fits and reports. The benchmark artifact retains labels and UMAP
coordinates without raw counts or neighbour graphs. Use `STUDIES=<list>` to
narrow the core run. Existing fixtures are preserved; use
`TEST_DATA_OVERWRITE=true` only to refresh them.

## Setup details

Input acquisition is a separate, strictly opt-in step: `make install`, tests,
and workflow runs never download data. Use [Download inputs](#download-inputs)
for acquisition commands and `INPUT_FILES.md` for paths and checksums. The
AIFI models require interactive Allen authentication. Other manual sources are
study-specific; an all-studies workflow skips a study when any required input
is unavailable.

For production, `make run STUDIES=all` analyses every study whose expression
input and declared preparation dependencies are present. It skips unavailable
studies with a warning and records requested, selected, and skipped studies in
`output/run_manifest.json`. An explicit list such as
`make run STUDIES=aifi,terekhova23` is strict: missing inputs are an error, so
the selected analysis population is never reduced silently. Test workflows
instead skip unavailable fixtures and fail only when no requested fixture is
available.

For local Python development or host-side tools such as linting and unit tests,
create the local environment with `make install`.
It is not a substitute for Docker for the complete workflow: in particular, it
does not provision the R, Seurat, reticulate, and anndata stack used for RDS
conversion.

## Current status

Current implementation work, validation status, and remaining scientific
review are tracked in [PLAN.md](PLAN.md). Run the focused target in
[Choosing a Make target](#choosing-a-make-target) before relying on a changed
workflow or report.

## Data flow

```text
input_data/ expression objects + metadata
                │
                ├── Seurat RDS ──> converted H5AD
                │
                ▼
 prepare canonical metadata plus a retained source-obs audit sidecar
                │
                ├──> output/prepared/<study>.cells.csv.gz
                └──> output/prepared/<study>.source_obs.csv.gz
                ▼
      harmonize one configured study
                │
                ├──> output/harmonized/<study>.h5ad
                ├──> output/reports/<study>.json
                └──> output/qc/<study>/{executed.ipynb,report.html}
                │
                ▼
  pseudobulk each study by sample/AIFI-L2
                │
                └──> output/pseudobulk/<study>.pseudobulk.h5ad
                ▼
       output/merged/pseudobulk_merged.h5ad
                │
                ├── optional: output/merged/single_cell_merged.h5ad
                └──> output/qc/merged/{executed.ipynb,report.html}

optional global integration benchmark
  single_cell_merged.h5ad ──> output/integration_benchmark/
                                └──> labels-and-UMAP-only H5AD + JSON + report.html

optional downstream cell-type workflow
  single_cell_merged.h5ad ──> split once by per-study AIFI-L2
                                ├──> output/cell_type_splits/<type>.h5ad
                                └──> analyse each type
                                      └──> analysis.h5ad + tables
                                            └──> render notebook/HTML
                                                  └──> output/cell_type_analysis/<type>/
                                                        ├── report.html
                                                        └── executed.ipynb + tables, including
                                                            fraction_model_diagnostics.tsv
                              output/cell_type_analysis/cell_type_manifest.json

optional pseudobulk differential-expression workflow
  pseudobulk_merged.h5ad ──> list per-study AIFI-L2 labels
                              └──> one parallel task per label, each loading and
                                   subsetting the same small pseudobulk H5AD in memory
                                     ├──> per-study fits with site adjustment
                                     └──> merged ~ study_site + age + sex + log10_total_counts fit
                                           └──> collector: gene-level CSVs + indexed JSON manifest
```

All immutable local inputs belong under `input_data/`. Generated conversions,
reports, and harmonized outputs are kept outside it and ignored by Git.

Per-cell preparation artifacts are generated under `<outdir>/prepared/`; they
are reproducible from the expression object plus the declared study dependencies.
Each study also publishes `<study>.source_obs.csv.gz`, the retained source
observation metadata keyed by `cell_id`. This sidecar is an audit trail for
adapter mappings; the harmonized H5AD continues to expose only the canonical
metadata schema.

### Annotation semantics

`aifi_l2_majority` is the retained per-study AIFI-L2 annotation and the sole
canonical label for pseudobulk grouping, cell-type splitting, and composition
denominators. The optional integration benchmark produces diagnostic
`benchmark_*_aifi_l2_majority` labels; those must not replace per-study L2 in
any downstream analysis.

## Add a study

Start with [templates/study_adapter.py](templates/study_adapter.py) and
[templates/study_config.json](templates/study_config.json). Copy the adapter to
`src/pbmc_pipeline/studies/<study_id>.py`, replace its rules with the source's
cell and sample metadata, then add the configuration object under `studies` in
`config/studies.json`. Add every expression object and supplementary metadata
file to `config/input_sources.json` and `INPUT_FILES.md`; paths must be relative
to the repository root. Both registry documents must satisfy their tracked JSON
Schemas as well as the pipeline's cross-reference checks. Declare supplementary files in
`preparation.dependencies` so Nextflow stages them and configuration changes are
tracked in task identity.

An adapter receives the expression object's `.obs` table and must return the
canonical metadata fields, with exactly one row for each retained cell. Use
`safe_left_join` for sample-level supplements: it rejects ambiguous joins and
preserves the expression-cell order. The usual route is to retain every input
cell. When an object deliberately mixes cohorts, filter it explicitly in the
adapter, explain the condition in `provenance.selection`, and set
`preparation.allow_cell_subset` to `true`. The preparation and harmonization
steps then verify that retained IDs are a strict subset of the source IDs and
drop the excluded cells from the count matrix before annotation. This prevents
metadata-only filtering from leaking excluded cells into downstream outputs.

Choose `counts_source` as `X`, `raw`, or `layer:<layer_name>` after checking
which matrix contains unnormalized integer counts. The harmonized output always
stores only those selected counts in `.X`; source `.raw` and every source layer
are removed.

When an H5AD uses Ensembl IDs as `.var_names` but carries gene symbols in a
feature metadata column, configure `features.source_column` to that symbol
column (usually `feature_name`). The harmonizer repairs `.var_names` from this
column and sums duplicate symbols according to `duplicate_policy`; it does not
use the Ensembl index for downstream gene joins. Fachrul26 and Perez22 follow
this pattern. Harmonization reports record the source column, symbol-like
feature count, identifier-like feature count, missing names, and duplicate-name
counts so source annotation quality is visible rather than implicit.

Configuration loading also checks adapter module existence, count-source and
feature-repair policies, conversion declarations, and that every study input
and preparation dependency is represented in `config/input_sources.json`.

Set `technology` to the most specific reported assay information. For example,
use `10X3'v3` rather than `10X3'` when the reagent-kit revision is available;
use the assay-family value only when no more specific source information is
reported. Record the source field in `provenance`, and retain a per-cell run or
library identifier in `batch_single_cell` when the source supplies one. Record
every inferred metadata value in `provenance` with its rationale; do not make
it look like a reported source field.

Before enabling the study in an all-study run, verify that the adapter's cell
selection, sample-level joins, count source, technology/provenance fields, and
input download declarations match the publication and supplements. Generate a
test fixture with
`make test-data`, then run `make run-test STUDIES=<study_id>` and inspect its
prepared-table report, harmonized H5AD, and QC report. Confirm that the output
has integer counts in `.X`, `raw is None`, and no layers before adding it to an
all-study run.

## Sample and metadata semantics

`subject` is a biological individual. `sample` is one biological specimen at
one collection timepoint; it is the unit aggregated into sample × AIFI-L2
pseudobulks. `batch_single_cell` is a technical library, run, pool, or well and
never defines a pseudobulk. Missing string metadata is written as
`not_provided`; missing numeric metadata is `NaN`.

| Study | Reference DOI | Subject | Sample | Sampling timepoint | Technical batch | Notes |
|---|---|---|---|---|---|---|
| AIDA25 | [10.1016/j.cell.2025.02.017](https://doi.org/10.1016/j.cell.2025.02.017) | `donor_id` | country-qualified donor ID | not provided | supplementary experimental batch | Source `obs.assay` reports 10x 5′ v2. Lonza material replicated across sites remains intentionally split until downstream exclusion. BMI uses the supplementary donor metadata. |
| AIFI | [10.1038/s41586-025-09686-5](https://doi.org/10.1038/s41586-025-09686-5) | specimen GUID prefix | specimen GUID | `sample.visitName` | `well_id` | A specimen GUID may span multiple wells; pseudobulk combines those technical partitions. The source H5AD has no assay field. |
| OneK1K | [10.1126/science.abf3041](https://doi.org/10.1126/science.abf3041) | `donor_id` | `donor_id` | not provided | `pool_number` | No separate specimen or timepoint key is present in the source H5AD; the sample-per-donor assumption remains under review. |
| Terekhova23 | [10.1016/j.immuni.2023.10.013](https://doi.org/10.1016/j.immuni.2023.10.013) | `Donor_id` | `Tube_id` | workbook visit | barcode-derived batch | Each tube maps to one donor and visit; a tube may span technical batches. The paper explicitly reports frozen input material. |
| Wang25 | [10.1038/s41590-024-02059-6](https://doi.org/10.1038/s41590-024-02059-6) | `Sample ID` | `Sample ID` | not provided | not provided | No separate donor, timepoint, or technical-batch identifier has yet been recovered. Samples are provisionally treated as fresh based on the public reference-sample selection; author confirmation is pending. |
| Nehar-Belaid26 | [10.1038/s41467-026-73729-2](https://doi.org/10.1038/s41467-026-73729-2) | Supplementary Data 1a `IDs` | Supplementary Data 1a `Names` | not provided | Supplementary Data 1a `runs_10x` | Reads the supplied GEO raw 10x tar directly (without extraction); `JB` subject IDs parsed from member names join to Supplementary Data 1a `IDs`, which supplies sample IDs and donor metadata. Published labels transfer only on a unique sample ID plus normalized barcode-sequence match. |
| Fachrul26 | [10.64898/2026.02.15.704933](https://doi.org/10.64898/2026.02.15.704933) | `donor_id` | `sample_id` | not provided | `library_id` | Paper methods state PBMCs were cryopreserved and thawed; this takes precedence over the conflicting CELLxGENE `obs.sample_preservation_method=fresh`. |
| Perez22 | [10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970) | `donor_id` | `sample_uuid` | not provided | `library_uuid` | Paper reports frozen input material. Uses embedded H5AD metadata; `feature_name` supplies gene symbols. `obs.disease` is normalized and only `normal`/healthy cells are retained; embedded SLE cells are excluded before harmonization. |

`tests/` contains only automated test code. Reproducible local intermediates
are collected under `cache/` for RDS-to-H5AD conversions. Downsampled smoke-test
fixtures are generated under the ignored `test_data/` directory.

## Download inputs

> **Strictly opt-in and potentially large.** Only `make download-inputs` and
> `make download-inputs-strict` retrieve data; `make install`, test targets,
> and production workflow targets never do. Before it starts, the downloader
> warns about pending artifacts and their recorded size where available. A
> selected study can require multi-gigabyte source files.

Specify the studies to acquire. The best-effort command downloads configured
public URLs without portal authentication or scraping and records outcomes in
ignored `input_data/download_manifest.json`. It returns a failure if a direct
download fails or fails its size/checksum check; manual-only inputs remain
reported for best-effort acquisition:

```bash
make download-inputs STUDIES=aifi
make download-inputs STUDIES=nehar_belaid26
```

The downloader prints each artifact's position in the run and its associated
study name. During a transfer it also shows bytes downloaded, transfer rate,
and percentage when the source reports a size or the manifest records one;
otherwise it reports downloaded bytes and rate.
Interrupted direct downloads retain a compatible `.part` file and transfer
metadata. The downloader makes up to three attempts with delay and resume
messages, and resumes with HTTP byte ranges only when the source confirms that
the saved partial is compatible. A final size and checksum check is still
required before publishing the requested file.

`make download-inputs-strict STUDIES=fachrul26` is the explicit acquisition
preflight. In addition to failed downloads and size/checksum mismatches, it
exits nonzero for a manual-only artifact that has not yet been placed at its
configured path. It is an acquisition preflight, not a routine test, and never
runs as part of analysis commands. Run it before production when download
integrity matters. Existing files are left unchanged unless
`DOWNLOAD_FORCE=true` is supplied deliberately.

The manifest records manual-only sources and their acquisition instructions.
Synapse expression objects require portal access, the AIDA publisher
attachment rejects unauthenticated programmatic requests, and the AIFI full
H5AD link currently redirects through Allen's Google-authenticated workspace.
The downloader does not implement these access flows. Some files currently
present at production paths are test-scale subsets. Recorded SHA-256 values and
CELLxGENE byte sizes are verified; confirm that acquired files are suitable for
production analysis.
See `config/input_sources.json` and [INPUT_FILES.md](INPUT_FILES.md).

## QC reports

Passing `--qc` executes the tracked `reports/qc_report.py` template after each
harmonized H5AD is written. Each study gets a separate directory:

```text
output/qc/<study>/
├── executed.ipynb
└── report.html
```

The report includes provenance, metadata completeness, cohort
composition, count-depth checks, UMAPs, canonical PBMC marker expression, and
study-label/AIFI-label comparisons. All plots use the complete harmonized study.

To generate a report for an existing harmonized output without rerunning the
study:

```bash
pbmc-qc \
  --input output/harmonized/onek1k.test.h5ad \
  --run-report output/reports/onek1k.test.json \
  --output-dir output/qc/onek1k.test \
  --study onek1k
```

The merge workflow renders one combined document at `output/qc/merged/`. It
covers pseudobulk composition and technical covariates; the pseudobulk outer
union retains all input counts. Shared-gene availability, hypothetical count
loss under the all-study intersection, and the highest-count excluded gene are
reported with the cross-cell-type expression atlas, where that intersection is
actually used. Global UMAPs and benchmark AIFI-L2 comparisons belong to the
separate integration benchmark.
The merged artifact retains per-study labels in `aifi_l2_majority` for all
downstream grouping and splitting. Run `make run-integration-benchmark` after
the core workflow to create a small diagnostic H5AD with global UMAP coordinates
and CellTypist labels from the Harmony and unintegrated-PCA graphs. It does not
duplicate the merged raw count matrix. Its `report.html` contains the global
UMAPs and per-study-versus-benchmark AIFI-L2 concordance matrices.
`output/run_manifest.json` records requested, selected, and skipped studies,
image references, configuration checksum, and checksums of the merged
deliverables.

## Nextflow and Docker

Nextflow runs the configured studies in parallel and caches completed tasks.
The Docker profile is the supported workflow environment and is selected by
default. The test profile uses isolated 200-cell fixtures from `test_data/`;
RDS fixtures still use the normal conversion process. `make test-data` creates
only missing fixtures and never changes `input_data/`; use
`TEST_DATA_OVERWRITE=true` to refresh them. Pass `STUDIES=<list>` to narrow a
test or production run.

## Choosing a Make target

Choose the smallest target that directly exercises the changed contract. The
repository instructions in [AGENTS.md](AGENTS.md) use this table when deciding
which validation to run.

| Change or goal | Preferred target | When to use a broader target |
|---|---|---|
| Python library, CLI, or deterministic test change | `make lint test-unit` | Run a focused workflow too if a task contract or published artifact changed. |
| Nextflow workflow or configuration change | `make workflow-lint` | Add `make run-test STUDIES=<study>` when task wiring must execute. |
| One cell-type report or analysis | `make run-cell-type-test CELL_TYPE=<slug>` (requires an existing split; create test splits with `make split-cell-types-test`) | Use `make render-cell-type-test CELL_TYPE=<slug>` to rerender one report from existing analysis artifacts, or `make render-cell-type-reports-test` to rerender all test reports. |
| Rerender reports from existing analysis artifacts | `make render-cell-type-reports[-test]` | Use `make render-cell-type[-test] CELL_TYPE=<slug>` to refresh only one report. |
| Cell-type splitting or shared analysis | `make run-test`, then `make run-cell-type-analysis-test` | Use `make run-all-test` when the complete ordered test workflow is needed. |
| Positive pseudobulk DE fitting | `make run-de-synthetic-test` and `make check-synthetic-de-results` | This fixture is independent of the core merge. To check report integration, also run `make run-test` followed by `make run-cell-type-analysis-test` and `make check-synthetic-de-reports`. |
| DE from an existing test pseudobulk merge | `make run-de-test` | Requires `make run-test` first; the small core fixture may not have enough independent samples for successful fits. Use the synthetic target above for positive fitting coverage. |
| Cross-cell-type trajectory report only | `make render-trajectory-report[-test]` | Reuses existing DE results and reruns only the combined trajectory analysis and notebook. |
| End-to-end trajectory path | `make run-trajectory-test` | Reuses an existing test merge, then runs the config-sized synthetic DE fit, both report levels, and artifact checks; skips the integration benchmark. Run `make run-test` first if the merge is missing. |
| Core workflow for one or a few studies | `make run-test STUDIES=<list>` | Omit `STUDIES` to request all configured studies. |
| Global integration/annotation diagnostics | `make run-integration-benchmark[-test]` | Run after an existing single-cell merge; `make run-all-test` also exercises its test form. |
| Ordered Docker integration test | `make run-all-test` | Use `make verify` when local lint, unit, and documentation checks are also needed. |
| Markdown documentation or Make target names | `make docs-check` | Also run the target described by changed commands. |
| Broad cross-cutting or release validation | `make verify` | This is deliberately non-resumed and includes all studies by default. |
| Production analysis | `make run` | It is resumable; use only after the relevant focused checks pass. |

Nextflow submission alone is not a passing workflow check: confirm each task's
`.exitcode` and the published deliverables. Development test targets request
all configured studies and skip unavailable test fixtures; pass an explicit
comma-separated `STUDIES` list to narrow a run.

Harmony integration uses the pinned Harmony2 C++ backend. It uses the backend's
native thread policy; the legacy OpenMP/OpenBLAS thread overrides that crashed
on large production runs have been removed.

Cell-type reports provide local embeddings, marker contrasts, study coverage,
descriptive fraction trends, and optional DE summaries. Their model and
interpretation details are documented in the downstream section below.

`make run-test` always creates the test single-cell merge. Production runs
create it only through `make run-all` or when explicitly enabled:

```bash
make run STUDIES=all MERGE_SINGLE_CELL=true
```

To benchmark the global Harmony and unintegrated-PCA label/UMAP variants without
rerunning that merge:

```bash
make run-integration-benchmark
```

Study-QC, merge-QC, integration-benchmark, and cell-type report templates are
Jupytext Python sources under `reports/`. Their runners materialize a temporary
notebook immediately before execution; generated templates are not tracked.
The published `executed.ipynb` is the audit artifact.

The merge QC report summarizes sample-level technology, intronic-read
inclusion, and frozen-status distributions with provenance status. OneK1K's
sample counts use `donor_id` until its biological collection identity is
confirmed.

### Report artifacts and web delivery

Each report run publishes a self-contained static `report.html`, which can be
opened directly in a browser or served from a static file host, plus an
`executed.ipynb` that retains code and outputs for technical reproduction. The
HTML report hides implementation cells and execution prompts, adds a compact
summary/navigation layer, and is the reader-facing analysis artifact; the
notebook is the audit artifact. These standalone pages are intentionally
separate from any future project-level website or custom frontend, which should
browse manifests and structured report tables rather than parse notebook HTML.

#### Automatic table of contents

Cell-type and merge QC reports share their reader-facing header and styles from
[`src/pbmc_pipeline/report_theme.py`](src/pbmc_pipeline/report_theme.py). Their
navigation lists are generated from the final rendered HTML after notebook
export. The templates provide a TOC placeholder, and
[`src/pbmc_pipeline/report_toc.py`](src/pbmc_pipeline/report_toc.py) collects
the rendered H2, H3, and H4 headings in document order, nests them by heading
level, and replaces the placeholder with linked entries. This also picks up
headings emitted conditionally or programmatically by notebook cells. When
editing a report, add or reorder document headings and let the renderer update
the TOC; do not maintain a separate list of TOC links. Give headings explicit
IDs when their anchors need to remain stable across heading text changes.

The cell-type and merge QC runners call these helpers after notebook execution
and HTML export. Other report templates do not use them unless their runner is
updated to add the shared header and TOC placeholder, then populate the TOC.

## Per-cell-type reports

The optional downstream workflow loads a completed merged single-cell H5AD once
and splits it by the retained per-study `aifi_l2_majority` call. Experimental
benchmark labels are never used for this split. For each type, `ANALYSE_CELL_TYPE`
creates the type-specific PCA/UMAP/Leiden embedding, marker tables, PC--age
correlation table, and `analysis.h5ad`. The PCA section shows each principal
component's individual share of variance explained, not cumulative variance.
`RENDER_CELL_TYPE_REPORT` then loads
that analysis artifact plus the primitive raw-count H5AD to create the executed
notebook and HTML. It also plots the retained per-study `aifi_l3_majority`
labels; no merged L3 prediction is created. Global neighbours are deliberately
not copied into the split artifacts.

The main adjusted-fraction forest plot fits adult-only OLS models separately by
study. Its age coefficient is per decade; sex, BMI, and CMV terms are shown when
estimable. Fractions are modeled on the 0–1 scale and effects, confidence
intervals, and residual SDs are displayed as absolute percentage-point changes
(0.01 fraction = 1 percentage point). Confidence intervals extending beyond
the displayed scale are marked with arrowheads. The title identifies the
adult-only analysis, current L2 cell type, and its L1 parent.

Secondary descriptive views plot the type's sample-level fraction across age,
coloured by study. They show two denominators, both recorded before splitting:

- all retained PBMCs in the same `study × sample` (`n_cells_in_sample`);
- all retained cells in the configured AIFI-L1 parent compartment in that
  `study × sample` (`n_cells_in_sample_l1_parent`).

For the latter, every retained cell is first mapped from its retained per-study AIFI-L2
label to the configured theoretical L1 parent, then counted by `study × sample
× derived parent`. The workflow intentionally does not use `aifi_l1_majority`:
that is an independent classifier result and may disagree with the refreshed
merged L2 annotation.

The L2→L1 parent mapping is manually curated in `cell_type_analysis.l2_parent_l1`
in `config/pipeline.json`. Its provenance records that it is inferred from AIFI
label semantics, checked against the shipped model vocabularies, and must be
validated against the AIFI atlas/reference taxonomy before inferential use.
Both plots are descriptive composition views, not age-effect tests. The split
code retains a comment for future per-sample compartment counts, allowing more
hierarchy views without extra source files.

Fraction plots exclude samples with fewer than 10 cells in the relevant
denominator (`cell_type_analysis.min_fraction_denominator_cells`): retained
PBMCs for the first view, the derived L1 parent for the second, and the split
cell type for local-cluster fractions. This guard prevents unstable fractions
from tiny sample compartments. The study facets show separate within-study
linear fits using all samples and using samples aged at least 20 years; points
are colored by study, with solid black and dashed gray fit lines. Their fit
legend appears once.

The all-study view puts all studies in each denominator panel: sample points
have alpha 0.2 and each study has its own descriptive lowess curve. It does not
fit a combined model or imply a pooled age effect. The within-study trend panels
are open on load; the all-study lowess figure remains collapsible after the
forest plot.

```bash
make run-cell-type-analysis
```

This downstream workflow requires an existing merged single-cell H5AD at
`<outdir>/merged/single_cell_merged.h5ad`. Override that location with
`MERGED_INPUT=/path/to/single_cell_merged.h5ad`.

For a complete production run, `make run-all` executes the core workflow with
single-cell merging enabled, then runs differential expression and cell-type
analysis/reporting in order. The workflow is resumable. It uses the configured
`STUDIES`, `OUTDIR`, and `WORK_DIR` values. For test output, run `make run-test`
first, then `make run-cell-type-analysis-test` to consume the test merge under
`output/test/`. All test targets publish only beneath `output/test/`.

The raw-count splits are published under `<outdir>/cell_type_splits/` so one
report can be rerun without re-splitting the merged H5AD:

```bash
make run-cell-type CELL_TYPE=cd14-monocyte
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type-reports
make render-cell-type-reports-test
```

To deliberately discard published production outputs while preserving the
ignored `output/test` fixture outputs, use one of:

```bash
make clean-core-output      # all production outputs, including stale downstream analysis
make clean-analysis-output  # differential expression and cell-type outputs
```

These are separate because re-running DE or report analyses normally reuses a
compute-heavy production merge. A core cleanup also removes downstream outputs
that no longer match a future merge. Neither target removes Nextflow's cache; use
`make clean-work` only when that cache itself should be discarded.

The targeted test command deliberately does not build a core test merge or all
reports; it requires the corresponding split to already exist. It runs the
analysis and rendering stages for that one type. Create only the reusable splits
from an existing merge with `make split-cell-types-test`, or run `make run-test`
followed by `make run-cell-type-analysis-test` for all types. A targeted rerun
refreshes that type directory but not `cell_type_manifest.json`; rerun the
all-type workflow before treating the published report set as a new
manifest-backed release. `cd14-monocyte` is the default targeted test type
because it is well represented in the bundled multi-study fixture.

For a report-only edit, `make render-cell-type[-test]` rerenders one type, while
`make render-cell-type-reports[-test]` rerenders every type. These targets
require the published splits and `analysis.h5ad` artifacts. They schedule only
rendering and overwrite the executed notebooks and HTML without rerunning local
PCA, Harmony, UMAP, clustering, or marker ranking.

Use `make check-cell-type-test-prerequisites CELL_TYPE=<slug>` before a
targeted test rerun when unsure whether its input exists. It prints the minimal
next command: create splits from the existing test merge, or create a fresh core
test merge and reports when that merge is absent. The production equivalent is
`make check-cell-type-prerequisites CELL_TYPE=<slug>`.

Each published type directory under `output/cell_type_analysis/` contains its
derived `analysis.h5ad`, HTML report, executed notebook, and report tables.
Reusable primitive raw-count splits are published separately under
`output/cell_type_splits/`; they are inputs for targeted report reruns, not
standalone final analyses. `cell_type_manifest.json` records the exact merged input, configuration,
reference-model checksums, analysis specification, and the status/checksums of
every expected type report. Its creation fails if a report directory is missing
its required deliverables, so a successful workflow cannot silently publish a
partial report set.
PC--age correlations in these reports are explicitly descriptive at the cell
level; use sample/subject-level pseudobulks or mixed models for inference.
The local neighbour graph, UMAP, and Leiden clustering use Harmony correction
by study; native (uncorrected) PCs remain the basis of the current PC--age
diagnostic. This separation is deliberate but should be revisited before any
inferential downstream use, since integration can also alter age-associated
structure. V(D)J genes matching `processing.exclude_vdj_regex` are excluded
from the local HVG/PCA input; the report records the count excluded.

Local embeddings are exploratory. The current eligibility threshold is 50 total
cells and does not require balanced study representation; small or highly
imbalanced types can yield visually structured UMAPs even after Harmony. In
particular, the configured 50-neighbour graph is relatively dense for rare
types. Inspect study-coloured UMAPs and per-study cluster composition before
interpreting local clusters as biological states.

Nextflow's `work/` directory is the cache for task intermediates. Published
outputs are copied to `output/`, including final reports, reusable cell-type
split H5ADs, and the derived per-type analysis H5ADs; do not use paths inside
`work/` as report inputs.

Workflow targets build their required content-tagged images automatically; use
`make images` only to build them without a run. For a production run whose
repository and task filesystem are separate, keep the repository readable by
Docker and choose external work/output locations:

```bash
make run NF_PROFILE=docker STUDIES=all \
  WORK_DIR=/path/to/run/work OUTDIR=/path/to/run/output
```

Nextflow stages declared inputs into the task directory. Keep the repository
readable to Nextflow and use a `WORK_DIR` that supports symbolic links.

## Pseudobulk differential expression

### Gene universes

The workflow intentionally uses different gene universes for different tasks:

| Context | Gene universe | Purpose |
|---|---|---|
| Merged pseudobulk | Outer union across studies | Retains each study's observed genes for per-study models. |
| Per-study DE | Genes available in that study | Fits only observed counts for that study. |
| Combined DE | Intersection across studies included in that fit | Prevents study-absent outer-join zeros from entering the model. |
| Study-balanced expression atlas | Intersection across all studies | Makes its equally study-weighted expression values comparable. |

The separate `differential_expression.nf` workflow consumes the merged
sample-level pseudobulk H5AD. It first lists retained per-study
`aifi_l2_majority` labels, then runs one independent task per label. Each task
loads the same small merged pseudobulk H5AD and subsets its label in memory; it
does not materialize a per-type pseudobulk H5AD. Per-study models adjust for
`study_site` when multiple sites remain in that study, along with age, sex, and
available optional covariates. The merged age-and-sex model uses
`~ study_site + age + sex + log10_total_counts` across studies. `log10(total_counts)`
from each pseudobulk is included as a technical covariate in all DE models and
its coefficient is shown in the combined-model results and volcano plot. When
present and variable, BMI and CMV each receive separate combined fits adjusted
for study site, age, sex, and log10 total counts.
DESeq2 estimates size factors to normalize library depth, and this additional
depth term is an empirical adjustment motivated by the observed data: pseudobulk
depth varies with age, especially in Naive CD8 T cells where cell-type abundance
declines strongly with age, and sample-level residuals retained depth-dependent
patterns after fitting. Including `log10(total_counts)` explicitly adjusts for
that age-associated technical variation beyond size-factor normalization. It
reduces this confounding risk; it does not establish that depth explains the
original volcano pattern or remove the need to inspect residual diagnostics.
Combined fits include studies with usable values for that covariate and report
the studies, sites, and sample counts used. The site term is dropped when only
one site remains. Age, sex, BMI, CMV, and log10(total_counts) effects each have
a combined-model volcano plot where estimable. Age and log10(total_counts) also
have volcano plots colored by log10 of the combined-fit baseMean. A final collector writes the result
directories and manifest. Per-study results also show significant-gene counts
inside and outside the per-covariate intersection of genes tested by all
available study fits. The gene-universe table above is the canonical definition
of which genes enter each fit. The parallel type-fitting task is configured for
one CPU, which is passed directly to PyDESeq2; all other CPU, time, and memory
requests use the executor defaults.
In addition to the continuous-age fit, each cell type gets an exploratory,
site- and covariate-adjusted age-trajectory analysis. Its report shows
significant trajectory shapes and clusters; a companion cross-cell-type report
summarizes shared trajectory and Pearson-residual structure. That report also
contains a descriptive, study-balanced expression atlas, including shared-gene
availability and count-loss diagnostics, expression maps, and technical-label
contrasts. These views complement rather than replace the per-cell-type DE
results, and they are not covariate-adjusted estimates of absolute RNA
abundance. Per-type artifacts are stored under each type's `combined/`
directory and cross-type artifacts under
`differential_expression/trajectory_analysis/`. See
[IMPLEMENTATION.md](IMPLEMENTATION.md#age-trajectory-settings) for the model,
support rules, clustering, output tables, and current configurable settings.

Each cell type also renders
one shared-age diagnostic figure: study-site age support, age against raw
pseudobulk depth in the tested gene intersection,
and a Q-Q plot of unadjusted p-values from the combined age model fitted across
all eligible studies for that cell type and adjusted for study site, sex, age,
and log10 total counts. Its dashed y=x line marks the null
reference, with independently scaled axes for readability. Combined-model
standard volcano point colors show the number of available per-study fits for the same
covariate and contrast where that gene passes the FDR threshold; point position
continues to show the combined-fit effect and adjusted p-value. On the age
volcano, outlined study-coloured diamonds additionally show each study's most
significant gene that passes the configured FDR threshold and lies outside the
combined all-study gene intersection, using that study's own effect and
adjusted p-value. Therefore a
diamond cannot duplicate a combined-fit point; it is a context marker for a
gene omitted from the combined model.

Samples must be age 20 or older and have at least 10 cells in that
sample × cell-type pseudobulk. Samples missing required metadata or covariates,
or with zero counts over the fit's gene universe, are also omitted. Optional
covariate missingness receives additional handling within each study: a
per-study model includes an optional covariate only when it varies in that
study, and samples missing a covariate are omitted from fits that include it.
Combined optional-covariate models use complete cases from studies that recorded
the covariate. All sample exclusions and model fits that cannot estimate the
configured covariates are recorded in the JSON manifest. CSV results and the manifest are published
beneath `<outdir>/differential_expression/`; each cell-type directory also
contains the task fit record (`cell_type_result.json`) used by the collector.

The DE manifest is the handoff contract for reports. Its
`results_by_cell_type` index lists the per-study and combined result paths for
each covariate (or an empty result set) for each cell-type slug. The cell-type
workflow reads this manifest and passes only the matching cell-type results to
a single report renderer. `make render-cell-type-reports` rerenders all existing
reports from the published splits and analysis artifacts without repeating
analysis. The renderer writes all single-cell artifacts beneath
`<outdir>/cell_type_analysis/`; it never publishes report files into the DE
directory. Named Nextflow subworkflows group splitting plus analysis, analysis
of an existing split, and report rendering while retaining the separate
single-cell and pseudobulk entry points.

The Nextflow `test` profile used by `make run-de-test` enables a test-only mode
that bypasses only the minimum-cell cutoff. It keeps the adult-age, metadata, and model-estimability
checks, so the small fixture can exercise PyDESeq2 without changing the
production inclusion rules. Every test-mode run is labeled `test_only`; its
root and per-celltype metadata, result CSVs, and rendered reports carry a clear
warning that the results are for software validation and must not be treated
as biological evidence. Use `make run-de` for the configured analysis rules.

The small core fixture may produce no DE fits because bypassing a cell cutoff
does not add independent samples or fix a rank-deficient design.
`make run-de-synthetic-test` instead generates a deterministic, test-only
pseudobulk fixture with balanced age support, planted covariate and
study-specific age signals, gene-availability gaps, and technical-label
contrasts. It validates fit and report wiring—not biological power or expected
real-data results. The fixture design and marker contract are documented in
[IMPLEMENTATION.md](IMPLEMENTATION.md#running-and-verification). Do not
upsample real test samples to simulate replication; more cells within a sample
do not add independent samples.

For the focused end-to-end trajectory path, use `make run-trajectory-test`
after `make run-test`; see [Choosing a Make target](#choosing-a-make-target)
for the smallest target for other test scenarios.

`make run-all-test` is the complete ordered Docker test. It runs the core and
integration benchmark on available downsampled real-data fixtures, then fits
the separate synthetic pseudobulk DE fixture and renders reports from the real
merge/splits with those DE results. `make verify` runs it after linting, unit
tests, Nextflow lint, and documentation checks. `STUDIES=<list>` narrows only
the core test inputs.

The generated input stays under `output/test/synthetic_de/`. DE outputs use the
standard `output/test/differential_expression/` directory, and the ordinary
cell-type report target renders them under
`output/test/cell_type_analysis/`.

Cell-type report generation checks for matching DE result files under that
directory. When available, it appends per-study significant-gene recurrence
within the genes tested in every per-study fit, total significant-gene counts,
and combined-model volcano plots for the available covariates at the bottom of
the report. Each plot is
omitted if its corresponding result files are absent.

Run the workflow after creating the pseudobulk merge:

```bash
make run-de
```

When only the cross-cell-type trajectory analysis or its notebook changes,
rerender it from the existing DE manifest and result directories without
refitting the models:

```bash
make render-trajectory-report
make render-trajectory-report-test
```

The test form requires `output/test/differential_expression/` from
`make run-de-synthetic-test`; it does not require the core test merge.

By default this reads `output/merged/pseudobulk_merged.h5ad`. Override the
input or output roots with `PSEUDOBULK_INPUT=/path/to/pseudobulk_merged.h5ad`
or `OUTDIR=/path/to/output`. Override where cell-type report generation looks
for DE results with `DE_RESULTS_DIR=/path/to/differential_expression`.
`make run-de-test` requires an existing test pseudobulk merge; create it first
with `make run-test`. It uses the explicitly labeled test mode described above.
Use `make run-de-synthetic-test` to exercise successful PyDESeq2 fits with the
production sample-inclusion filters enabled.

Each study has one maximal DE model: `study_site` when it varies, age, sex,
log10 total counts, and every optional covariate with at least two observed
values. All coefficients for that study therefore come from one formula and
complete-case sample set.
A missing BMI or CMV value
excludes that pseudobulk only from models containing that term; one missing
value removes one row. Combined models target one covariate at a time, using
all complete-case samples from studies that recorded it and adjusting for
study site, age, sex, and log10 total counts. The site term is omitted when one
site remains. CMV
negative/positive values are harmonized
to no/yes. Categorical volcano plots and summaries name their comparison
groups, including male versus female and yes CMV versus no CMV.

## Further reading

See [INPUT_FILES.md](INPUT_FILES.md) for exact external paths and checksums,
[IMPLEMENTATION.md](IMPLEMENTATION.md) for architecture and output contracts,
and [PLAN.md](PLAN.md) for open scientific and reproducibility work.

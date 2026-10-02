# PBMC Ageing Harmonization

Configuration-driven processing of PBMC ageing scRNA-seq studies. It produces
one harmonized H5AD per study, per-study sample × AIFI-L2 pseudobulks, and a
cross-study pseudobulk matrix. An explicitly enabled pathway also
creates a jointly embedded, freshly AIFI-L2-annotated single-cell merge.

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

By default this requests every configured study and
skips any unavailable fixture. It runs the core pipeline on available downsampled
real-data fixtures, then uses a separate generated pseudobulk fixture for
positive DE fits, and finally renders reports from the core test merge with
those synthetic DE results. Use an explicit `STUDIES=<list>` to narrow the core
run.
Existing fixtures are preserved; use `TEST_DATA_OVERWRITE=true` only when they
need refreshing.

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

The last recorded `make verify` passed with all eight configured test studies.
It validated preparation, RDS conversion, harmonization, QC, pseudobulk and
single-cell merges, positive synthetic DE fits, and 21 cell-type reports. The
full production workflow has completed with the real Wang25 RDS, and a
fresh-clone run completed input acquisition and `make run-all`. Those runs
predate the latest metadata corrections; rerun production after the final
manual metadata check. Remaining provenance and follow-up work is listed in
[PLAN.md](PLAN.md).

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
                                     ├──> per-study ~ age + sex fits
                                     └──> merged ~ study + age + sex fit
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
denominators. The optional whole-dataset merge produces only diagnostic
`experimental_aifi_l2_*` labels; those must not replace per-study L2 in any
downstream analysis.

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
starts with pseudobulk composition and a gene-presence UpSet plot; if the
single-cell branch is enabled it appends its embedding, depth checks, UMAPs of
log(UMIs per cell) and mitochondrial percentage (capped at 15%), and the matrix
comparing experimental merged AIFI L2 calls to each cell's per-study label.
The report also accounts for raw counts retained or discarded by each
gene-symbol join. The pseudobulk outer union retains all input counts; the
technical-covariate plots show the number of studies for each technology,
intronic-count, and frozen-status value. The optional single-cell inner
intersection includes a per-study bar chart of the fraction of counts discarded,
alongside a table with count totals and the gene symbol with the largest loss per
study.
The merged artifact retains per-study labels in `aifi_l2_majority` for all
downstream grouping and splitting. Diagnostic re-annotations are stored in
`experimental_aifi_l2_majority` (Harmony graph) and
`experimental_aifi_l2_unintegrated_majority` (unintegrated PCA graph); the QC
report shows their UMAPs and label-concordance heatmaps.
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
| One cell-type report template | `make run-cell-type-test CELL_TYPE=cd14-monocyte` | Use `make render-cell-type-test` or `make render-cell-type-reports-test` when existing analysis artifacts are sufficient. |
| Rerender reports from existing analysis artifacts | `make render-cell-type-reports[-test]` | Use `make render-cell-type[-test] CELL_TYPE=<slug>` to refresh only one report. |
| Cell-type splitting or shared analysis | `make run-test`, then `make run-cell-type-analysis-test` | Use `make run-all-test` when the complete ordered test workflow is needed. |
| Positive pseudobulk DE fitting | `make run-test`, then `make run-de-synthetic-test` and `make check-synthetic-de-results` | Follow with `make run-cell-type-analysis-test` when checking that DE results reach reports. |
| Core workflow for one or a few studies | `make run-test STUDIES=<list>` | Omit `STUDIES` to request all configured studies. |
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

Study-QC, merge-QC, and cell-type report templates are maintained as Jupytext
Python sources (`reports/qc_report.py`, `reports/merge_qc_report.py`, and
`reports/cell_type_report.py`). Their report runners materialize a temporary
notebook immediately before execution; the generated template notebooks are
intentionally not tracked. The published `executed.ipynb` is the executed
report artifact.

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

## Per-cell-type residual-variation reports

The optional downstream workflow loads a completed merged single-cell H5AD once
and splits it by the retained per-study `aifi_l2_majority` call. Experimental
merged labels are never used for this split. For each type, `ANALYSE_CELL_TYPE`
creates the type-specific PCA/UMAP/Leiden embedding, marker tables, PC--age
correlation table, and `analysis.h5ad`. `RENDER_CELL_TYPE_REPORT` then loads
that analysis artifact plus the primitive raw-count H5AD to create the executed
notebook and HTML. It also plots the retained per-study `aifi_l3_majority`
labels; no merged L3 prediction is created. Global neighbours are deliberately
not copied into the split artifacts.

Each report also plots the type's sample-level fraction across age, coloured by
study. It shows two descriptive denominators, both recorded before splitting:

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

The combined view fits a weighted least-squares model with study-specific
intercepts and a shared age slope. Each study receives equal total weight; the
plot shows the prediction averaged equally across studies, with a 95% HC3
confidence interval. Shading marks the age range observed in every study;
predictions outside that range extrapolate for at least one study. This is a
study-adjusted descriptive summary, separate from the per-study fraction models.
The HC3 interval treats sample rows as independent; repeat samples from the same
subject are not clustered, which may understate uncertainty.

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

The separate `differential_expression.nf` workflow consumes the merged
sample-level pseudobulk H5AD. It first lists retained per-study
`aifi_l2_majority` labels, then runs one independent task per label. Each task
loads the same small merged pseudobulk H5AD and subsets its label in memory; it
does not materialize a per-type pseudobulk H5AD. Each task runs one `~ age +
sex` fit per study and one shared age-and-sex `~ study + age + sex` fit across
studies. When present and variable, BMI and CMV each receive separate per-study
and combined fits, with age and sex as adjustment covariates. Combined fits
include studies with usable values for that covariate and report the studies
and sample counts used. If only one study contributes, the fit drops the study
term and reports that study's adjusted association. Age, sex, BMI, and CMV
effects each have a combined-model volcano plot where estimable. A final collector writes the result
directories and manifest. Per-study results also show significant-gene counts
inside and outside the per-covariate intersection of genes tested by all
available study fits. Per-study fits use genes available in that study; each combined fit
uses the intersection of genes available across its included studies, so
study-absent genes' synthetic outer-join zeros are excluded. The parallel
type-fitting task is configured for one CPU, which is passed directly to
PyDESeq2; all other CPU, time, and memory requests use the executor defaults.

Samples must be age 20 or older and have at least 10 cells in that
sample × cell-type pseudobulk. Samples with missing/unknown sex metadata or
zero counts over the fit's gene universe are also omitted. All sample
exclusions and model fits that cannot estimate the configured covariates are
recorded in the JSON manifest. CSV results and the manifest are published
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

The small end-to-end core fixture is intentionally allowed to produce no DE
fits: bypassing the cell cutoff does not create additional independent samples
or fix a rank-deficient age/sex design. `make run-de-synthetic-test` creates a
small deterministic pseudobulk H5AD directly (the fixture is generated, not
checked in, and does not pass through single-cell pseudobulking). It contains
three synthetic studies, eight independent samples per study and cell type,
and at least 15 cells per pseudobulk. Its negative-binomial counts include
known synthetic age, sex, BMI, and CMV signals and study-specific gene-
availability flags. BMI is present in two studies; CMV is present in one.
Missing values exercise covariate-specific study selection and complete-case
sample counts, including the single-study combined-fit path. This exercises the
normal age and cell-count filters and produces positive per-study and combined
fits without copying the large real pseudobulk file. Every result is labeled
synthetic and test-only. This fixture is realistic in workflow shape, not in
statistical complexity: it has balanced sample counts, a small number of planted
effects and gene-availability gaps, and no modeled cross-gene correlation or
study-specific effect heterogeneity. It validates fitting and report wiring,
not biological power or expected real-data results. Do not upsample real test
samples to simulate replication; more cells within a sample do not add
independent samples.

To render the standard test reports with the synthetic DE results, run:

```bash
make run-test
make run-de-synthetic-test
make run-cell-type-analysis-test
```

`make run-all-test` is the complete ordered Docker test. Its core and cell-type
analysis stages use available downsampled real-data fixtures; its positive DE
stage uses the separate synthetic pseudobulk. The final reports combine the
real test merge/splits with those synthetic DE results to check the handoff and
plots. `make verify` runs this same target after linting, unit tests, Nextflow
lint, and documentation checks. Pass `STUDIES=<list>` to narrow the core test
inputs; it does not change the synthetic DE fixture.

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

By default this reads `output/merged/pseudobulk_merged.h5ad`. Override the
input or output roots with `PSEUDOBULK_INPUT=/path/to/pseudobulk_merged.h5ad`
or `OUTDIR=/path/to/output`. Override where cell-type report generation looks
for DE results with `DE_RESULTS_DIR=/path/to/differential_expression`.
`make run-de-test` requires an existing test pseudobulk merge; create it first
with `make run-test`. It uses the explicitly labeled test mode described above.
Use `make run-de-synthetic-test` to exercise successful PyDESeq2 fits with the
production sample-inclusion filters enabled.

DE fits estimate age and sex together per study and in a study-adjusted combined
model. BMI and CMV are fitted separately when metadata support them, adjusting
for age and sex; missing values are excluded from the relevant fit. A single-
study optional-covariate fit omits the study term and still produces its
volcano plot. CMV negative/positive values are harmonized to no/yes. Combined
volcano plots list the studies and samples that contributed. Cell-type fraction
forest plots label categorical contrasts explicitly: male versus female, and
yes CMV versus no CMV.

## Further reading

See [INPUT_FILES.md](INPUT_FILES.md) for exact external paths and checksums,
[IMPLEMENTATION.md](IMPLEMENTATION.md) for architecture and output contracts,
and [PLAN.md](PLAN.md) for open scientific and reproducibility work.

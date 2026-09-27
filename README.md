# PBMC ageing scRNA-seq pipeline

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

## Current status

- The complete Docker test workflow has passed with every configured
  preparation adapter, per-study harmonization/QC, pseudobulks, and both merge
  branches (29 successful processes).
- Executed QC notebooks can produce self-contained per-study HTML reports.
- The full datasets, including the real Wang25 RDS conversion, have not been
  run in this workspace.
- A DSL2 Nextflow workflow connects conversion, preparation, harmonization,
  merge, and QC, with a test profile for the 200-cell inputs.
- Both Docker images build successfully. A synthetic Seurat object passed the
  RDS-to-H5AD conversion with counts, names, and metadata intact.
- Pseudobulk aggregation and merge are part of the production Nextflow graph;
  merge QC reports include study/sample contribution, AIFI-L2 overlap, depth,
  gene coverage, and (when applicable) merged UMAP checks.
- Per-cell-type analysis and report rendering are separate workflow stages.
  The focused CD14-monocyte test has produced the derived analysis H5AD, report
  tables, executed notebook, and self-contained HTML report.
- Python image dependencies are pinned in `requirements.lock`; preparation also
  publishes retained source-`obs` metadata sidecars for adapter auditing.

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
                                                        └── executed.ipynb + tables
                              output/cell_type_analysis/cell_type_manifest.json
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

| Study | Subject | Sample | Sampling timepoint | Technical batch | Notes |
|---|---|---|---|---|---|
| AIDA25 | `donor_id` | country-qualified donor ID | not provided | supplementary experimental batch | Lonza material replicated across sites remains intentionally split until downstream exclusion. BMI uses the supplementary donor metadata. |
| AIFI | specimen GUID prefix | specimen GUID | `sample.visitName` | `well_id` | A specimen GUID may span multiple wells; pseudobulk combines those technical partitions. |
| OneK1K | `donor_id` | `donor_id` | not provided | `pool_number` | No separate sampling-timepoint identifier is available. |
| Terekhova23 | `Donor_id` | `Tube_id` | workbook visit | barcode-derived batch | Each tube maps to one donor and visit; a tube may span technical batches. |
| Wang25 | `Sample ID` | `Sample ID` | not provided | not provided | No separate donor, timepoint, or technical-batch identifier has yet been recovered. |
| Nehar-Belaid26 | Supplementary Data 1a `IDs` | Supplementary Data 1a `Names` | not provided | Supplementary Data 1a `runs_10x` | Reads the supplied GEO raw 10x tar directly (without extraction); `JB` subject IDs parsed from member names join to Supplementary Data 1a `IDs`, which supplies sample IDs and donor metadata. Published labels transfer only on a unique sample ID plus normalized barcode-sequence match. |
| Fachrul26 | `donor_id` | `sample_id` | not provided | `library_id` | Uses embedded H5AD metadata; `feature_name` supplies gene symbols. |
| Perez22 | `donor_id` | `sample_uuid` | not provided | `library_uuid` | Uses embedded H5AD metadata; `feature_name` supplies gene symbols. `obs.disease` is normalized and only `normal`/healthy cells are retained; embedded SLE cells are excluded before harmonization. |

`tests/` contains only automated test code. Reproducible local intermediates
are collected under `cache/` for RDS-to-H5AD conversions. Downsampled smoke-test
fixtures are generated under the ignored `test_data/` directory.

## Local Python setup

`requirements.lock` pins the Python 3.12 runtime and report environment used
by the Docker image. `requirements.dev.lock` extends that exact base with the
test and lint tools used in the local `.venv`. Docker installs the runtime lock
before the project package, so source-only image rebuilds reuse the dependency
layer. Update the locks only as a deliberate dependency change and review their
complete diffs alongside `pyproject.toml`. Use `make install` for the local
development environment.

```bash
python -m venv .venv
. .venv/bin/activate
make install
```

## Run harmonization

```bash
# Complete all downsampled test studies and generate QC reports
make run-test STUDIES=all

# Recreate downsampled inputs from locally available full inputs
python scripts/create_test_data.py --cells 200
```

Each study is prepared before harmonization. `pbmc-prepare` materializes a
gzip-compressed CSV with exactly one canonical metadata row for every expression
cell and a retained-source-`obs` sidecar with the same `cell_id` index.
`pbmc-harmonize --prepared-obs <study>.cells.csv.gz` then attaches only the
canonical artifact to raw counts. This is the debugging boundary used by
Nextflow; the sidecar exists solely to audit source-to-canonical mappings.

## Download public inputs

`make download-inputs STUDIES=aifi` runs the opt-in, best-effort downloader.
It never runs automatically with analysis commands, does not authenticate or
scrape portals, and records outcomes in ignored
`input_data/download_manifest.json`. Sources requiring terms acceptance or a
manual portal download are reported with instructions rather than treated as a
pipeline failure. The Nehar-Belaid26 raw 10x archive, label H5AD, and Supplementary Data 1
are also direct, verified downloads, so a fresh input can be acquired with:

```bash
make download-inputs STUDIES=nehar_belaid26
```

The remaining study inputs are explicitly manual-only until stable, verified
file URLs are recorded. See `config/input_sources.json` and [INPUT_FILES.md](INPUT_FILES.md).

## Convert a Seurat RDS

The converter requires R, Seurat, reticulate, and Python anndata. It processes exactly one object
and retains only raw counts plus observation metadata:

```bash
Rscript scripts/convert_rds.R \
  --input input_data/wang25/scRNA-seqProcessedLabelledObject.rds \
  --output cache/converted/scRNA-seqProcessedLabelledObject.h5ad
```

Nextflow invokes this for studies with a `conversion` entry in
`config/studies.json`. At present that is Wang25. Test runs use its downsampled
RDS fixture and therefore exercise the same conversion process.

## QC reports

Passing `--qc` executes the tracked `reports/qc_report.ipynb` template after each
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
The merged artifact retains per-study labels in `aifi_l2_majority` for all
downstream grouping and splitting. Diagnostic re-annotations are stored in
`experimental_aifi_l2_majority` (Harmony graph) and
`experimental_aifi_l2_unintegrated_majority` (unintegrated PCA graph); the QC
report shows their UMAPs and label-concordance heatmaps.
`output/run_manifest.json` records selected studies, image references,
configuration checksum, and checksums of the merged deliverables.

## Nextflow and Docker

Nextflow >=25.04 is the canonical workflow runner. It reads studies and inputs from
`config/studies.json`, runs independent studies in parallel, caches successful
processes, and publishes the final H5AD, JSON, and QC artifacts. Limit a run
with a comma-separated study list:

```bash
make run-test STUDIES=wang25
make run STUDIES=aifi,onek1k
```

The `test` profile lowers resource requests and uses isolated 200-cell fixtures
from `test_data/`. RDS-based studies retain an RDS fixture and therefore still
pass through the same `CONVERT_RDS` process as a production run. Create or refresh fixtures
from the production inputs once with `make test-data`; this never modifies
`input_data/`. Override the fixture size with `make test-data TEST_CELLS=500`
when needed. `--skip_qc` stops after the H5AD merge artifacts are created (it
skips only notebook/HTML QC). By default the workflow does not create a
whole-dataset single-cell merge, because all studies must fit in RAM at once.

For faster development, `make run-test` and `make run-cell-types-test` default
to every configured study except `nehar_belaid26`, whose much larger fixture
dominates test runtime. This is Make-level test selection only: it does not
change the configured studies or production runs. Pass `STUDIES=all` (or an
explicit comma-separated list) to override it. `make verify` deliberately
restores all configured studies by default because it is now an occasional
broad check; pass an explicit `STUDIES` list if a narrower verification is
needed.

Choose the smallest check that covers a change during development. Python-only
changes normally need `make lint test-unit`; Nextflow or configuration changes
need `make workflow-lint` plus a focused `make run-test STUDIES=<study>` when
the wiring needs execution. For one cell-type report, reuse its split with
`make run-cell-type-test CELL_TYPE=cd14-monocyte`; when only the report
template changed, reuse its published analysis with
`make render-cell-type-test CELL_TYPE=cd14-monocyte`. Run `make verify` for
broad cross-cutting or release-level validation, rather than after every
focused change. Run `make docs-check` after changing tracked Markdown files or
Make target names.

Harmony integration uses the pinned Harmony2 C++ backend. It uses the backend's
native thread policy; the legacy OpenMP/OpenBLAS thread overrides that crashed
on large production runs have been removed.

Cell-type reports show signed marker log-fold changes (higher and lower genes),
direct cluster labels on the local UMAP, compact study metadata coverage, and
up to ten local PCA loading panels. They include native local PC1-versus-PC2
and PC3-versus-PC4 scatterplots, with the standard study palette as colour and
the reported `technology` field as marker shape. They also include cumulative
variance explained and descriptive age correlations for the first ten PCs.
PCA loading figures use explicit subplot spacing because Scanpy's composite
figure is not compatible with Matplotlib `tight_layout()`.

The deterministic `make run-test` workflow enables that branch by default;
production runs remain opt-in.
Enable that separate path on a runner sized for the selected inputs:

```bash
make run STUDIES=all MERGE_SINGLE_CELL=true
```

The merge and cell-type report templates are maintained as Jupytext Python
sources (`reports/merge_qc_report.py` and `reports/cell_type_report.py`). Their
report runners materialize a temporary notebook immediately before execution;
the generated template notebooks are intentionally not tracked. The published
`executed.ipynb` is instead the executed report artifact. The separate
`reports/qc_report.ipynb` study-QC template remains hand-maintained.

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
from tiny sample compartments. Study-specific linear trends are shown as small
multiples; pooled quadratic and LOWESS curves are additional, explicitly
study-unadjusted descriptive views.

```bash
make run-cell-types
```

This is a downstream workflow: it requires an existing merged single-cell H5AD
at `output/merged/single_cell_merged.h5ad` and does not run the merge automatically.
Override that location with `MERGED_INPUT=/path/to/single_cell_merged.h5ad`.
`make run-cell-types-test` first
creates the small core test merge and then runs the downstream workflow under
`output/test/`. All test targets publish only beneath `output/test/`; production
outputs remain under `output/` unless `OUTDIR` is explicitly overridden.

The raw-count splits are published under `<outdir>/cell_type_splits/` so one
report can be rerun without re-splitting the merged H5AD:

```bash
make run-cell-type CELL_TYPE=cd14-monocyte
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
```

The targeted test command deliberately does not build a core test merge or all
reports; it requires the corresponding split to already exist. It runs the
analysis and rendering stages for that one type. Create only the reusable splits
from an existing merge with `make split-cell-types-test`, or run the full
`make run-cell-types-test` when a fresh merge is needed. A targeted rerun
refreshes that type directory but not `cell_type_manifest.json`; rerun the
all-type workflow before treating the published report set as a new
manifest-backed release. `cd14-monocyte` is the default targeted test type
because it is well represented in the bundled multi-study fixture.

For a template-only report edit, `make render-cell-type[-test]` instead
requires the published split plus `analysis.h5ad` in that type's analysis
directory. It schedules only rendering and overwrites the executed notebook and
HTML without rerunning local PCA, Harmony, UMAP, clustering, or marker ranking.

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

Build and use both containers with:

```bash
make images
make test-data
make run-test STUDIES=wang25
make run STUDIES=all
make run-cell-types
make run-cell-types-test
make split-cell-types-test
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
make docs-check
```

`make run`, `make run-test`, and `make run-no-qc` build the required
content-tagged local images automatically whenever `NF_PROFILE` contains
`docker`; use `make images` to build both without launching a workflow. Python
and R images are tagged independently from the runtime files copied by their
respective Dockerfiles. Report templates and models are staged as task inputs,
so report-only or documentation-only edits do not retag the Python image or
invalidate unrelated workflow cache entries. Set `IMAGE_TAG` explicitly to use
one common release tag for both images.
They also select the matching Nextflow resource profile (`core` or
`cell_type_analysis`), keeping process selectors scoped to the workflow that
defines them. For direct Nextflow invocation, include the corresponding profile.

For a production run whose repository and task filesystem are separate, keep
the repository path readable by Docker and choose an external work/output
location:

```bash
make run NF_PROFILE=docker STUDIES=all \
  WORK_DIR=/path/to/run/work OUTDIR=/path/to/run/output
```

Nextflow stages each declared task input into the Docker work directory. Keep
the repository readable to Nextflow, and ensure `WORK_DIR` supports symbolic
links; do not rely on arbitrary host repository paths being available inside a
container.

The Python image is defined by `docker/python.Dockerfile`; the conversion image is defined by
`docker/r-conversion.Dockerfile`. The latter starts from a digest-pinned Seurat
image and adds the Python/anndata bridge. Metadata and externally
supplied CellTypist models remain outside the images and are declared as
Nextflow task dependencies. Expression files and generated artifacts are staged
through Nextflow. The Python image and a complete Wang25 test workflow have
been verified. The R conversion image is defined separately because it is much
larger and only needed for RDS inputs.

The Python image installs the reviewed `requirements.lock` before application
source, then installs the project with `--no-deps --no-build-isolation`.
It contains package code and workflow configuration; task-specific report
templates and models are staged by Nextflow. This preserves the pinned
dependency layer during source-only changes without allowing the project install
to resolve a different environment.

The content-derived image tag participates in Nextflow task identity, so a
changed build context automatically receives separate resumable cache entries.

## Make shortcuts

The `Makefile` is intentionally a thin convenience layer rather than a second
workflow. Run `make help` to see its targets. Common commands include:

```bash
make lint
make workflow-lint
make docs-check
make test
make validate-test STUDIES=wang25
make validate-full STUDIES=terekhova23
make run-test STUDIES=wang25
make run-no-qc STUDIES=wang25
make run-test STUDIES=all MERGE_SINGLE_CELL=true
```

Use `make verify` when a full, fresh Docker validation is warranted; it is not
the routine command for an isolated edit.

For a non-Docker Nextflow run, activate `.venv` so the installed
`pbmc-harmonize`, `pbmc-qc`, `pbmc-merge`, `pbmc-merge-qc`, and `pbmc-manifest`
commands are on `PATH`. The
Docker profile does not require the Python environment on the host.

See `input_data/README.md` for the expected input layout and `PLAN.md` for the
remaining reproducibility work.

# PBMC Ageing Harmonization

Configuration-driven processing of PBMC ageing scRNA-seq studies. It produces
one harmonized H5AD per study, per-study sample × AIFI-L2 pseudobulks, and a
cross-study pseudobulk matrix. An explicitly enabled pathway also creates a
raw-count cross-study single-cell merge. A separate, opt-in integration
benchmark compares unintegrated PCA, Harmony, and BBKNN representations and
their method-specific cell graphs.

Study-specific inputs and adapter dependencies live in `config/studies.json`;
the tracked `config/studies.schema.json` and
`config/input_sources.schema.json` enforce their document structures before a
workflow starts. Explicit adapters under `src/pbmc_pipeline/studies/` own
source-specific joins and reshaping. The output metadata contract lives in
`config/harmonized_obs_schema.json`. See `IMPLEMENTATION.md` for the high-level
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
project downloader does not support. AIDA25, Terekhova23, Terekhova26, and
Wang25 also have manual-only inputs; `STUDIES=all` skips studies with
incomplete inputs and records the decision in `output/run_manifest.json`.

After the required AIFI models and the inputs for the studies you want to
include are in their configured paths, run:

```bash
make run-all STUDIES=all
```

This builds the workflow images when needed, runs the core harmonization and
single-cell merge, fits pseudobulk differential-expression models, profiles the
top age-associated genes for each cell type, and renders separate per-type DE
and single cell analysis reports.
It uses the Docker profile by default; no local Python or R installation is
required.

## Test quick start

Run this before a production workflow: it confirms that Docker, Nextflow,
models, and the available inputs work before committing to the much longer full
analysis. It creates missing test fixtures for locally available inputs, then
runs the complete ordered Docker test workflow:

```bash
make test-data
make run-all-test
```

By default this uses the representative `aida25,cima,terekhova26` set and
skips unavailable fixtures. It covers a supplementary-metadata adapter, a
counts-layer input, and a raw-counts input. Use `STUDIES=all` when a test must
exercise every available fixture.
It runs the core pipeline on available downsampled real-data fixtures, checks
the integration benchmark, fits and checks a generated pseudobulk fixture for
positive DE coverage, profiles two top age-associated genes per cell type,
smoke-tests a focused gene report, and renders single cell analysis reports
from the real merge plus per-type DE reports from the synthetic fixture. The
benchmark artifact retains cell labels, UMAP coordinates, and neighbour graphs
without raw counts. Use
`STUDIES=<list>` to select a different core-test set. Existing fixtures are preserved; use
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

Run the focused target in [Choosing a Make target](#choosing-a-make-target)
before relying on a changed workflow or report.

## Data flow

```text
input_data/ expression objects + metadata
                │
                ├── Seurat RDS ──> converted H5AD
                │
                ▼
 prepare canonical metadata plus a retained source-obs audit sidecar
                │
                ├──> output/harmonization/prepared/<study>.cells.csv.gz
                └──> output/harmonization/prepared/<study>.source_obs.csv.gz
                ▼
      harmonize one configured study
                │
                ├──> output/harmonization/harmonized/<study>.h5ad
                ├──> output/harmonization/provenance/<study>.json
                └──> output/harmonization/harmonization_qc/<study>/{executed.ipynb,report.html}
                │
                ▼
  pseudobulk each study by sample/AIFI-L2
                │
                └──> output/harmonization/pseudobulk/<study>.pseudobulk.h5ad
                ▼
       output/harmonization/merged/pseudobulk_merged.h5ad
                │
                ├── optional: output/harmonization/merged/single_cell_merged.h5ad
                └──> output/harmonization/harmonization_qc/merged/{executed.ipynb,report.html}

optional global integration benchmark
  single_cell_merged.h5ad ──> output/harmonization/integration_benchmark/
                                └──> embeddings + graphs + labels H5AD + JSON + report.html

optional downstream cell-type workflow
  single_cell_merged.h5ad ──> backed split by per-study AIFI-L2
                                (materialize and write one type at a time)
                                ├──> output/analysis/per_type/<type>/split.h5ad
                                └──> analyse each type
                                      └──> analysis.h5ad + tables
                                            └──> render notebook/HTML
                                                  └──> output/analysis/per_type/<type>/analysis/
                                                        ├── report.html
                                                        └── executed.ipynb + tables
                              output/analysis/per_type/cell_type_manifest.json
                              output/index.html (links to published harmonization and analysis reports)

optional pseudobulk differential-expression workflow
  pseudobulk_merged.h5ad ──> list per-study AIFI-L2 labels
                              └──> one parallel task per label, each loading and
                                   subsetting the same small pseudobulk H5AD in memory
                                     ├──> per-study fits with site adjustment
                                     └──> merged ~ study_site + age + sex + log10_total_counts fit
                                           ├──> output/analysis/per_type/<type>/differential_expression/report.html
                                           ├──> output/analysis/combined/differential_expression/
                                           └──> output/analysis/combined/trajectory/report.html

Focused gene profiles are stored under output/analysis/per_gene/<type>/<gene>/
```

All immutable local inputs belong under `input_data/`. Generated conversions,
reports, and harmonized outputs are kept outside it and ignored by Git.
This includes the pipeline-wide GO enrichment files in `input_data/go/`, as
well as study-specific data.

Published outputs are organized under `output/harmonization/` for preparation,
harmonized data, pseudobulk, merged inputs, QC reports, and integration
benchmarks. `output/analysis/` contains `combined/`, `per_type/`, and
`per_gene/` results.

Per-cell preparation artifacts are generated under
`<outdir>/harmonization/prepared/`; they are reproducible from the expression
object plus the declared study dependencies.
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
are removed. For CellTypist studies, harmonization writes this count object to
the final H5AD before normalizing `.X` in memory, then updates metadata and
embeddings on disk while preserving the saved counts.

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
reported. Record the source field in `provenance`, and use the lowest-level
run, pool, library, or processing-group identifier that assigns one value per
sample and groups multiple samples where possible. Record
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

`subject` is a biological individual. `sample` is one biological specimen;
it is the unit aggregated into sample × AIFI-L2
pseudobulks. `batch_single_cell` records the lowest available technical
grouping above samples that can be assigned consistently: each sample maps to
one value, and a useful grouping is shared by multiple samples. Examples
include multiplexed 10x processing runs and processing batches. A sample may
not map to multiple technical batches. If a source value is missing or
ambiguous for a sample, that sample receives a sample-specific fallback. If a
source only provides one study-wide value or a unique value per sample,
samples are used as their own batches so the field does not simply mirror the
study. The preparation report records how many batch groups and sample-specific
fallbacks were written. Study prefixes keep source batch IDs distinct in merged
data. This field describes technical correlation groups; it never defines a
pseudobulk. Missing string metadata in other fields is written as
`not_provided`; missing numeric metadata is `NaN`.

`study_celltype` preserves the highest-resolution published or author-provided
cell-type annotation available per cell. An unsupervised cluster is not used
merely because it has more values. The selected fields below are source labels
for audit and reporting; cross-study analyses use the separately harmonized
AIFI labels.

Every harmonized and merged H5AD also contains `obs.study_doi`. It is a
categorical per-cell DOI drawn from the required `publication_doi` study
configuration field; it is intentionally not duplicated in `uns`.

| Study | `study_celltype` source | Selection basis |
|---|---|---|
| AIDA25 | `Annotation_Level4` | Highest level of the published AIDA hierarchy. |
| AIFI | `AIFI_L3` | Highest AIFI level supplied in the source H5AD. |
| CIMA | `cell_type_l4` | Highest CIMA hierarchy level. |
| Fachrul26 | `cell_type` | One-to-one with `author_cell_type`; no finer author annotation is supplied. |
| Nehar-Belaid26 | `LS_L4` | Highest published LS hierarchy level. |
| OneK1K | `cell_type` | One-to-one with the supplied `predicted.celltype.l2`; no finer level is supplied. |
| Perez22 | `cell_type` | One-to-one with `author_cell_type`; no finer author annotation is supplied. |
| Terekhova23 | `Cluster_names` | Only supplied per-cell type-label field. |
| Terekhova26 | `author_celltype` | Highest author label; identical to `annotation`. |
| Wang25 | `secondary_type` | Finer of the supplied `primary_type`/`secondary_type` fields. |

| Study          | Reference DOI                                                                | Subject                     | Sample                        | Technical batch                  | Notes                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| -------------- | ---------------------------------------------------------------------------- | --------------------------- | ----------------------------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| AIDA25         | [10.1016/j.cell.2025.02.017](https://doi.org/10.1016/j.cell.2025.02.017)     | `donor_id`                  | country-qualified donor ID    | supplementary experimental batch | Source `obs.assay` reports 10x 5′ v2. Lonza controls (`DCP_ID` prefixed `LONZA`) are reported as European and retained in harmonized and pseudobulk data for calibration, but excluded from DE models by subject ID. Repeated donor-country records are deduplicated in the supplementary metadata join and their batch is marked unclear. BMI uses the supplementary donor metadata.                                                                                                                                                                                                                                                        |
| AIFI           | [10.1038/s41586-025-09686-5](https://doi.org/10.1038/s41586-025-09686-5)     | specimen GUID prefix        | specimen GUID                 | `pool_id` (`aifi_` prefix)       | The configured Immune Health Atlas reference has one sample per donor. Each `pool_id` combines hashed samples and is distributed across chips and wells; a sample maps to one pool, while it can span wells. We use this lowest shared grouping for technical batch. A broader `batch_id` groups the processing occasion (up to 24 samples across multiple pools), not an exact timestamp; when the same aliquots are processed for scRNA-seq and flow cytometry on the same day, they share a Batch ID. The source fixture confirms one pool per specimen and multiple specimens per pool. Although the source includes `sample.visitName`, it is not retained as a separate harmonized timepoint field. The source H5AD has no assay field. ([AIFI pipeline](https://apps.allenimmunology.org/aifi/pipelines/multiplexed-scrna-seq/)) |
| CIMA           | [10.1126/science.adt3130](https://doi.org/10.1126/science.adt3130)           | `obs.sample`                | `obs.sample`                  | `obs.RNA_experiment`             | The H5AD embeds the retained sample and experiment fields; `RNA_experiment` groups samples by date-coded RNA processing experiment, while `library` can vary within a sample. Each of the 408 retained sample IDs maps to exactly one of 51 RNA experiment groups in local metadata. The separately supplied metadata table is deliberately not a workflow input. Native `cell_type_l4` (70 labels) is retained as the per-study label. Paper methods report DNBelab C Series, STAR v2.7.2b, GRCh38, and exclusion of intronic reads; no demultiplexing method is described, so `none` is inferred. Frozen status is not reported.                                                            |
| OneK1K         | [10.1126/science.abf3041](https://doi.org/10.1126/science.abf3041)           | `donor_id`                  | `donor_id`                    | `pool_number`                    | No separate specimen identifier is present in the source H5AD; the sample-per-donor assumption remains under review.                                                                                                                                                                                                                                                                                                             |
| Terekhova23    | [10.1016/j.immuni.2023.10.013](https://doi.org/10.1016/j.immuni.2023.10.013) | `Donor_id`                  | `Tube_id`                     | barcode-derived batch            | Only each donor's earliest available visit is retained before harmonization; later-visit cells are excluded. Tubes spanning multiple barcode-derived batches receive sample-specific fallback values to preserve one batch assignment per tube. The paper explicitly reports frozen input material.                                                                                                                                                                                                                                    |
| Terekhova26    | [10.1016/j.immuni.2026.09.012](https://doi.org/10.1016/j.immuni.2026.09.012) | `obs.donor_id`              | `obs.donor_id`                | `obs.Pool`                       | The supplied author-reannotated KCL H5AD embeds all runtime metadata. Native `author_celltype` (59 labels) is retained as the per-study label. Genetic demultiplexing (`souporcell_donor`) is paired with HTO antibodies that map genotypes to samples. Paper methods report cryopreserved/thawed PBMCs, Cell Ranger multi v9.0.1, and GRCh38-2020-A; intronic reads are inferred included from the Cell Ranger v6+ default arguments. |
| Wang25         | [10.1038/s41590-024-02059-6](https://doi.org/10.1038/s41590-024-02059-6)     | `Sample ID`                 | `Sample ID`                   | sample-specific fallback        | No separate donor or technical-batch identifier has yet been recovered, so each sample is its own batch. Fresh/frozen status is not provided in the paper or supplementary metadata; the author was contacted but did not respond.                                                                                                                                                                                                                           |
| Nehar-Belaid26 | [10.1038/s41467-026-73729-2](https://doi.org/10.1038/s41467-026-73729-2)     | Supplementary Data 1a `IDs` | Supplementary Data 1a `Names` | Supplementary Data 1a `runs_10x` | Reads the supplied GEO raw 10x tar directly (without extraction); `JB` subject IDs parsed from member names join to Supplementary Data 1a `IDs`, which supplies sample IDs and donor metadata. Published labels transfer only on a unique sample ID plus normalized barcode-sequence match.                                                                                                                                            |
| Fachrul26      | [10.64898/2026.02.15.704933](https://doi.org/10.64898/2026.02.15.704933)     | `donor_id`                  | `sample_id`                   | `library_id`                     | Paper methods state PBMCs were cryopreserved and thawed; this takes precedence over the conflicting CELLxGENE `obs.sample_preservation_method=fresh`.                                                                                                                                                                                                                                                                                  |
| Perez22        | [10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970)           | `donor_id`                  | `sample_uuid`                 | `library_uuid`                   | Paper reports USA frozen PBMCs processed with 10x 3' v2 and Cell Ranger (version unspecified); intronic-read handling is not reported. Uses embedded H5AD metadata; `feature_name` supplies gene symbols. `obs.disease` is normalized and only `normal`/healthy cells are retained; embedded SLE cells are excluded before harmonization.                                                                                                  |

`tests/` contains only automated test code. Reproducible local intermediates
are collected under `cache/` for RDS-to-H5AD conversions. Downsampled smoke-test
fixtures are generated under the ignored `test_data/` directory.

## Download inputs

> **Strictly opt-in and potentially large.** Only `make download-inputs`,
> `make download-go-inputs`, and `make download-inputs-strict` retrieve data; `make install`, test targets,
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
make download-go-inputs
```

Study selections include the small pipeline-wide GO enrichment inputs in
`input_data/go/`; use `make download-go-inputs` to retrieve only those files.
Their source files can change upstream, so the
registry intentionally has no fixed checksum or release requirement; the
download manifest and each report's enrichment provenance record the exact
SHA-256 values used, along with the GO ontology `data-version`.

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
present at production paths are test-scale subsets. Entries explicitly marked
as allowing a local subset retain their published full-file size as download
documentation; a normal download check reports `local_subset` rather than a
size failure, while the strict preflight still rejects it. Confirm that acquired
files are suitable for production analysis.
See `config/input_sources.json` and [INPUT_FILES.md](INPUT_FILES.md).

## QC reports

Passing `--qc` executes the tracked `reports/qc_report.py` template after each
harmonized H5AD is written. Each study gets a separate directory:

```text
output/harmonization/harmonization_qc/<study>/
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
  --input output/harmonization/harmonized/onek1k.test.h5ad \
  --run-report output/harmonization/provenance/onek1k.test.json \
  --output-dir output/harmonization/harmonization_qc/onek1k.test \
  --study onek1k
```

The merge workflow renders one combined document at
`output/harmonization/harmonization_qc/merged/`. It
covers pseudobulk composition, technical covariates, and a covariate-completeness
heatmap showing the fraction of pseudobulk profiles missing each covariate per study. String
metadata uses the canonical `not_provided` value for both source missingness and
explicitly unreported fields, while numeric metadata uses `NaN`. The pseudobulk outer
union retains all input counts. Shared-gene availability, hypothetical count
loss under the all-study intersection, and the highest-count excluded gene are
reported with the cross-cell-type expression atlas, where that intersection is
actually used. Cross-study UMAP and benchmark AIFI-L2 comparisons belong to the
separate integration benchmark; they are computed on the cells selected for
that benchmark and do not provide a full-data embedding if the configured cap
is exceeded.
The merged artifact retains per-study labels in `aifi_l2_majority` for all
downstream grouping and splitting. Run `make run-integration-benchmark` after
the core workflow to create a diagnostic H5AD with unintegrated PCA, Harmony,
and BBKNN representations, their neighbor graphs, and CellTypist labels
majority-voted on each graph. Scanorama is temporarily disabled in the active
benchmark configuration because its correction step can require large dense
temporary arrays. By default, the benchmark processes the
full merged H5AD unless it exceeds `integration_benchmark.max_cells`, currently
20,000,000. Above the cap, it selects cells uniformly at random without
replacement using `processing.random_seed` (42 by default); that selection is
shared by all methods. Set the cap in `config/pipeline.json` to change this
behavior. The report and JSON sidecar record source and selected cell counts,
fraction, strategy, and seed. The input is opened in backed mode, and only the
selected rows are materialized. The published artifact does not duplicate the
merged expression matrix. Graph edge overlap and per-study versus benchmark
AIFI-L2 label concordance are diagnostic comparisons, not acceptance
thresholds.
Timestamped benchmark stage progress is captured in the Nextflow task log
(`.command.out` and `.command.err` in that task's `work/` directory), like
other process output. The terminal shows Nextflow task progress. Set
`INTEGRATION_BENCHMARK_LOG_LEVEL=DEBUG` for additional diagnostics. Production
and test targets use the standard Python image.
Its `report.html` has separate embedding-integration and type-calling sections,
including graph edge overlap and per-study-versus-benchmark AIFI-L2 concordance.
`output/run_manifest.json` records requested, selected, and skipped studies,
image references, configuration checksum, and checksums of the merged
deliverables.

## Nextflow and Docker

Nextflow runs the configured studies in parallel and caches completed tasks.
The Docker profile is the supported workflow environment and is selected by
default. The test profile uses isolated 200-cell fixtures from `test_data/`;
RDS fixtures still use the normal conversion process. `make test-data` creates
only missing fixtures and never changes `input_data/`; use
`TEST_DATA_OVERWRITE=true` to refresh them. Tests default to
`aida25,cima,terekhova26`; pass `STUDIES=<list>` (including `STUDIES=all`) to
change that selection. Production `make run` still defaults to all studies.

## Choosing a Make target

Choose the smallest target that directly exercises the changed contract. The
repository instructions in [AGENTS.md](AGENTS.md) use this table when deciding
which validation to run.

| Change or goal | Preferred target | When to use a broader target |
|---|---|---|
| Python library, CLI, or deterministic test change | `make lint test-unit` | Run a focused workflow too if a task contract or published artifact changed. |
| Nextflow workflow or configuration change | `make workflow-lint` | Add `make run-test STUDIES=<study>` when task wiring must execute. |
| One single cell analysis report | `make run-cell-type-test CELL_TYPE=<slug>` (requires an existing split; create test splits with `make split-cell-types-test`) | Use `make render-cell-type-test CELL_TYPE=<slug>` to rerender one report from existing analysis artifacts, or `make render-cell-type-reports-test` to rerender all test reports. |
| Rerender single cell analysis reports from existing analysis artifacts | `make render-cell-type-reports[-test]` | Use `make render-cell-type[-test] CELL_TYPE=<slug>` to refresh one report. |
| Rerender per-type pseudobulk DE reports from existing model results | `make render-de-reports[-test]` | Runs one render task per cell type. Requires current per-type model outputs and `abundance_input.tsv`; rendering also executes abundance models, residual clustering, and WGCNA. It does not refit PyDESeq2. Rerun `make run-de[-test]` to regenerate missing abundance inputs. |
| Rebuild the central report index | `make render-report-index` | It is also refreshed automatically after core, integration-benchmark, DE, cell-type, and age-clock workflow targets. |
| Cell-type splitting or shared analysis | `make run-test`, then `make run-cell-type-analysis-test` | Use `make run-all-test` when the complete ordered test workflow is needed. |
| Positive pseudobulk DE fitting | `make run-de-synthetic-test` and `make check-synthetic-de-results` | This fixture is independent of the core merge and renders per-type DE plus cross-cell-type reports. To check single cell analysis reports too, run `make run-test`, `make run-cell-type-analysis-test`, and `make check-synthetic-de-reports`. |
| Cell-type transcriptomic age clocks | `make run-age-clock` | Fits sample-weighted, study-balanced, and within-sample rank elastic-net clocks for one type selected with `CELL_TYPE` (default `cd14-monocyte`; label or slug). CV candidates run with `AGE_CLOCK_N_JOBS=20` by default; lower it to limit CPU use. Use `make run-age-clocks` to fit every cell type present in the pseudobulk. It logs and skips types with no eligible samples in any configured training or external study; an unexpected type-specific fit failure is recorded with its traceback in `analysis/per_type/age_clock_all_types_summary.json`, successful types are still reported, and the target then exits nonzero. Defaults to holding out `nehar_belaid26,perez22`; these relatively small studies retain as much training data as possible while spanning a broad validation age range. Override with any comma-separated list via `AGE_CLOCK_EXTERNAL_STUDIES=study_a,study_b`. Per-method reports and a merged comparison go to `analysis/per_type/<slug>/ageing_clock/`. Smoke-test with `make run-age-clock-synthetic-test` and check that expected output files exist with `make check-age-clock-synthetic-test`; both are included in `make run-all-test` and `make verify`. |
| Focused gene-by-cell-type profile | `make profile-gene CELL_TYPE=<label-or-slug> GENE=<gene-id>` | Use `make run-gene-profile-synthetic-test` to fit and render a planted synthetic gene and check the report end to end. |
| Top merged age-DE gene profiles | `make profile-top-age-genes-test` | Uses the synthetic DE fixture and profiles two genes per cell type; included in `make run-all-test`. |
| DE from an existing test pseudobulk merge | `make run-de-test` | Requires `make run-test` first; the small core fixture may not have enough independent samples for successful fits. Use the synthetic target above for positive fitting coverage. |
| Focused profiles of top merged age-DE genes | `make profile-top-age-genes` | Runs ten profiles per cell type by adjusted p-value; set `TOP_AGE_GENE_COUNT` to change the number. Requires completed merged DE results and pseudobulk input. |
| Core workflow for one or a few studies | `make run-test STUDIES=<list>` | Omit `STUDIES` for the representative AIDA25/CIMA/Terekhova26 set; use `STUDIES=all` for every available fixture. |
| Global integration/annotation diagnostics | `make run-integration-benchmark[-test]` | Run after an existing single-cell merge; `make run-all-test` also exercises its test form. |
| Ordered Docker integration test | `make run-all-test` | Use `make verify` when local lint, unit, and documentation checks are also needed. |
| Markdown documentation or Make target names | `make docs-check` | Also run the target described by changed commands. |
| Broad cross-cutting or release validation | `make verify` | This is deliberately non-resumed and uses the representative test set by default; pass `STUDIES=all` for every available fixture. |
| Production analysis | `make run` | It is resumable; use only after the relevant focused checks pass. |

Nextflow submission alone is not a passing workflow check: confirm each task's
`.exitcode` and the published deliverables. Development test targets use the
representative AIDA25/CIMA/Terekhova26 fixtures and skip unavailable inputs;
pass an explicit comma-separated `STUDIES` list to change a run.

Harmony integration uses the pinned Harmony2 C++ backend. It uses the backend's
native thread policy; the legacy OpenMP/OpenBLAS thread overrides that crashed
on large production runs have been removed.

Single cell analysis reports provide local embeddings, marker contrasts, study
coverage, and local-cluster composition views. Separate per-type DE reports and a
combined cross-cell-type analysis report are generated by `run-de`; those DE
reports need only the merged pseudobulk input.

`make run-test` and `make run` create the single-cell merge by default. To skip
the production merge when it is not needed, set `MERGE_SINGLE_CELL=false`:

```bash
make run STUDIES=all MERGE_SINGLE_CELL=false
```

To benchmark global embedding and type-calling methods without rerunning that
merge:

```bash
make run-integration-benchmark
```

Study-QC, merge-QC, integration-benchmark, cell-type, DE, and gene-profile
report templates are Jupytext Python sources under `reports/`. Their runners
materialize a temporary notebook immediately before execution; generated
templates are not tracked. The published `executed.ipynb` is the audit artifact.

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

Cell-type, merge-QC, integration-benchmark, DE, combined trajectory, and gene-profile reports share
their reader-facing header and styles from
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

The report runners call these helpers after notebook execution and HTML export.

## Single cell analysis reports

The optional downstream workflow opens the completed merged single-cell H5AD in
read-only backed mode and splits it by the retained per-study
`aifi_l2_majority` call. It materializes and writes one type at a time, keeping
the full sparse expression matrix on disk during splitting. AnnData still loads
the observation metadata into memory, and the largest individual type must fit
in memory while its split is written. Experimental benchmark labels are never
used for this split. For each type, `ANALYSE_CELL_TYPE`
creates the type-specific PCA/UMAP/Leiden embedding, marker tables, PC--age
correlation table, and `analysis.h5ad`. The PCA section shows each principal
component's individual share of variance explained, not cumulative variance.
Analysis writes the raw-count `analysis.h5ad` before normalization, then updates
its metadata and embeddings on disk after analysis finishes. This avoids
holding a second full counts matrix in memory. Shared feature selection and PCA
preparation live in `pbmc_pipeline.embedding.prepare_sparse_pca`; per-type graph
construction, Harmony correction, clustering, and marker ranking are in
`pbmc_pipeline.cell_type.analyse_cell_type`.
`RENDER_CELL_TYPE_REPORT` then loads
that analysis artifact plus the primitive raw-count H5AD to create the executed
notebook and HTML. It also plots the retained per-study `aifi_l3_majority`
labels; no merged L3 prediction is created. Global neighbours are deliberately
not copied into the split artifacts.
During report execution, the notebook also computes an unintegrated UMAP from
the retained native PC scores for comparison; this embedding is not written
back to `analysis.h5ad`.

Sample-level abundance and fraction views are generated with the pseudobulk
DE reports, where cell counts are summed from the merged pseudobulk object.
The cell-type report retains local-cluster composition views and excludes
samples below `cell_type_analysis.min_cluster_fraction_denominator_cells`.
The configured `cell_type_analysis.l2_parent_l1` map supplies the parent label
in `report.json` and the cell-type manifest, and defines the parent denominator
for DE abundance. The split and analysis H5ADs do not carry a derived parent
column or per-sample abundance totals.

```bash
make run-cell-type-analysis
```

This downstream workflow requires an existing merged single-cell H5AD at
`<outdir>/harmonization/merged/single_cell_merged.h5ad`. Override that location with
`MERGED_INPUT=/path/to/single_cell_merged.h5ad`.

For a complete production run, `make run-all` executes the core workflow with
single-cell merging enabled, then runs differential expression, detailed
profiles for the ten lowest-adjusted-p-value age-DE genes in each cell type,
and cell-type analysis/reporting in order. The profiles are written under
`<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/`. Run
`make profile-top-age-genes TOP_AGE_GENE_COUNT=5` to change the number or to
rerun profiling from existing merged DE results and pseudobulk data. The
workflow is resumable. It uses the configured `STUDIES`, `OUTDIR`, and
`WORK_DIR` values. For test output, run `make run-test`
first, then `make run-cell-type-analysis-test` to consume the test merge under
`output/test/`. All test targets publish only beneath `output/test/`.

The central HTML report index is written to `<outdir>/index.html`. It links
to the merged and per-study harmonization QC reports, integration benchmark,
combined cross-cell-type analysis, and currently published residual, DE, and
age-clock reports for each type. Age-clock entries include the merged method
comparison and each method's report, based on the clock manifest. Core,
integration-benchmark, DE, cell-type, and age-clock Make targets refresh it
after publishing. The index checks report files independently, and adds known
types from the cell-type and age-clock manifests when available. To rebuild it
manually after moving or removing reports, run
`make render-report-index OUTDIR=<outdir>`.

The raw-count splits are published under `<outdir>/analysis/per_type/<cell-type-slug>/split.h5ad` so one
report can be rerun without re-splitting the merged H5AD:

```bash
make run-cell-type CELL_TYPE=cd14-monocyte
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type-reports
make render-cell-type-reports-test
```

Published outputs and the Nextflow work cache are managed directly in the
filesystem. Remove only the directories you intend to discard; production and
test artifacts live under separate output roots, while Nextflow's cache is
stored separately under the configured work directory.

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
report execution and overwrite the executed notebooks and HTML without
rerunning local PCA, Harmony integration, Leiden clustering, or marker ranking.
The report execution does calculate its comparison UMAP from saved native PCs.

Use `make check-cell-type-test-prerequisites CELL_TYPE=<slug>` before a
targeted test rerun when unsure whether its input exists. It prints the minimal
next command: create splits from the existing test merge, or create a fresh core
test merge and reports when that merge is absent. The production equivalent is
`make check-cell-type-prerequisites CELL_TYPE=<slug>`.

Each published type directory under `output/analysis/per_type/` contains its
raw-count `split.h5ad`, `analysis/analysis.h5ad`, HTML report, executed notebook,
and report tables. `cell_type_manifest.json` records the exact merged input, configuration,
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
structure. Each report now shows one study-coloured UMAP from native PCs before
the existing Harmony-integrated UMAP plots. V(D)J genes matching
`processing.exclude_vdj_regex` are excluded from the local HVG/PCA input; the
report records the count excluded. HVGs use Scanpy's configured mean-dispersion
selection on normalized, log-transformed expression.

Harmonization and per-cell-type analysis use the same PCA preparation function,
`pbmc_pipeline.embedding.prepare_sparse_pca`; the selected HVG thresholds and
maximum PC count are set in `config/pipeline.json` under `processing.hvg` and
`processing.neighbors`. It selects mean-dispersion HVGs using the configured
minimum mean and dispersion cutoffs, without a maximum mean cutoff, and excludes
genes matching `processing.exclude_vdj_regex`
(case-insensitively), falls back to all non-V(D)J genes if fewer than two HVGs
remain, casts the selected matrix to float32, scales without centering to retain
sparsity, and runs ARPACK PCA with implicit centering. The PCA component count
is capped at `processing.neighbors.n_pcs` and at the dimensions available in the
selected matrix, so analysis does not compute PCs that the configured neighbor
graph will not use. Stage-specific embedding orchestration is in
`pbmc_pipeline.harmonize._compute_embedding` and
`pbmc_pipeline.cell_type.analyse_cell_type`; Harmony correction is shared via
`pbmc_pipeline.harmonize._harmony_integrate`.

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

### Covariate policy across analyses

The general rule is to use the broadest adjustment set that the analysis
population can support. For a covariate-focused estimate, use the largest
study set with usable observations and enough variation to estimate the target
contrast. Missing or unclear values are excluded from that model; they are not
treated as another biological level. Each model records its design, included
studies, and sample counts so the actual support can be checked.

| Analysis | Covariate use and sample selection |
|---|---|
| Per-study DE | Uses age, sex, library depth, and study site when estimable, plus each optional covariate (BMI, CMV, smoking) that varies in that study. Samples missing a term in that model are excluded. This is the maximal varying optional-covariate set within each study. |
| Combined primary age DE | Uses the configured core terms: study site when estimable, age, sex, and library depth. It does not add every optional covariate, which would unnecessarily restrict the shared age-analysis cohort when those fields are sparsely recorded. |
| Combined covariate-focused DE | Fits one optional covariate at a time on top of the core terms. It uses complete cases and retains studies with usable observations; smoking is explicitly `yes vs no` with `no` as the reference. Unlike abundance models, DE has no 95% coverage rule. The current fit checks variation and estimability across the combined design, but does not impose a minimum sample count in each study or exposure group; inspect the recorded study counts before interpreting a focused result. |
| Developmental versus adult age effects | Fits the developmental model on ages [0,20) and compares its gene-level age coefficients with the adult model (ages ≥20). The report shows developmental age support by study site and sex before the effect scatterplots; the two ranges meet at age 20 without overlap. |
| Age-trajectory DE | Replaces continuous age with age bins and retains the other core terms. It uses the configured age and bin-support filters; optional covariates are not added. |
| Abundance changes | Fits adult-only, study-specific OLS models. Within each study, an optional covariate is included only if at least 95% of adult samples have a valid value and it varies. The model then uses complete cases across its selected terms. Smoking accepts only `yes` and `no`; other values are missing. These are study-specific estimates, not a pooled abundance effect. |
| Gene profiles | The PyDESeq2 age model uses the combined primary design. The exploratory age-bin model uses the same core adjustments; the per-study-slope and mixed models adjust for sex and library depth. These profile models do not add BMI, CMV, or smoking. |
| Expression atlas and residual WGCNA | Atlas technology contrasts are descriptive, not covariate-adjusted DE. WGCNA's metadata associations are post-fit associations with module eigengenes, not covariate adjustment in the gene-level fit. |

For a future covariate-focused analysis with a within-study support rule, define
the minimum sample and group counts explicitly before fitting. The current
combined DE implementation uses complete cases and global design
estimability, so it can retain a study even when that study does not contain
both levels of a categorical contrast.

The separate `differential_expression.nf` workflow consumes the merged
sample-level pseudobulk H5AD. It first lists retained per-study
`aifi_l2_majority` labels, then runs one independent task per label. Each task
loads the same small merged pseudobulk H5AD and subsets its label in memory; it
does not materialize a per-type pseudobulk H5AD. Per-study models adjust for
`study_site` when multiple sites remain in that study, along with age, sex, and
available optional covariates. The merged age-and-sex model uses
`~ study_site + age + sex + log10_total_counts` across studies. `log10(total_counts)`
from each pseudobulk is included as a technical covariate in all DE models and
its coefficient is shown in the combined-model result tables. When
present and variable, BMI, CMV, and smoking status each receive separate
combined fits adjusted for study site, age, sex, and log10 total counts.
Before PyDESeq2 fitting, continuous design columns are centered within each
fit. This changes the intercept parameterization, while preserving the
interpretation of reported slopes and contrasts and improving optimizer stability.
For each merged age contrast, the report also runs preranked enrichment against
the downloaded GO Biological Process 2023 collection. It ranks all tested genes by
their signed Wald statistic when complete; otherwise it uses signed
`−log10(p-value)`. Positive normalized enrichment scores indicate terms among
genes increasing with age and negative scores terms among genes decreasing with
age. The collapsible summary uses the native GSEA FDR q-value and shows each
term's leading-edge gene count plus the first twelve ranked leading-edge genes;
the complete untruncated list remains in the CSV.
DE fits have no minimum covariate-coverage percentage: an optional term enters
a per-study design when it varies, and each fit uses complete cases for its
own design. Smoking status is accepted only as `yes` or `no`; other labels are
treated as missing. Its categorical contrast is `yes vs no`, with `no` as the
reference.
DESeq2 estimates size factors to normalize library depth, and this additional
depth term adjusts for the age-associated pseudobulk-depth pattern observed in
Naive CD8 T cells, in addition to size-factor normalization. It is included in
the per-study and combined DE models.
Each per-type DE task writes
`analysis/per_type/<type>/differential_expression/abundance_input.tsv` and its
report includes sample-level abundance plots and adult-only fraction models.
The report task executes these models and its trajectory/residual clustering
and WGCNA analyses when rendering; rerendering does not refit the upstream
PyDESeq2 models.

Pseudobulk `n_cells` values are summed within `study × sample` for the
retained-PBMC denominator and within `study × sample × configured L1 parent`
for the parent denominator. The abundance input includes every observed sample
for each reported type, with zero cells for sample/type combinations absent
from the pseudobulk table. Each fraction view applies its own denominator
cutoff (currently 10 cells); the per-type DE pseudobulk cell filter does not
exclude zero-count abundance observations. The configured L2-to-L1 map remains
in `cell_type_analysis.l2_parent_l1`; it is applied to retained L2 labels, not
the independent `aifi_l1_majority` classifier output. The existing OLS,
bootstrap, and plot implementation moved from the cell-type report section.
Denominator and bootstrap settings are configured under
`differential_expression.abundance`; diagnostics are published alongside each
per-type DE report.

The report fits adult-only (age ≥20), study-specific OLS models. A candidate
covariate is included only when at least 95% of that study's adult samples have
valid values and it varies; the fit then uses complete cases across all
included covariates. Smoking status is accepted only as `yes` or `no`; other
labels are treated as missing, and its contrast is `yes vs no` with `no` as
the reference. This 95% rule applies to the abundance models, not the DE fits.
The site term is dropped when only one site remains. A final collector writes
the result directories and manifest.
Per-type DE reports show model status and result tables, with up to 25 rows per
fit sorted by adjusted p-value. The gene-universe table above is the canonical
definition of which genes enter each fit. The parallel type-fitting task is
configured for one CPU, which is passed directly to PyDESeq2; all other CPU,
time, and memory requests use the executor defaults.
In addition to the continuous-age fit, each cell type gets an exploratory,
site- and covariate-adjusted age-trajectory model. Per-type DE reports show
model status and result tables for available per-study and combined fits.
`run-de` also renders the combined cross-cell-type analysis report. It
summarizes shared age trajectories, cross-type recurrence, and Pearson-residual
structure, and includes a descriptive, study-balanced expression atlas with
gene availability, count-loss diagnostics, expression maps, and technical-
label contrasts. It also reports unique sample counts per technical batch
within each study site to show batch-size imbalance and batch coverage. These
are exploratory outputs, not covariate-adjusted
estimates of absolute RNA abundance. Per-type DE results are stored under each type's
`analysis/per_type/<cell-type-slug>/differential_expression/`; combined DE and
cross-type artifacts are under `analysis/combined/`. See
[IMPLEMENTATION.md](IMPLEMENTATION.md#age-trajectory-settings) for the model,
support rules, clustering, output tables, and current configurable settings.

The single cell analysis reports are generated separately by
`run-cell-type-analysis` and contain local embeddings, markers, coverage, and
fraction views. They do not consume the DE manifest. The residual WGCNA helper
is run by each per-type DE report from the fitted age-bin residual matrix. The
report includes soft-threshold and sample-connectivity diagnostics, eigengene
structure and covariate associations, hub genes, CMV evidence, and
Module GSEA ranks all selected network
genes by correlation with each module eigengene and tests GO Biological Process
2023 terms from the downloaded `input_data/go/` library. It publishes complete results as
`module_gsea.csv`, a hierarchy-aware specific-term summary, and a provenance
manifest recording the library and ontology files actually used. The summary uses the
native GSEA FDR q-value within each module and collapses a GO parent only when
a same-direction descendant has at least 0.5 leading-edge Jaccard overlap.
The strongest CMV- and smoking-associated modules are
also plotted on the residual UMAP and by observed yes/no/not-provided status.
Binary covariate heatmap values are
signed point-biserial correlations; multi-level categorical values are η², the
fraction of eigengene variance explained by category membership. These are
descriptive post-fit associations, not adjusted effects. CMV and smoking are
checked because the residual model adjusts for age bins, sex, study site, and
depth but omits these two plausible sources of biological residual structure.
The checks ask whether omitted traits remain associated with modules; they do
not assume that such effects are present. WGCNA color names are module IDs,
not biological classes. Grey is the special unassigned-gene label, so it is
excluded from module GSEA because those genes do not form a coherent module.
Module colors are assigned independently when a network is fit, so color names
are not stable across fits.
The tables are published beside the report under `residual_wgcna/`. To keep large networks readable, module-indexed
figures show at most the 25 largest modules and state when that limit applies;
the complete module-level tables remain available in that directory.

The current enrichment scope is GO Biological Process only. Molecular Function
and Cellular Component are not tested or combined into the reported FDR family.
The reusable runner accepts a caller-provided gene-set collection when a future
analysis needs one, and records that collection's content hash in its output
provenance.

Samples must be age 20 or older and have at least 10 cells in that
sample × cell-type pseudobulk. Samples missing required metadata or covariates,
or with zero counts over the fit's gene universe, are also omitted. Optional
covariate missingness receives additional handling within each study: a
per-study model includes an optional covariate only when it varies in that
study, and samples missing a covariate are omitted from fits that include it.
Combined optional-covariate models use complete cases from studies that recorded
the covariate. All sample exclusions and model fits that cannot estimate the
configured covariates are recorded in the JSON manifest. Combined results and
the manifest are published beneath
`<outdir>/analysis/combined/differential_expression/`; per-cell-type results
and task fit records are under
`<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/`.

The DE manifest's `results_by_cell_type` index lists the per-study and combined
result paths for each covariate (or an empty result set) for each cell-type
slug. `run-de` uses it to render separate per-type DE reports. The cell-type
workflow does not read this manifest: `make render-cell-type-reports` rerenders
single cell analysis reports from published splits and analysis artifacts
without repeating analysis. That renderer writes single-cell artifacts beneath
`<outdir>/analysis/per_type/<cell-type-slug>/analysis/`; it never publishes
report files into the DE directory. Named Nextflow subworkflows group splitting
plus analysis, analysis of an existing split, and report rendering while
retaining the separate single-cell and pseudobulk entry points.

### Focused gene profiling

`make profile-gene CELL_TYPE=<label-or-slug> GENE=<gene-id>` fits the configured
merged DESeq2 model and a naive normalization linear model for one cell type
from the pseudobulk H5AD, then renders a focused report for one gene inside
the Python Docker image. For example:

```bash
make profile-gene CELL_TYPE=cd14-monocyte GENE=IL7R
```

The executed notebook and HTML report give each model its own top-level
heading. The PyDESeq2 section includes the gene's position in the fitted
mean–dispersion plot, its age coefficient and adjusted p-value, sample-level
size-factor normalized counts, and Pearson residuals against age and
`log10(total_counts)`, each faceted by study and colored by study site. The
mean–dispersion plot marks the selected gene's gene-wise dispersion before
shrinkage alongside the final dispersion used in the model. Correlated-gene
scatter plots color points by study. The
residual section shows a table with the selected gene followed by its ten
strongest residual correlates, including each gene's mean normalized count and
final dispersion, and plots the four strongest correlates. A separate Naive
normalization section documents geometric-mean total-count scaling and the
`log(normalized count + 1)` transform once, shows the normalized-count
histogram, and compares naive expression with PyDESeq2 VST. Continuous-age
models use a geometric mean from their included samples; the age-bin model
calculates its geometric mean across all samples in that cell type before
applying model filters, so lower-age samples contribute to its scale. The naive
linear age model uses the PyDESeq2 design covariates and reports its coefficient
table and residual plots against age and library depth. The age-bin model uses
the configured decade width, includes samples below age 20, adjusts for the
other varying design terms, and plots standardized fitted levels alongside the
sample expression values, with residual plots against age and library depth.

A further model allows a separate linear age slope and intercept per study,
adjusting for sex and library depth. `age * study` estimates study-specific
intercept contrasts and age-by-study slope differences, without pooling; its
residuals are also plotted against age and library depth. The mixed model keeps
study intercepts fixed and estimates random age slopes around a shared slope.
It reports fixed effects, random-slope and residual standard deviations, fit
convergence, and partially pooled conditional slopes for each study, plus
residual plots against age and library depth. Treat slope-variation estimates
as exploratory when only a few studies are available. Sample-level expression
points use 0.3 opacity so overlapping samples remain visible. The displayed
geometric means use scientific notation. The HTML table
of contents includes each implemented model heading, deferred `limma-voom` and
R DESeq2 headings marked `TODO: implement R execution`, and rendered
subheadings. The
target also writes a compressed, full fit-gene-set
VST matrix and the included sample metadata under
`<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/`. VST values are
transformed data for downstream exploration and modeling; differential
expression inference should continue to use raw integer counts.

The focused profile refits the combined primary age DESeq2 model from the
pseudobulk input, then caches that fit under
`<outdir>/analysis/per_gene/<cell-type-slug>/.model_fit_cache/` for other gene
profiles of the same cell type. It does not currently reuse a fitted model
from `run-de`; that workflow publishes coefficient tables and fit metadata,
not a reusable fitted PyDESeq2 object. A shared fit artifact would be cleaner
for the matching combined primary design, provided it records and matches the
exact included samples, gene universe, design, and software version. Per-study
fits and covariate-focused fits have different samples or designs and must not
be substituted for that age model.

When completed `run-de` results are available, each focused profile also shows
the selected gene's existing age-bin trajectories across cell types. These are
the adjusted log2 fold changes relative to each cell type's own reference age
bin, so they describe within-type trajectories and do not compare expression
abundance between cell types. The report publishes the displayed rows and a
per-cell-type availability table beside the profile.

To smoke-test the complete fit and report with planted synthetic expression,
run `make run-gene-profile-synthetic-test`. It writes the pseudobulk fixture to
`output/test/analysis/combined/gene_profile_synthetic/` and the checked report artifacts to
`output/test/analysis/per_gene/cd14-monocyte/synth-age-marker-1/`. This target is
included in `make run-all-test`; it is separate from pytest because it fits
DESeq2 and executes a notebook.

`make profile-top-age-genes` applies the same detailed report to the top age-DE
genes from each cell type's merged model results. Ranking uses finite adjusted
p-values in `<outdir>/analysis/per_type/<cell-type>/differential_expression/merged.csv`; it
profiles up to ten genes per type by default, using the merged pseudobulk H5AD.
The report paths follow the same `analysis/per_gene/<cell-type>/<gene>/` layout as
the single-gene target.

The Nextflow `test` profile used by `make run-de-test` enables a test-only mode
that bypasses only the minimum-cell cutoff. It keeps the adult-age, metadata, and model-estimability
checks, so the small fixture can exercise PyDESeq2 without changing the
production inclusion rules. DE test-mode runs are labeled `test_only`; their
DE metadata, result CSVs, and DE reports carry a warning that the results are
for software validation and must not be treated as biological evidence. Use
`make run-de` for the configured analysis rules.

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

### Synthetic age-clock experiments

The synthetic pseudobulk generator is deliberately compact: its single
implementation, [`src/pbmc_pipeline/synthetic_de.py`](src/pbmc_pipeline/synthetic_de.py),
defines the studies, assay-availability flags, planted effects, and count
sampling in `create_synthetic_pseudobulk`. For one-off clock diagnostics,
modify that source directly rather than adding a scenario framework; record
the command and seed, and revert temporary source edits afterwards.

The default fixture is a software smoke test, not evidence that a clock
transfers between real sites. In particular, when `synthetic_study_c` is the
external study, development contains only A and B. Their study-only genes are
unavailable in the opposite leave-one-study-out fold and deliberately have
opposite age slopes. This is a useful assay-coverage stress case, but can
produce large opposite fold biases and should not be interpreted as clock
performance. To isolate that mismatch, temporarily set the two
`study_specific_age_effects` slopes to zero. The sample count is allocated
equally over nine age bins, so `--samples-per-study` must be divisible by nine;
for example, 360 creates the tenfold diagnostic fixture.

Keep exploratory inputs and reports in a dedicated ignored workspace directory
such as `output/experiments/<experiment-name>/`, not only under system-wide
`/tmp`: files outside the repository may be inaccessible to browsers or the
shared workspace UI. Inspect gene-intersection loss, missing-gene coverage,
fold bias, and top coefficients before interpreting feature selection. Any
retained fixture change must preserve the default smoke-test contract, update
`adata.uns["synthetic_test_data"]`, and have focused unit-test coverage.

After `make run-de-synthetic-test`, run `make profile-top-age-genes-test` to
render two genes per cell type, ranked by merged age-DE adjusted
p-value. The outputs are written under
`output/test/analysis/per_gene/<cell-type>/<gene>/`; this target is also
included in `make run-all-test` and `make verify`.

`make run-all-test` is the complete ordered Docker test. It runs the core and
integration benchmark on available downsampled real-data fixtures, then fits
the separate synthetic pseudobulk DE fixture and checks its results. The
synthetic age-clock smoke test creates its own synthetic pseudobulk and checks
that model, prediction, CV, report, and index files exist, then the suite
profiles its top two merged age-DE genes per cell type, runs the focused
synthetic gene-profile smoke test, renders per-type DE and combined
cross-cell-type reports from that synthetic fixture, and renders single cell
analysis reports from the real merge/splits.
`make verify` runs it after linting, unit
tests, Nextflow lint, and documentation checks. `STUDIES=<list>` changes only
the core test inputs; use `STUDIES=all` for all available fixtures. The synthetic WGCNA fixture exercises the ordinary
network and report path, but its
two-module result does not exercise the 25-module figure cap or production-scale
runtime.

The generated input stays under `output/test/analysis/combined/synthetic_de/`.
DE outputs and the combined manifest are written under
`output/test/analysis/combined/differential_expression/`; `run-de-synthetic-test`
renders per-type DE reports under `output/test/analysis/per_type/<slug>/differential_expression/`.
The separate single cell analysis reports are rendered from the test single-cell merge by
`run-cell-type-analysis-test` under each type's `analysis/` directory.

`run-cell-type-analysis` renders the single cell analysis report and requires a
merged single-cell input. `run-de` uses the pseudobulk merge to render the
separate per-type DE reports and cross-cell-type analysis report, so those
reports do not require a single-cell merge.

Run the workflow after creating the pseudobulk merge:

```bash
make run-de
```

`make run-de` renders one per-type DE report task for each completed fit,
alongside the separate combined cross-cell-type analysis report. If a report
template changes, Nextflow resume reuses completed DE fits and reruns the
affected downstream report tasks.

By default this reads `output/harmonization/merged/pseudobulk_merged.h5ad`. Override the
input or output roots with `PSEUDOBULK_INPUT=/path/to/pseudobulk_merged.h5ad`
or `OUTDIR=/path/to/output`.
`make run-de-test` requires an existing test pseudobulk merge; create it first
with `make run-test`. It uses the explicitly labeled test mode described above.
Use `make run-de-synthetic-test` to exercise successful PyDESeq2 fits with the
production sample-inclusion filters enabled.

Each study has one maximal DE model: `study_site` when it varies, age, sex,
log10 total counts, and every optional covariate with at least two observed
values. All coefficients for that study therefore come from one formula and
complete-case sample set.
A missing BMI, CMV, or smoking-status value
excludes that pseudobulk only from models containing that term; one missing
value removes one row. Combined models target one covariate at a time, using
all complete-case samples from studies that recorded it and adjusting for
study site, age, sex, and log10 total counts. The site term is omitted when one
site remains. CMV negative/positive values are harmonized to no/yes. Smoking
status is modeled only as yes/no, with no as the reference group; other labels
are treated as missing. Categorical result tables name their
comparison groups, including male versus female, yes CMV versus no CMV, and
yes versus no smoking status.

## Further reading

See [INPUT_FILES.md](INPUT_FILES.md) for exact external paths and checksums and
[IMPLEMENTATION.md](IMPLEMENTATION.md) for architecture and output contracts.

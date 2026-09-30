# PBMC ageing pipeline plan

## Goal of the current milestone

Produce one independently reproducible, validated H5AD file per study. Each file
must contain raw counts, homogeneous observation metadata, and AIFI L1/L2/L3
cell-type labels. The workflow also builds study/sample/AIFI-L2 pseudobulks and
a merged pseudobulk matrix; an opt-in whole-dataset single-cell
merge is available.

## Decisions

- `config/studies.json` is the source of truth for study-specific inputs,
  preparation-adapter dependencies, assumptions, and processing choices.
- `config/harmonized_obs_schema.json` defines the output metadata contract.
- Installed `pbmc-*` package commands are the only processing entry points;
  Make invokes Nextflow for complete workflows.
- Raw counts remain in `X` in final outputs. Normalized values are temporary.
- Source-study QC is authoritative. The pipeline does not apply additional
  manual cell filtering or QC thresholds, except for an adapter's explicitly
  declared cohort-exclusion rule when a source object combines studies.
- `subject` identifies a biological individual; `sample` identifies one
  specimen at one collection timepoint; and `batch_single_cell` identifies a
  technical processing unit. Pseudobulks aggregate all technical partitions
  of one sample. AIDA Lonza material is intentionally retained as separate
  site-specific samples until it is excluded downstream.
- Legacy notebooks are retained as provenance but are not pipeline dependencies.
- Exploratory notebooks are kept locally but excluded from the initial Git
  history because their embedded outputs dominate file size and contain
  machine-specific paths. A source-only archival snapshot can be added later.
- Inferred values and unresolved questions are recorded explicitly in config.
- Test data and full data use the same code path; only the input root changes.
- Use Nextflow as the workflow layer. Its process-level container support cleanly
  separates the RDS conversion environment from Python harmonization and gives
  resumable, parallel per-study execution. Snakemake would also work, but offers
  less benefit for this deliberately mixed-language/containerized portfolio.
- Keep the `Makefile` thin: it provides memorable lint, test, build, and focused
  run commands, while Nextflow remains the only production dependency graph.
- All processing commands accept explicit input/output paths. Nextflow owns
  input staging, per-study parallelism, caching, resources, and output publishing.
- Full and test Wang25 runs enter through their respective configured RDS
  inputs, so smoke tests exercise the RDS conversion process.
- Nextflow's `work/` is the canonical cache for workflow intermediates. The
  existing `cache/converted` remains useful only for direct/manual conversion.
- The Python image installs the reviewed runtime/report `requirements.lock`
  before copying `src/`, then installs the package with
  `--no-deps --no-build-isolation`. `requirements.dev.lock` extends the same
  base with test and lint tools for the local `.venv`; source-only image rebuilds
  reuse the runtime dependency layer.
- SciPy is pinned to 1.16.0 to prevent a segmentation fault in `sc.tl.umap`.

## Milestones

- [x] Inventory the repository and legacy notebooks.
- [x] Define the initial study configuration and homogeneous metadata schema.
- [x] Add a parameterized harmonization command and validation/reporting modules.
- [x] Add dependency metadata and configuration-level tests.
- [x] Install the scientific Python environment and run full smoke tests on all
  configured downsampled H5AD files.
  - All configured studies complete metadata harmonization, embeddings, CellTypist
    where required, H5AD writing, report generation, and round-trip loading.
- [x] Add the AIDA and Terekhova auxiliary metadata and validate their
  preparation adapters.
- [x] Recover and document the Wang metadata join and RDS conversion that
  produced the legacy intermediate `2_wang.h5ad`.
- [ ] Review all `needs_review` provenance entries against publications/source
  metadata.
- [x] Run the full production pipeline successfully on production datasets,
  including the real Wang25 RDS conversion (confirmed by the user).
- [ ] Review the production per-study QC reports and confirm their scientific
  outputs are suitable for downstream analysis.
- [x] Add legacy-compatible per-study sample × AIFI-L2 pseudobulks,
  outer-joined pseudobulk merge with gene-availability flags, a combined merge
  report, and a gene-presence UpSet plot.
- [x] Add an opt-in whole-dataset single-cell merge with a
  shared-gene embedding and experimental merged AIFI-L2 diagnostics.
- [ ] Review the full production merge artifacts and merged QC reports.
- [ ] Add one shared cross-study design report, rather than duplicating
  donor-level covariate/confounding plots in every AIFI-L2 report. It should
  cover covariate support and correlation across studies, plus high-level
  AIFI-L1 composition and age-stability summaries.
- [x] Add the initial Python harmonization Docker image definition.
- [x] Add an executable QC notebook and self-contained HTML report generation.
- [x] Replace notebook-based RDS conversion with a one-input/one-output R script.
- [x] Add and test the R/Seurat/reticulate/anndata conversion image.
  - A synthetic Seurat v5 object round-tripped to H5AD with raw counts, feature
    and cell names, and observation metadata intact.
  - `sceasy` was removed from the runtime path because version 0.0.7 uses the
    defunct Seurat 5 `slot=` API. The converter now constructs AnnData directly.
- [x] Build the Python harmonization/QC image and execute a fresh containerized
  Wang25 test workflow through both harmonization and QC.
- [x] Add a lint-clean DSL2 Nextflow workflow for conversion, harmonization, QC,
  per-study selection, test inputs, Docker execution, and resumable caching.
- [x] Run the containerized Nextflow test profile for all configured studies,
  including preparation, pseudobulk/single-cell merges, and QC (29 successful
  processes).
- [x] Run the real Wang25 RDS conversion through the successful production
  workflow.
- [ ] Compare its production H5AD against the prior converted representation.
- [x] Remove the redundant Docker project mount in `nextflow.config`; Nextflow
  now stages the project inputs without a duplicate Docker mount.
- [x] Enforce the study and input-source registry structures with tracked JSON
  Schemas, in addition to semantic cross-reference checks.
- [x] Publish retained source-`obs` CSV sidecars keyed by `cell_id` for every
  prepared study, while keeping the harmonized H5AD metadata contract compact.
- [x] Split per-cell-type computation from notebook/HTML rendering so analysis
  artifacts are explicit workflow outputs.
- [x] Pin the complete Python Docker runtime in `requirements.lock` and make it
  part of the Python image's content-derived tag.
- [x] Add the separate pseudobulk PyDESeq2 workflow with per-study and merged
  age models, configured sample filters, task-level cell-type fan-out without
  materialized per-type pseudobulk inputs, and a manifest indexed by cell type.
- [x] Allow standard cell-type reports to include matching DE results through
  that manifest, keeping DE files and single-cell reports in separate folders.
- [x] Generate a deterministic positive-fit DE fixture on demand and include
  fit, planted-marker, and report-artifact assertions in `make verify`.

## Local development state

- The scientific dependencies are installed in the ignored `.venv`.
- All configured studies pass the complete test-data workflow.
- Ruff passes, 52 Python tests pass, and `nextflow lint` reports no errors for
  the core, cell-type, and differential-expression workflows (Nextflow 26.04.6;
  minimum declared version is 25.04).
- Revision-tagged Python and R images build successfully.
- Harmony integration uses the pinned `harmonypy==2.0.2` C++ backend through a
  direct adapter that validates the cells-by-PC output orientation and records
  the implementation version in embedding provenance. Legacy OpenMP/OpenBLAS
  overrides were removed after the backend upgrade.
- Cell-type reports now use signed marker contrasts, direct UMAP cluster labels,
  a compact study-coverage representation, PC1/2 and PC3/4 score plots coloured
  by study and shaped by technology, a companion intronic-read/3′--5′ PCA row,
  up to ten PCA loading panels, and a cumulative variance-explained plot.
  Sample-fraction views include within-study linear trends with and without
  samples under age 20, separate combined linear and LOWESS plots with unshaded
  per-study fits and sample-share-weighted means, plus adult-only adjusted
  fraction-model forest plots and published residual/sampling diagnostics.
  Local PCA age correlations are limited to the first ten PCs. A focused
  CD14-monocyte run completed both the separate analysis and report-rendering
  tasks.
- A fresh `docker,test` Wang25 Nextflow run published its H5AD, JSON run report,
  executed QC notebook, and HTML report under an ignored temporary output
  directory.
- The AIFI Docker test run exercised harmonization and QC; timestamped
  application progress logs are captured in each Nextflow task's `.command.err`.
- Test fixtures live under ignored `test_data/`; conversion intermediates live
  in Nextflow's ignored `work/` directory. `tests/` contains test code only.
- A narrowed `make verify TEST_STUDIES=wang25` passed, including the generated
  DE fixture, planted age-marker checks, and standard report artifact checks.
- The full production pipeline has since completed successfully on production
  datasets, including real Wang25 RDS conversion. Production QC and merge
  artifacts still need explicit review and their run records should be retained.

## Expected commands

```bash
make install
make lint test-unit
make workflow-lint
make docs-check
make run-test STUDIES=wang25
make run-de-synthetic-test
make verify TEST_STUDIES=wang25
make run STUDIES=all MERGE_SINGLE_CELL=true
```

## Known blockers and scientific questions

- AIDA, Terekhova, and Wang source metadata are local ignored inputs in their
  respective `input_data/<study_id>/` directories. The 424 MB Terekhova
  table has been losslessly reduced for pipeline purposes to a 9.2 MB compressed
  cell-to-tube lookup; its original should be archived externally.
- OneK1K's sample-per-donor interpretation and intronic-read setting need review.
- AIDA's technology version and ambiguous Lonza batch assignments need review.
- Terekhova freezing and demultiplexing descriptions need review.
- Wang's freezing and intronic-read interpretations need review. Its published
  CMV IgM field is retained as a qualitative serostatus and is negative for all
  61 workbook records.
- Decide whether final study files should retain UMAP/PCA artifacts or only labels.

## Planned shared cross-study design and AIFI-L1 report

Create a single reader-facing report from the merged inputs, published beside
the merge-level QC rather than repeated for every cell type. It should use one
row per retained study × sample (or biological individual where that distinction
is material) and make the study design visible before downstream interpretation.

- Summarize availability, within-study variation, and pairwise association of
  age, sex, BMI, CMV, technology, intronic-read inclusion, and other model
  covariates. Display missingness and study-level confounding explicitly; do
  not imply that a cross-study association identifies an independent effect.
- Add a high-level AIFI-L1 composition section. For each L1 compartment, show
  study-specific mean sample fractions alongside the underlying sample values,
  so the potentially large between-study differences in mean composition are
  immediately visible. Keep biological replication at the sample level rather
  than cell-weighting the comparison.
- For each AIFI-L1 compartment, summarize age-dependent variation using the
  same study-aware sample-level framing: show within-study trends and their
  uncertainty or a clearly labelled descriptive age model. The working
  expectation is that these broad L1 fractions are comparatively stable with
  age, but the report must show the observed effect sizes and precision rather
  than encode that expectation as a conclusion.
- Link this global design/AIFI-L1 report from per-cell-type HTML reports when
  it exists. It is context for interpretation, not a replacement for the
  adjusted per-type fraction models or pseudobulk DE workflow.

## Primary analysis: age-associated pseudobulk DE

The first implementation uses **PyDESeq2** on sample × cell-type raw-count
pseudobulks. It runs one model per study (`~ age + sex`) and one shared-slope
merged model (`~ study + age + sex`) for each AIFI L2 type. The age coefficient
is a linear effect per year, tested with a Wald test. Nextflow lists cell types
then runs one task per type; each task loads and subsets the same small merged
pseudobulk H5AD in memory, without materializing a per-type pseudobulk input.
The final collector writes the standard result layout and root manifest.
Per-study fits use genes available in that study; the merged fit uses the
intersection of genes available in all studies contributing samples to that
fit, avoiding synthetic outer-join zeros. Each parallel type task is configured
for one CPU; other CPU, time, and memory requests are left to executor defaults.

The DE workflow remains independent of per-celltype analysis; celltype reports
consume matching DE result directories through the root manifest's explicit
cell-type result index. The single-cell and pseudobulk workflows stay separate;
named Nextflow subworkflows group split/analyse/render stages, and one report
renderer accepts an optional matching DE result. Their output products remain
separate under `cell_type_analysis/` and `differential_expression/`.

The test profile bypasses only the minimum-cell cutoff so the small real-data
fixture can exercise filtering and graceful skip behavior while retaining the
adult and model-estimability safeguards. A separate deterministic synthetic
pseudobulk fixture provides the positive-fit path: it generates independent
sample rows with adequate age/sex variation and at least 15 cells per
pseudobulk, then applies the ordinary production filters. Synthetic results
carry a test-only warning throughout the CSVs, manifest, and reports. Do not
upsample real test samples to imitate replication; more sampled cells in one
sample do not increase the number of independent observations. Neither kind of
test-mode result is biological evidence.

The first pass includes samples aged at least 20 years with at least 10 cells
in that sample × cell-type pseudobulk. It also requires complete sex metadata;
all sample exclusions are recorded. Per-cell-type reports optionally show the
per-study recurrence of FDR-significant genes and the merged volcano plot when
those DE results are present. BMI and CMV are not model covariates yet. Later
analysis will probably add all available covariates, including BMI and CMV,
to per-study plots. Keep nonlinear age and an age-20 sensitivity analysis as
future decisions. R DESeq2, edgeR quasi-likelihood, and limma-voom comparisons
are outside this first implementation.

## Planned exploratory analysis: residual structure after pseudobulk DE

After fitting the pseudobulk differential-expression model separately for each
cell type, retain a sample × gene residual matrix: expression not explained by
the prespecified model covariates. This is an exploratory companion to the
primary age-DE analysis, not a replacement for it.

- Prespecify the model formula and explicitly decide whether age is removed
  before generating residuals. If age is in the model, residual structure for
  age-DE genes represents heterogeneity beyond their fitted age trend; if it is
  not, age-associated covariance remains in the residual matrix. Keep these
  interpretations separate.
- Inspect sample--sample correlation, clustering, and low-dimensional views of
  the residual matrix to identify reproducible participant-level patterns or
  possible ageing phenotypes. Test whether apparent groups replicate across
  studies and are not driven by library size, batch, study, or sparse samples.
- For primary age-DE genes, inspect gene--gene residual correlations and
  cluster/module structure. This can reveal co-varying programmes whose mixed
  contributions may underlie a single marginal age-DE list, and lets each
  participant be described by a combination of module scores rather than a
  single ageing label.
- Report module scores and their associations with available metadata only as
  exploratory results, with appropriate multiple-testing and cross-study
  validation. Do not infer discrete ageing subtypes from an unreplicated
  clustering.
- As a targeted biological example, score a documented CMV-response signature
  in memory T-cell pseudobulks and test whether it aligns with residual modules
  or participant clusters. Establish the signature source, direction, and
  scoring method in advance; distinguish measured CMV status from inferred
  signature activity.

## Production follow-up and architecture risks

1. **Memory remains deployment-specific.** The full production run succeeded,
   but this repository does not yet retain a task-level RAM/CPU profile for it.
   Harmonization and QC load a full study object; single-cell merge streams the
   concat to disk, then loads the merged matrix into RAM for normalization and
   integration. Record production peak usage before changing resource requests
   or sizing a different runner.
2. **Keep study transformation code explicit.** `studies.json` now selects a
   preparation adapter and declares only its auxiliary inputs. Each adapter
   materializes exactly one canonical metadata row per retained expression cell
   before the generic harmonizer runs. A subset is permitted only when the
   adapter explicitly declares a cohort-exclusion rule in configuration. The
   registry documents are now Schema-validated; keep their versioning and the
   semantic cross-reference checks in step as new study features are added.
3. **Preserve source metadata for auditability.** Preparation now publishes a
   retained `<study>.source_obs.csv.gz` sidecar keyed by `cell_id`, containing
   source `obs` rows only for retained cells. It is an audit artifact rather
   than part of the harmonized metadata contract. Revisit a columnar format if
   sidecars become too large for convenient inspection.
4. **Freeze remaining dependencies for production.** The complete Python image
   environment is now pinned in `requirements.lock`, and the lock participates
   in the content-derived image tag. The remaining reproducibility work is to
   publish immutable image digests for production releases and, where practical,
   record package artifact hashes in addition to version pins.
5. **Complete run provenance.** `run_manifest.json` records the Git revision,
   image references, configuration checksum, selected studies, and hashes of
   artifacts and reports. Add the input and CellTypist model checksums, random
   seeds, and production image digests so the successful run can be reproduced
   and audited.
6. **Atomic writes are low priority inside Nextflow.** Failed processes remain in
   isolated work directories and outputs are published only after success. A
   temporary-write/validate/rename helper remains useful for large H5AD files
   written directly via the CLI, but no larger transaction system is warranted.
7. **Test edge cases explicitly.** Random 200-cell subsets may omit rare labels,
   batches, or join cases. Keep them for end-to-end smoke tests and add small
   synthetic fixtures for metadata and conversion edge cases.
8. **Add optional resource profiling.** Start with Nextflow's built-in execution
   report, task trace, and timeline on a production run; summarize per-process
   peak RSS, CPU use, and elapsed time with run metadata. The merge already logs
   phase timings. Add in-process memory sampling only if task-level metrics do
   not explain its peaks. Keep these measurements observational rather than
   making noisy resource estimates CI pass/fail criteria.

## Recommended next sequence

1. Review the successful production QC, merge, and DE outputs; archive the run
   manifest, task logs, and execution reports with the production artifacts.
2. Resolve the scientific `needs_review` entries against source publications
   and metadata, preserving any remaining inferences explicitly.
3. Complete the run manifest with input/model checksums and random seeds; record
   immutable Python and R image digests for the production run.
4. Add opt-in Nextflow execution reports and traces, then record per-task RAM,
   CPU, and elapsed-time baselines from production.
5. Add targeted synthetic tests for metadata joins, rare labels, and conversion
   edge cases that random smoke-test subsets may omit.
6. Continue with the exploratory residual-structure analysis only after the
   production outputs and provenance have been reviewed.
7. Polish the existing static HTML reports for reader-facing use (without
   sacrificing the executed notebooks as audit artifacts), then add a static
   manifest-driven report index before considering a custom web frontend.

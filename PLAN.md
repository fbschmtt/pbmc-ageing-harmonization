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
- [ ] Run the full datasets and review per-study QC reports.
- [x] Add legacy-compatible per-study sample × AIFI-L2 pseudobulks,
  outer-joined pseudobulk merge with gene-availability flags, a combined merge
  report, and a gene-presence UpSet plot.
- [x] Add an opt-in whole-dataset single-cell merge with a
  shared-gene embedding and experimental merged AIFI-L2 diagnostics.
- [ ] Run and review merge artifacts on the full frozen studies.
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
- [ ] Run the real Wang25 RDS conversion and compare its H5AD against the
  existing converted/test representation.
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

## Local development state

- The scientific dependencies are installed in the ignored `.venv`.
- All configured studies pass the complete test-data workflow.
- Ruff passes, 46 Python tests pass, and both `nextflow lint main.nf` and
  `nextflow lint cell_type_analysis.nf` report no errors (verified with
  Nextflow 26.04.6; minimum declared version is 25.04).
- Revision-tagged Python and R images build successfully.
- Harmony integration uses the pinned `harmonypy==2.0.2` C++ backend through a
  direct adapter that validates the cells-by-PC output orientation and records
  the implementation version in embedding provenance. Legacy OpenMP/OpenBLAS
  overrides were removed after the backend upgrade.
- Cell-type reports now use signed marker contrasts, direct UMAP cluster labels,
  a compact study-coverage representation, PC1/2 and PC3/4 score plots coloured
  by study and shaped by technology, up to ten PCA loading panels, and a
  cumulative variance-explained plot. Sample-fraction views include within-study
  linear trends with and without samples under age 20, plus separate combined
  linear and LOWESS plots with unshaded per-study fits and sample-share-weighted
  means. Local PCA age correlations are limited to the first ten PCs. A focused
  CD14-monocyte run completed both the separate analysis and report-rendering
  tasks.
- A fresh `docker,test` Wang25 Nextflow run published its H5AD, JSON run report,
  executed QC notebook, and HTML report under an ignored temporary output
  directory.
- The AIFI Docker test run exercised harmonization and QC; timestamped
  application progress logs are captured in each Nextflow task's `.command.err`.
- Test fixtures live under ignored `test_data/`; conversion intermediates live
  in Nextflow's ignored `work/` directory. `tests/` contains test code only.

## Expected commands

```bash
make install
make lint test-unit
make workflow-lint
make docs-check
make run-test STUDIES=wang25
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

## Architecture risks to address before full-data runs

1. **Memory remains the principal full-data constraint.** Harmonization reads a
   complete H5AD, copies counts, and performs selected-HVG scaling/PCA. Full
   Terekhova (~1.9 million cells) needs a measured memory budget before
   execution. QC also loads the complete object.
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
5. **Improve provenance.** `run_manifest.json` now records the Git revision,
   image references, configuration checksum, selected studies, and hashes of
   merged artifacts. Add input/model checksums and random seeds next.
6. **Atomic writes are low priority inside Nextflow.** Failed processes remain in
   isolated work directories and outputs are published only after success. A
   temporary-write/validate/rename helper remains useful for large H5AD files
   written directly via the CLI, but no larger transaction system is warranted.
7. **Test edge cases explicitly.** Random 200-cell subsets may omit rare labels,
   batches, or join cases. Keep them for end-to-end smoke tests and add small
   synthetic fixtures for metadata and conversion edge cases.

## Recommended next sequence

1. Populate verified direct input URLs and checksums; manual-only sources are
   recorded in `config/input_sources.json`.
2. Measure peak memory on the largest testable/full study before full execution.
3. Run and validate the real Wang25 conversion through Nextflow.
4. Resolve the scientific `needs_review` entries and add complete input hashes.
5. Record immutable Python and R image digests for production releases.
6. Run each full study independently before launching all studies together.

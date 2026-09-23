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
- The Python image uses a temporary minimal package to cache dependency
  installation before copying `src/`; the final install removes its setuptools
  build artifacts and uses `--no-deps`. This is safe while packaging metadata
  remains static in `pyproject.toml`.
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

## Local development state

- The scientific dependencies are installed in the ignored `.venv`.
- All configured studies pass the complete test-data workflow.
- Ruff passes, 23 Python tests pass, and `nextflow lint main.nf` reports no
  errors (verified with Nextflow 26.04.6; minimum declared version is 25.04).
- Revision-tagged Python and R images build successfully.
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
make verify
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

## Architecture risks to address before full-data runs

1. **Memory remains the principal full-data constraint.** Harmonization reads a
   complete H5AD, copies counts, and performs selected-HVG scaling/PCA. Full
   Terekhova (~1.9 million cells) needs a measured memory budget before
   execution. QC also loads the complete object.
2. **Keep study transformation code explicit.** `studies.json` now selects a
   preparation adapter and declares only its auxiliary inputs. Each adapter
   materializes exactly one canonical metadata row per retained expression cell
   before the generic harmonizer runs. A subset is permitted only when the
   adapter explicitly declares a cohort-exclusion rule in configuration. Add
   formal JSON Schema or Pydantic validation as the configuration grows.
3. **Preserve source metadata for auditability.** The harmonized output currently
   replaces `obs` with the standardized schema. Consider a namespaced source
   metadata sidecar (preferably Parquet) so mappings can be debugged later.
4. **Freeze remaining dependencies for production.** SciPy is pinned to 1.16.0
   because later resolution caused a `sc.tl.umap` segmentation fault; the R base
   image is digest-pinned and key Python packages in that bridge are pinned.
   The main Python environment otherwise still uses version ranges. Add a lock
   file and use immutable/versioned image references in real runs. Rebuilding
   mutable image tags is avoided by revision-derived local tags. A lock file
   remains desirable for the main Python image.
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
5. Lock the Python environment and record immutable image tags/digests.
6. Run each full study independently before launching all studies together.

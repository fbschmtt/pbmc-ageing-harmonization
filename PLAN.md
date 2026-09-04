# PBMC ageing pipeline plan

## Goal of the current milestone

Produce one independently reproducible, validated H5AD file per study. Each file
must contain raw counts, homogeneous observation metadata, and AIFI L1/L2/L3
cell-type labels. Merging and pseudobulk generation are deliberately deferred.

## Decisions

- `config/studies.json` is the source of truth for study-specific inputs,
  metadata mappings, assumptions, and processing choices.
- `config/harmonized_obs_schema.json` defines the output metadata contract.
- `scripts/harmonize_study.py` is the only per-study processing entry point.
- Raw counts remain in `X` in final outputs. Normalized values are temporary.
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

## Milestones

- [x] Inventory the repository and legacy notebooks.
- [x] Define the initial study configuration and homogeneous metadata schema.
- [x] Add a parameterized harmonization command and validation/reporting modules.
- [x] Add dependency metadata and configuration-level tests.
- [x] Install the scientific Python environment and run full smoke tests on all
  five downsampled H5AD files.
  - All five studies complete metadata harmonization, embeddings, CellTypist
    where required, H5AD writing, report generation, and round-trip loading.
- [x] Add the AIDA and Terekhova auxiliary metadata and validate their joins.
- [x] Recover and document the Wang metadata join and RDS conversion that
  produced the legacy intermediate `2_wang.h5ad`.
- [ ] Review all `needs_review` provenance entries against publications/source
  metadata.
- [ ] Run the full datasets and review per-study QC reports.
- [ ] Freeze the five harmonized files before designing the merge stage.
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
- [ ] Run the containerized Nextflow test profile for all five studies.
- [ ] Run the real Wang25 RDS conversion and compare its H5AD against the
  existing converted/test representation.

## Local development state

- The scientific dependencies are installed in the ignored `.venv`.
- All five studies pass the complete test-data workflow.
- Ruff passes, all seven Python tests pass, and `nextflow lint main.nf` reports no
  errors (verified with Nextflow 26.04.6; minimum declared version is 25.04).
- Both `pbmc-ageing-python:local` and `pbmc-ageing-r:local` build successfully.
- A fresh `docker,test` Wang25 Nextflow run published its H5AD, JSON run report,
  executed QC notebook, and HTML report under an ignored output directory.
- Test fixtures live under ignored `test_data/`; conversion intermediates live
  in Nextflow's ignored `work/` directory. `tests/` contains test code only.

## Expected commands

```bash
python -m pip install -e '.[dev,qc]'
python scripts/harmonize_study.py --study onek1k --test
python scripts/harmonize_study.py --study all --test --validate-only
pytest
nextflow lint main.nf
make docker-build
make pipeline-test STUDIES=wang25
```

## Known blockers and scientific questions

- AIDA, Terekhova, and Wang metadata are local ignored inputs under
  `input_data/metadata/`. The 424 MB Terekhova
  table has been losslessly reduced for pipeline purposes to a 9.2 MB compressed
  cell-to-tube lookup; its original should be archived externally.
- OneK1K's sample-per-donor interpretation and intronic-read setting need review.
- AIDA's technology version and ambiguous Lonza batch assignments need review.
- Terekhova freezing and demultiplexing descriptions need review.
- Wang's ancestry, freezing, intronic-read, and CMV interpretations need review.
- Decide whether final study files should retain UMAP/PCA artifacts or only labels.

## Architecture risks to address before full-data runs

1. **Memory is the immediate risk.** Harmonization currently reads a complete
   H5AD, copies the count matrix, and scales an HVG matrix for PCA. QC also loads
   the complete object. This is especially unsafe for Terekhova (~1.9 million
   cells). Refactor before attempting all full datasets: keep counts backed or
   sparse, avoid dense scaling of all cells, load only marker-gene columns for
   QC, and calculate summaries in chunks. Plotting every cell does not require
   loading the full expression matrix.
2. **Do not let study JSON become a transformation language.** It is appropriate
   for mappings, constants, and join declarations. New complex study-specific
   reshaping should live in a small preparation adapter that emits normalized
   metadata for the generic harmonizer. Add formal JSON Schema or Pydantic
   validation as the configuration grows.
3. **Preserve source metadata for auditability.** The harmonized output currently
   replaces `obs` with the standardized schema. Consider a namespaced source
   metadata sidecar (preferably Parquet) so mappings can be debugged later.
4. **Lock remaining dependencies.** The R base image is digest-pinned and key
   Python packages in that bridge are pinned, but the main Python environment
   still uses version ranges. Add a lock file and use immutable/versioned image
   references in real runs. Rebuilding a mutable `:local` tag does not invalidate
   Nextflow cache entries; run without `-resume` after rebuilding or change the
   image tag.
5. **Improve provenance.** Reports currently include a configuration digest and
   package versions. Add the Git commit, input/model checksums, image digest,
   command invocation, and random seeds.
6. **Atomic writes are low priority inside Nextflow.** Failed processes remain in
   isolated work directories and outputs are published only after success. A
   temporary-write/validate/rename helper remains useful for large H5AD files
   written directly via the CLI, but no larger transaction system is warranted.
7. **Test edge cases explicitly.** Random 200-cell subsets may omit rare labels,
   batches, or join cases. Keep them for end-to-end smoke tests and add small
   synthetic fixtures for metadata and conversion edge cases.

## Recommended next sequence

1. Run `make pipeline-test STUDIES=all` and inspect all QC.
2. Refactor harmonization and QC memory behavior, then measure peak memory on
   the largest testable inputs.
3. Run and validate the real Wang25 conversion through Nextflow.
4. Resolve the scientific `needs_review` entries and add complete input hashes.
5. Lock the Python environment and record immutable image tags/digests.
6. Run each full study independently before launching all studies together.

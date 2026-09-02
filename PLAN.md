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
- [ ] Recover or document the preprocessing that produced Wang's intermediate
  `2_wang.h5ad`.
- [ ] Review all `needs_review` provenance entries against publications/source
  metadata.
- [ ] Run the full datasets and review per-study QC reports.
- [ ] Freeze the five harmonized files before designing the merge stage.
- [x] Add the initial Python harmonization Docker image definition.
- [x] Replace notebook-based RDS conversion with a one-input/one-output R script.
- [ ] Add and test the R/Seurat/sceasy conversion image.
- [ ] Build and execute the Python image on all five test inputs.
- [ ] Add a Nextflow workflow after the per-study input contracts are settled.

## Environment limitations in the current session

- The workspace exposes `.git` as an empty read-only mount. It cannot be removed,
  written, or initialized, including through the approved elevated operation.
  The requested initial commit must therefore be created once that mount is
  removed by the workspace host.
- The scientific dependencies are installed in the ignored `.venv`, and all five
  studies pass test-data schema validation.

## Expected commands

```bash
python -m pip install -e '.[dev]'
python scripts/harmonize_study.py --study onek1k --test
python scripts/harmonize_study.py --study all --test --validate-only
pytest
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

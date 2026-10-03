# PBMC Ageing Harmonization plan

## Current stage

The production workflow has completed successfully on production inputs, and
a fresh-clone run completed input acquisition and `make run-all` (confirmed by
the user). Those runs predate the latest source-metadata corrections. After
the final manual metadata check, rerun production to publish the corrected
harmonized metadata and merge-report covariates.

The delivered workflow produces one harmonized H5AD per study with raw counts,
canonical observation metadata, and AIFI L1/L2/L3 labels. It also builds
study/sample/AIFI-L2 pseudobulks, a merged pseudobulk matrix, and an opt-in
whole-dataset single-cell merge.

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
- [ ] Resolve the remaining provenance questions: OneK1K sample identity and
  Wang25 fresh-status confirmation from the authors.
- [ ] Complete one final manual per-study metadata check against each paper and
  its supplements. Verify the harmonized metadata and technical covariates,
  record paper/source evidence and discrepancies, and keep unresolved values
  explicitly marked as inferred or pending confirmation.
- [x] Run the full production pipeline successfully on production datasets,
  including the real Wang25 RDS conversion (confirmed by the user).
- [x] From a fresh clone, acquire inputs and complete `make run-all`
  successfully (confirmed by the user).
- [ ] Rerun the production workflow after the final manual metadata check so
  corrected study metadata reaches the published H5ADs and merge report.
- [x] Add legacy-compatible per-study sample × AIFI-L2 pseudobulks,
  outer-joined pseudobulk merge with gene-availability flags, a combined merge
  report, and a gene-presence UpSet plot.
- [x] Add an opt-in raw-count whole-dataset single-cell merge, with global
  Harmony and AIFI-L2 diagnostics in a separate benchmark task.
- [x] Move global UMAPs and benchmark AIFI-L2 concordance matrices from the
  merge QC notebook into the integration-benchmark report.
- [x] Include the integration benchmark and report artifact check in the
  ordered test workflow after the fresh core merge.
- [ ] Extend the integration benchmark beyond Harmony and unintegrated PCA once
  comparison methods and quantitative acceptance criteria are specified.
- [ ] Add one shared cross-study design report, rather than duplicating
  donor-level covariate/confounding plots in every AIFI-L2 report. It should
  cover age, sex, BMI, and CMV support and correlation across studies, plus
  high-level AIFI-L1 composition and age-stability summaries. The merge QC
  report already gives sample-level distributions for technology,
  intronic-read handling, and frozen status.
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
  including preparation, pseudobulk/single-cell merges, and QC.
- [x] Run the real Wang25 RDS conversion through the successful production
  workflow.
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
  fit, planted-marker, and report-artifact assertions in `make run-all-test`
  and `make verify`.

## Current implementation state

- The ignored `.venv` supplies local checks; revision-tagged Python and R
  images supply workflows.
- The ordered test workflow runs a fresh real-fixture core merge, integration
  benchmark and report check, synthetic positive DE fit, and cell-type reports.
  Fixtures live under ignored `test_data/`; Nextflow intermediates live in
  ignored `work/`.
- The core single-cell merge preserves only raw counts, canonical metadata, and
  harmonization-stage QC fields. Global normalization, Harmony, graphs, UMAP,
  and benchmark labels are computed only by the optional integration workflow.
- Harmony uses the pinned `harmonypy==2.0.2` backend and records embedding
  provenance. Cell-type reports and DE results are separate workflow outputs;
  their detailed contracts are in `IMPLEMENTATION.md` and `README.md`.
- Production, including Wang25 RDS conversion, and fresh-clone input acquisition
  plus `make run-all` completed successfully. Both predate the final metadata
  corrections and require the planned production rerun.

## Expected commands

```bash
make install
make lint test-unit
make workflow-lint
make docs-check
make run-test STUDIES=wang25
make run-de-synthetic-test
make run-all-test
make verify
make run STUDIES=all MERGE_SINGLE_CELL=true
```

## Open provenance questions and follow-up work

- AIDA, Terekhova, and Wang source metadata are local ignored inputs in their
  respective `input_data/<study_id>/` directories. The 424 MB Terekhova
  table has been losslessly reduced for pipeline purposes to a 9.2 MB compressed
  cell-to-tube lookup; its original should be archived externally.
- OneK1K has no sample or collection-timepoint key in its source H5AD. The
  adapter currently uses `donor_id` as both `sample` and `subject`; confirm
  that each donor contributes one biological collection before interpreting
  the pseudobulk sample count.
- Wang25 is provisionally marked fresh based on the paper's selection of fresh
  reference samples; author confirmation is pending.
- Include OneK1K and Wang25 intronic-read handling, AIDA25 Lonza batch
  assignments, and Terekhova23 demultiplexing in the final metadata review.
  Keep unresolved assumptions labelled as inferred.
- Wang25's published CMV IgM field is retained as qualitative serostatus and is
  negative for all 61 workbook records.
- Decide whether final study files should retain UMAP/PCA artifacts or only labels.

## Planned shared cross-study design and AIFI-L1 report

Create a single reader-facing report from the merged inputs, published beside
the merge-level QC rather than repeated for every cell type. It should use one
row per retained study × sample (or biological individual where that distinction
is material) and make the study design visible before downstream interpretation.

- Summarize availability, within-study variation, and pairwise association of
  age, sex, BMI, CMV, and other model covariates. Technology, intronic-read
  inclusion, and frozen status are covered in the merge QC report. Display
  missingness and study-level confounding explicitly; do not imply that a
  cross-study association identifies an independent effect.
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

The experimental branch `experiment/study-site-de-adjustment` uses **PyDESeq2**
on sample × cell-type raw-count pseudobulks. It runs one maximal model per
study, adding `study_site` when multiple sites remain, and one shared-slope
merged model (`~ study_site + age + sex`) for each AIFI L2 type. Combined
covariate models also adjust for `study_site`, age, and sex; the site term is
dropped when only one site remains. The age coefficient is a linear effect per
year, tested with a Wald test. Nextflow lists cell types
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
pseudobulk, plants age, sex, BMI, and CMV signals, and applies the ordinary
production filters. BMI is available in two synthetic studies and CMV in one,
with missing values to exercise complete-case selection; one synthetic study
contains multiple sites to exercise site adjustment. Synthetic results
carry a test-only warning throughout the CSVs, manifest, and reports. This is
workflow-path coverage, not biological realism: sample counts are balanced,
gene-availability differences are sparse, and effects do not vary by study.
Do not
upsample real test samples to imitate replication; more sampled cells in one
sample do not increase the number of independent observations. Neither kind of
test-mode result is biological evidence.

The first pass includes samples aged at least 20 years with at least 10 cells
in that sample × cell-type pseudobulk. It also requires complete sex metadata;
all sample exclusions are recorded. Per-cell-type reports show per-study,
per-covariate significance inside/outside the shared tested-gene intersection, and combined
volcano plots for age, sex, BMI, and CMV where metadata support the fits. BMI
and CMV use separate complete-case models adjusted for study site, age, and
sex; per-study models include site when it varies within the study. Each model
records its included studies, sites, and sample counts. Categorical coefficient labels show the
reference group (female for sex; no/negative for CMV). Keep
nonlinear age and an age-20 sensitivity analysis as future decisions. R DESeq2,
edgeR quasi-likelihood, and limma-voom comparisons are outside this first
implementation.

### Priority investigation: unusual Naive CD8 T-cell age signal

The production age volcano for Naive CD8 T cells has an unusually strong,
asymmetric pattern. Diagnose it against the matching production pseudobulk and
DE outputs before treating it as a biological result:

- Reproduce the fit and inspect its sample-level age, study-site, technology,
  sex, and pseudobulk cell-count distributions. Compare the pooled age
  coefficient with within-study and within-site fits, study/site-specific age
  slopes, and leave-one-study/site-out fits. Check count normalization,
  size-factor behavior, filtering, and the result-to-plot mapping as part of
  the reproduction.
- Assess whether conditioning on the minimum 10 cells per sample × cell-type
  pseudobulk induces selection/collider bias, especially for cell types whose
  abundance falls with age. Compare inclusion against age, site, abundance,
  and expression-related measurements, and run prespecified threshold
  sensitivity analyses before deciding whether the cutoff contributes to the
  signal.
- Evaluate the experimental `study_site` batch adjustment against the current
  study-level design before adopting it in production. This uses source labels
  for AIDA sites/countries, Fachrul26 villages, and Nehar-Belaid26 collection
  sites, while studies with one site retain one level. Check label uniqueness,
  site-level sample sizes and age support, design rank, and changes in
  estimability. In Nehar-Belaid26, the adult metadata include 32 UCHC and 9
  HSNRI samples before cell-type cell-count filtering, so site-specific support
  may be limited.

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
   concat to disk, then loads the joined raw matrix to write its artifact. The
   optional integration benchmark owns full-matrix normalization and integration.
   Record production peak usage before changing resource requests or sizing a
   different runner.
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
4. **Revisit merged cross-study cell-type fractions.** The current combined
   age-fraction view is descriptive, with a shared age slope, equal study
   weights, and unclustered sample-level uncertainty. Further work should
   assess study-specific age effects, uneven age support, and repeated subjects
   before choosing an inferential cross-study summary.
5. **Freeze remaining dependencies for production.** The complete Python image
   environment is now pinned in `requirements.lock`, and the lock participates
   in the content-derived image tag. The remaining reproducibility work is to
   publish immutable image digests for production releases and, where practical,
   record package artifact hashes in addition to version pins.
6. **Complete run provenance.** `run_manifest.json` records the Git revision,
   image references, configuration checksum, selected studies, and hashes of
   source inputs, metadata dependencies, annotation models, artifacts, and
   reports. Add random seeds and production image digests so the successful run
   can be reproduced and audited.
7. **Atomic writes are low priority inside Nextflow.** Failed processes remain in
   isolated work directories and outputs are published only after success. A
   temporary-write/validate/rename helper remains useful for large H5AD files
   written directly via the CLI, but no larger transaction system is warranted.
8. **Test edge cases explicitly.** Random 200-cell subsets may omit rare labels,
   batches, or join cases. Keep them for end-to-end smoke tests and add small
   synthetic fixtures for metadata and conversion edge cases.
9. **Add optional resource profiling.** Start with Nextflow's built-in execution
   report, task trace, and timeline on a production run; summarize per-process
   peak RSS, CPU use, and elapsed time with run metadata. The merge already logs
   phase timings. Add in-process memory sampling only if task-level metrics do
   not explain its peaks. Keep these measurements observational rather than
   making noisy resource estimates CI pass/fail criteria.
10. **Keep report computation separate from presentation.** Cell-type and
    integration reports already use workflow artifacts and tested helpers where
    practical. Extract further reusable data preparation and diagnostics when a
    second report needs them; leave Jupytext sources as reader-facing layouts.

## Recommended next sequence

1. Reproduce and diagnose the unusual Naive CD8 T-cell age result, including
   study/site-specific fits and sensitivity to the 10-cell inclusion cutoff.
2. Compare the experimental `study_site` DE models with the study-adjusted
   baseline, checking site support, design estimability, and changes in which
   cell types and covariates can be fitted.
3. Complete the final manual metadata review for every study against its paper
   and supplements, keeping unresolved assumptions explicitly labelled.
4. If the site-adjusted fits are supported by the diagnostics, rerun production
   so published outputs match the reviewed inputs and model specification.
5. Add random-seed and immutable image-digest provenance to the production run
   record; input, metadata, model, artifact, and report checksums are already
   included.
6. Add opt-in Nextflow execution reports and traces, then record per-task RAM,
   CPU, and elapsed-time baselines from production.
7. Add targeted synthetic tests for metadata joins, rare labels, and conversion
   edge cases that random smoke-test subsets may omit.
8. Continue with exploratory residual-structure analysis after its model and
   validation design is prespecified.
9. Add a static manifest-driven report index before considering a custom
   project-level web frontend.

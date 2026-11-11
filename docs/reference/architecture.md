# Architecture

[Documentation](../README.md)

## Purpose

PBMC Ageing Harmonization harmonizes PBMC ageing studies into a common raw-count H5AD
contract, creates comparable sample-level pseudobulks, and produces merge-aware
QC reports. Nextflow is the workflow engine; the study JSON is the source of
truth for source-specific metadata handling.

## Data flow

```text
source H5AD / Seurat RDS
        │
        ├─ RDS conversion (only configured studies)
        ▼
per-study preparation adapter
  exactly one canonical metadata row per retained source cell
        ▼
per-study harmonization
  raw counts + standardized obs + AIFI labels
        ├─ per-study QC
        ▼
sample × AIFI-L2 pseudobulk per study
        ▼
outer-gene merged pseudobulk ─┐
                              ├─ one combined merge QC report
optional: all harmonized cells → raw joined single-cell H5AD
                                      ├─ integration benchmark → method labels/UMAPs, graph metrics, report
                                      └─ cell-type split → per-type analysis and reports

The core merge writes no global normalization, embedding, graph, or prediction.
```

## Orchestration

`pipeline.nf` composes the named scientific workflows under `workflows/`.
Standalone entry scripts use the same components. Make provides commands and
operational setup, and its shared runner records task failures and refreshes the
workspace report index once. See [workflow orchestration](../development/orchestration.md)
for dependency ownership, failure isolation, and test composition.

## Components

- `config/studies.json`: input locations, feature repair, preparation-adapter
  identity, auxiliary input dependencies, and annotation settings. Its
  structure is enforced by `config/studies.schema.json` before workflow work
  is scheduled.
- `config/input_sources.json`: the manifest of expression and supplementary
  inputs. Its schema and the semantic configuration checks verify that every
  configured expression input and adapter dependency appears in this manifest.
- `src/pbmc_pipeline/studies/`: explicit per-study adapters. They own unusual
  source joins and reshaping, then write one identity-checked metadata row for
  every retained input cell. An adapter may retain a strict cohort subset only
  when its configuration explicitly permits it.
- `src/pbmc_pipeline/preparation.py`: validates exact expression-cell coverage,
  or an explicitly permitted source-cell subset, and writes the portable
  prepared-cell artifact. The canonical schema includes `disease_status`;
  current non-Perez cohorts are declared healthy, while Perez22 normalizes its
  embedded disease field and excludes SLE cells before this artifact is written.
  It also writes a retained-source-`obs` CSV sidecar keyed by `cell_id`, keeping
  source metadata available for adapter audit without expanding the harmonized
  H5AD contract.
- `config/pipeline.json`: common processing parameters, model locations, and
  merge policy. Pseudobulk grouping is `sample` and `aifi_l2_majority`, matching
  the legacy notebook; study-qualified row IDs prevent cross-study collisions.
- `src/pbmc_pipeline/harmonize.py`: creates the validated per-study H5AD.
  Its run report records feature-label provenance metrics, including the counts
  of symbol-like labels, Ensembl-like identifiers, missing labels, and duplicate
  labels before repair.
- `src/pbmc_pipeline/embedding.py`: the shared `prepare_sparse_pca` function
  owns normalized mean-dispersion HVG selection, V(D)J exclusion and fallback,
  sparse scaling, and implicitly centered ARPACK PCA for harmonization and
  per-cell-type analysis.
- `src/pbmc_pipeline/merge.py`: creates raw-count pseudobulks, outer-joins
  pseudobulk genes, and runs the optional single-cell merge.
- `modules/*.nf`: conversion, harmonization, merge, QC, and manifest process
  modules; `main.nf` is limited to study selection and DAG composition.
- `cell_type_analysis.nf`: a separate downstream entry workflow. It receives an
  existing merged single-cell H5AD, opens it in read-only backed mode, and
  materializes one retained per-study AIFI-L2 type at a time for splitting.
  It then fans out separate per-type analysis and report-rendering tasks.
- `modules/cell_type_analysis.nf`: the split, analysis, and report-rendering processes.
- `integration_benchmark.nf` and `modules/integration_benchmark.nf`: the
  optional global-integration workflow. It runs each configured method in a
  separate sequential task, then aggregates successful outputs and renders a
  report from the existing raw merged single-cell H5AD.
- `src/pbmc_pipeline/integration_benchmark.py` and
  `reports/integration_benchmark_report.py`: benchmark fitting and its
  reader-facing notebook source.
- `src/pbmc_pipeline/age_clock_cross_type.py`: combines same-method
  cell-type clock predictions for donors shared across types and reports
  prediction metrics and pairwise residual correlations.
- `reports/cell_type_report.py`: Jupytext source for the per-type notebook, rendering visualizations from the completed type-specific analysis artifact.
- `reports/differential_expression_report.py`: standalone Jupytext source for
  the per-type pseudobulk DE report, consuming prepared diagnostic artifacts.
- `src/pbmc_pipeline/de_diagnostics.py` and `abundance.py`: prepare abundance
  models, trajectory and residual analyses, WGCNA and enrichment in a separate
  task, with checksummed inputs and artifacts recorded in `diagnostics.json`.
- `src/pbmc_pipeline/notebook.py`: shared Jupytext execution and HTML export,
  with isolated caches, notebook configuration, titles, logging and TOCs.
- `src/pbmc_pipeline/cell_type_assets.py`: publishes versioned named images and
  site metadata alongside the single-cell report. The exporter reads these
  artifacts without interpreting notebook headings or code.
- `differential_expression.nf` and `modules/differential_expression.nf`: the
  independent pseudobulk DE workflow: label listing, one in-memory type subset
  and fit per task, parallel per-type report rendering from each fit result,
  result/manifest collection, and a separate combined report.
- `src/pbmc_pipeline/differential_expression.py`: configured sample filtering,
  per-study and merged PyDESeq2 fits, task-result collection, CSV output, and
  the result manifest. It delegates age-bin model preparation and trajectory
  statistics to `src/pbmc_pipeline/age_trajectory.py`.
- `src/pbmc_pipeline/config.py`: configuration validation and the typed,
  required `AgeTrajectorySettings` view used by the DE fit and report stages.
- `src/pbmc_pipeline/trajectory_analysis.py` and
  `src/pbmc_pipeline/trajectory_report.py`: support-filtered per-type and
  shared-bin cross-type hierarchical clustering of gene trajectories,
  cross-type residual PCA/UMAP and sample clustering, recurrence summaries,
  and data for the combined cross-cell-type report.
- `reports/trajectory_analysis_report.py`: source notebook for the combined
  cross-cell-type analysis report rendered by `run-de`.
- `scripts/render_report_index.py`: builds the lightweight HTML index from
  harmonization QC, integration benchmark, the combined cross-cell-type report,
  and available per-type residual, DE, and age-clock reports. When present, the
  run and cell-type manifests restrict the study and per-type lists.
- `src/pbmc_pipeline/synthetic_de.py`: deterministic on-demand positive-fit
  pseudobulk fixture generation (36 samples per study, balanced across two
  developmental and seven adult age decades); the generated H5AD is ignored
  output, not a checked-in fixture.
- `scripts/profile_gene.py`, `src/pbmc_pipeline/gene_profile_report.py`, and
  `reports/gene_profile_report.py`: standalone focused gene profiling from a
  merged pseudobulk, refitting PyDESeq2 and then fitting naive total-count
  normalization models: a covariate-adjusted linear age model, categorical
  age-bin model, unpooled study-specific intercept and slope model, and a
  fixed-study-intercept/random-slope mixed model. Each model has sample-level
  residual plots against age and library depth, faceted by study. The age-bin
  normalization reference uses all cell-type samples before model filters.
  When completed DE results are supplied, the report also reads the existing
  age-bin fits for the selected gene across cell types, preserving a support
  table for types where that comparison is unavailable.
  `scripts/profile_top_age_genes.py` defines the ranking;
  `workflows/gene_profiles.nf` selects genes from explicit DE outputs, prepares one
  immutable model per type, and dispatches calculation and rendering tasks.
- `scripts/check_gene_profile.py`: output-contract and notebook smoke checks
  for the focused gene report. `make run-gene-profile-synthetic-test` fits and
  renders a planted synthetic gene profile; `make profile-top-age-genes-test`
  renders two top merged age-DE genes per cell type from the synthetic DE
  fixture. The composed `make run-all-test` shares its synthetic fixture across these checks.
- `scripts/check_synthetic_de.py`: output-contract checks for positive
  synthetic DE and trajectory reports. They cover planted age, sex, BMI, CMV,
  and smoking-status effects; covariate-specific study selection and
  complete-case counts; trajectory clustering exports; and report content.
  The fixture uses 36 independent samples per study, balanced over ages 0–90;
  the seven adult trajectory decades meet the production support threshold and
  two developmental decades support the `[0,20)` comparison.
- `src/pbmc_pipeline/reporting.py`: testable report data transformations.
- `reports/merge_qc_report.py`: Jupytext source for the generated merge
  notebook. Report runners materialize these Python sources in a temporary
  location immediately before execution; generated templates are not tracked.
- `reports/qc_report.py`: the Jupytext source for the study-QC notebook.
  The runner materializes it in a temporary location before execution.

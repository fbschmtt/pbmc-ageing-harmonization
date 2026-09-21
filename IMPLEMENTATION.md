# Pipeline architecture

## Purpose

This pipeline harmonizes PBMC ageing studies into a common raw-count H5AD
contract, creates comparable sample-level pseudobulks, and produces merge-aware
QC reports. Nextflow is the workflow engine; the study JSON is the source of
truth for source-specific metadata handling.

## Data flow

```text
source H5AD / Seurat RDS
        │
        ├─ RDS conversion (only configured studies)
        ▼
per-study harmonization
  raw counts + standardized obs + AIFI labels
        ├─ per-study QC
        ▼
sample × AIFI-L2 pseudobulk per study
        ▼
outer-gene merged pseudobulk ─┐
                              ├─ one combined merge QC report
optional: all harmonized cells → shared-gene embedding → AIFI-L2 ─┘

optional: all harmonized cells → shared-gene embedding → AIFI-L2 → merge QC
```

## Components

- `config/studies.json`: input locations, feature repair, joins, metadata
  mapping, and per-study annotation settings.
- `config/pipeline.json`: common processing parameters, model locations, and
merge policy. Pseudobulk grouping is `sample` and `aifi_l2_majority`, matching
  the legacy notebook; study-qualified row IDs prevent cross-study collisions.
- `src/pbmc_pipeline/harmonize.py`: creates the validated per-study H5AD.
- `src/pbmc_pipeline/merge.py`: creates raw-count pseudobulks, outer-joins
  pseudobulk genes, and runs the optional single-cell merge.
- `modules/*.nf`: conversion, harmonization, merge, QC, and manifest process
  modules; `main.nf` is limited to study selection and DAG composition.
- `src/pbmc_pipeline/reporting.py`: testable report data transformations.
- `reports/merge_qc_report.py`: Jupytext source for the generated merge
  notebook; use `make notebook-sync` after editing it.
- `reports/qc_report.ipynb` and `reports/merge_qc_report.ipynb`: executable
  templates rendered to self-contained HTML reports.

## Outputs

For `--outdir <outdir>`, the normal outputs are:

- `<outdir>/harmonized/<study>.h5ad`
- `<outdir>/pseudobulk/<study>.pseudobulk.h5ad`
- `<outdir>/merged/pseudobulk_merged.h5ad`
- `<outdir>/reports/*.json`
- `<outdir>/qc/<study>/report.html`
- `<outdir>/qc/merged/report.html` (one report for both merge branches)
- `<outdir>/run_manifest.json`

Merged pseudobulk `.var` records per-study gene availability. A gene absent
from an input study is represented by an outer-join zero only for that study;
use `available_in_<study>`, `n_studies_with_gene`,
`synthetic_zero_filled_studies`, and `has_synthetic_zeros` to distinguish those
synthetic zeros from observed zero counts.

The optional `--merge_single_cell` branch additionally writes
`merged/single_cell_merged.h5ad`. It is disabled by default because it loads
all selected studies simultaneously. It preserves per-study predictions in
`obs['aifi_l2_study_majority']` before calculating the merged prediction in
`obs['aifi_l2_majority']`.
Its embedding uses only the gene intersection across studies: it selects HVGs
from that intersection, scales only those selected genes, runs PCA, applies
Harmony over `obs['study']`, and builds the neighbor graph/UMAP from
`X_pca_harmony`. The backward-compatible `aifi_l2_majority` prediction uses
that Harmony-derived neighbor graph; a second majority-voting result based on
the unintegrated `X_pca` graph is retained in
`aifi_l2_unintegrated_majority`. The saved graph and UMAP remain Harmony-based.

## Running and verification

Use Make targets rather than raw Nextflow or Docker commands:

```bash
make help
make images
make run-test STUDIES=all
make verify
make run STUDIES=all OUTDIR=/path/to/output WORK_DIR=/path/to/work
make run STUDIES=all MERGE_SINGLE_CELL=true
```

`make verify` runs linting, unit tests, Nextflow lint, cached image builds, and
a non-resumed all-study Docker test workflow. `make run` is the resumable
production entry point; `make run-no-qc` omits only report rendering.
`MERGE_SINGLE_CELL=true` enables the optional single-cell branch.

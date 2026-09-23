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
optional: all harmonized cells → shared-gene embedding → AIFI-L2 ─┘

separate downstream workflow:
merged single-cell H5AD → one split task → one notebook/report task per AIFI-L2 type
```

## Components

- `config/studies.json`: input locations, feature repair, preparation-adapter
  identity, auxiliary input dependencies, and annotation settings.
- `config/input_sources.json`: the manifest of expression and supplementary
  inputs. Configuration loading verifies that every configured expression input
  and adapter dependency appears in this manifest.
- `src/pbmc_pipeline/studies/`: explicit per-study adapters. They own unusual
  source joins and reshaping, then write one identity-checked metadata row for
  every retained input cell. An adapter may retain a strict cohort subset only
  when its configuration explicitly permits it.
- `src/pbmc_pipeline/preparation.py`: validates exact expression-cell coverage,
  or an explicitly permitted source-cell subset, and writes the portable
  prepared-cell artifact. The canonical schema includes `disease_status`;
  current non-Perez cohorts are declared healthy, while Perez22 normalizes its
  embedded disease field and excludes SLE cells before this artifact is written.
- `config/pipeline.json`: common processing parameters, model locations, and
merge policy. Pseudobulk grouping is `sample` and `aifi_l2_majority`, matching
  the legacy notebook; study-qualified row IDs prevent cross-study collisions.
- `src/pbmc_pipeline/harmonize.py`: creates the validated per-study H5AD.
  Its run report records feature-label provenance metrics, including the counts
  of symbol-like labels, Ensembl-like identifiers, missing labels, and duplicate
  labels before repair.
- `src/pbmc_pipeline/merge.py`: creates raw-count pseudobulks, outer-joins
  pseudobulk genes, and runs the optional single-cell merge.
- `modules/*.nf`: conversion, harmonization, merge, QC, and manifest process
  modules; `main.nf` is limited to study selection and DAG composition.
- `cell_type_analysis.nf`: a separate downstream entry workflow. It receives an
  existing merged single-cell H5AD, loads it once to split by retained per-study AIFI-L2,
  then fans out one executed report notebook per type.
- `modules/cell_type_analysis.nf`: the split and per-type report processes.
- `reports/cell_type_report.py`: Jupytext source for the per-type notebook. It
  plots the inherited global UMAP once, then creates the type-specific analysis.
- `src/pbmc_pipeline/reporting.py`: testable report data transformations.
- `reports/merge_qc_report.py`: Jupytext source for the generated merge
  notebook. Report runners materialize these Python sources in a temporary
  location immediately before execution; generated templates are not tracked.
- `reports/qc_report.ipynb`: the separately maintained study-QC template
  rendered to self-contained HTML reports.

## Outputs

For `--outdir <outdir>`, the normal outputs are:

- `<outdir>/harmonized/<study>.h5ad`
- `<outdir>/prepared/<study>.cells.csv.gz`
- `<outdir>/prepared/<study>.prepare.json`
- `<outdir>/pseudobulk/<study>.pseudobulk.h5ad`
- `<outdir>/merged/pseudobulk_merged.h5ad`
- `<outdir>/reports/*.json`
- `<outdir>/qc/<study>/report.html`
- `<outdir>/qc/merged/report.html` (one report for both merge branches)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/analysis.h5ad` (optional derived embeddings and clusters)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/report.html` (optional downstream report)
- `<outdir>/cell_type_analysis/cell_type_manifest.json` (optional downstream provenance and completeness contract)
- `<outdir>/run_manifest.json`

### Harmonized expression contract

Each harmonized H5AD has exactly one expression representation: `.X` contains
the configured unnormalized integer count matrix. `.raw` is `None` and
`.layers` is empty. Normalized matrices may be created transiently for
annotation or embeddings, but are never retained in the output. This prevents
an alternate normalized matrix from being mistaken for analysis input.

Merged pseudobulk `.var` records per-study gene availability. A gene absent
from an input study is represented by an outer-join zero only for that study;
use `available_in_<study>`, `n_studies_with_gene`,
`synthetic_zero_filled_studies`, and `has_synthetic_zeros` to distinguish those
synthetic zeros from observed zero counts.

The optional `--merge_single_cell` branch additionally writes
`merged/single_cell_merged.h5ad`. It is disabled by default because it loads
all selected studies simultaneously. It preserves per-study predictions in
`obs['aifi_l2_majority']` for downstream grouping and splitting, with a
redundant `obs['aifi_l2_study_majority']` comparison alias.
Its embedding uses only the gene intersection across studies: it selects HVGs
from that intersection, scales only those selected genes, runs PCA, applies
Harmony over `obs['study']`, and builds the neighbor graph/UMAP from
`X_pca_harmony`. Experimental majority-voting results from the Harmony-derived
and unintegrated graphs are retained in `experimental_aifi_l2_majority` and
`experimental_aifi_l2_unintegrated_majority`; neither replaces per-study L2.
The saved graph and UMAP remain Harmony-based.

### Cell-type analysis contract

`make run-cell-types` consumes the existing
`<outdir>/merged/single_cell_merged.h5ad`; it does not trigger the core merge.
The internal split task retains raw counts, canonical metadata, and only the inherited
global UMAP. It computes `n_cells_in_sample` and
`n_cells_in_sample_l1_parent` before splitting, so the reports need no auxiliary
inputs for their sample-level composition plots. The report notebook computes the
type-specific PCA/neighbours/UMAP/clusters in its first non-rendering cell. It
excludes configured V(D)J genes from local feature selection, uses study-Harmony
PCs for the local graph/UMAP/clusters, retains native PCs for the descriptive
PC--age table, and plots the inherited per-study L3 labels without creating a
merged L3 prediction. Its published directory contains the resulting H5AD,
report, executed notebook, and tables; the primitive split object remains only
in Nextflow's cache.

The 50-cell local-analysis threshold does not guarantee enough cells per study
for reliable Harmony correction. The fixed 50-neighbour graph is also dense for
rare types. Local UMAP geometry and clusters must therefore be checked against
study colouring and composition, especially for small or imbalanced types.

The downstream manifest records checksums for its input, configuration,
reference models, and each published report. It is written only after validating
that every report emitted by the split fan-out has its H5AD, HTML, executed
notebook, and recognized completion status. `analysis_version` pins the
configuration-defined downstream analysis contract.

The configured `l2_parent_l1` mapping is an inferred taxonomy with recorded
provenance, rather than the independent `aifi_l1_majority` predictions. It is
applied only to retained per-study L2 labels; the experimental merged L2 calls
are not used in this hierarchy. Validate the mapping against the AIFI atlas
before using within-L1 fractions for inference.

## Running and verification

Use Make targets rather than raw Nextflow or Docker commands:

```bash
make help
make images
make run-test STUDIES=all
make verify
make run STUDIES=all OUTDIR=/path/to/output WORK_DIR=/path/to/work
make run STUDIES=all MERGE_SINGLE_CELL=true
make run-cell-types
make run-cell-types-test
```

`make verify` runs linting, unit tests, Nextflow lint, cached image builds, and
non-resumed core and downstream all-study Docker test workflows. `make run` is the resumable
production entry point; `make run-no-qc` omits only report rendering.
`MERGE_SINGLE_CELL=true` enables the optional single-cell branch.

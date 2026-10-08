# PBMC Ageing Harmonization architecture

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
                                      ├─ integration benchmark → labels, UMAP, report
                                      └─ cell-type split → per-type analysis and reports

The core merge writes no global normalization, embedding, graph, or prediction.
```

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
- `src/pbmc_pipeline/merge.py`: creates raw-count pseudobulks, outer-joins
  pseudobulk genes, and runs the optional single-cell merge.
- `modules/*.nf`: conversion, harmonization, merge, QC, and manifest process
  modules; `main.nf` is limited to study selection and DAG composition.
- `cell_type_analysis.nf`: a separate downstream entry workflow. It receives an
  existing merged single-cell H5AD, loads it once to split by retained per-study AIFI-L2,
  then fans out separate per-type analysis and lightweight report-rendering tasks.
- `modules/cell_type_analysis.nf`: the split, analysis, and report-rendering processes.
- `integration_benchmark.nf` and `modules/integration_benchmark.nf`: the
  optional global-integration task and its report-rendering task, both reading
  an existing raw merged single-cell H5AD.
- `src/pbmc_pipeline/integration_benchmark.py` and
  `reports/integration_benchmark_report.py`: benchmark fitting and its
  reader-facing notebook source.
- `reports/cell_type_report.py`: Jupytext source for the per-type notebook. It
  renders visualizations from the completed type-specific analysis artifact,
  fraction-model covariate and noise
  diagnostics, plus optional DE summaries and plots when the matching cell-type
  slug is present in the DE manifest. Its reader HTML leads with support
  metrics, a grouped table of contents, and a study-specific fraction-model
  evidence summary; it places coverage and the adjusted-fraction forest plot
  before secondary age trends, embedding, and marker detail. Forest-plot
  effects and residual SDs are displayed in percentage points, with confidence
  intervals extending beyond the axis marked at its boundary.
- `differential_expression.nf` and `modules/differential_expression.nf`: the
  independent pseudobulk DE workflow: label listing, one in-memory type subset
  and fit per task, then result/manifest collection.
- `src/pbmc_pipeline/differential_expression.py`: configured sample filtering,
  per-study and merged PyDESeq2 fits, task-result collection, CSV output, and
  the result manifest. It delegates age-bin model preparation and trajectory
  statistics to `src/pbmc_pipeline/age_trajectory.py`.
- `src/pbmc_pipeline/config.py`: configuration validation and the typed,
  required `AgeTrajectorySettings` view used by the DE fit and report stages.
- `src/pbmc_pipeline/trajectory_analysis.py` and
  `src/pbmc_pipeline/trajectory_report.py`: support-filtered per-type and
  shared-bin cross-type hierarchical clustering, UMAP coordinates, recurrence
  summaries, and the trajectory notebook/HTML report.
- `reports/trajectory_analysis_report.py`: source notebook for the combined
  age-trajectory report.
- `src/pbmc_pipeline/synthetic_de.py`: deterministic on-demand positive-fit
  pseudobulk fixture generation (28 samples per study, balanced across seven
  age decades below 90); the generated H5AD is ignored output, not a checked-in fixture.
- `scripts/check_synthetic_de.py`: assertions used by `make run-all-test` and
  `make verify` for age, sex, BMI, and CMV marker detection; covariate-specific
  study selection; complete-case counts; and report output. The fixture uses
  28 independent samples per study, balanced over the seven age decades below
  90, so it meets the production trajectory support threshold.
- `src/pbmc_pipeline/reporting.py`: testable report data transformations.
- `reports/merge_qc_report.py`: Jupytext source for the generated merge
  notebook. Report runners materialize these Python sources in a temporary
  location immediately before execution; generated templates are not tracked.
- `reports/qc_report.py`: the Jupytext source for the study-QC notebook.
  The runner materializes it in a temporary location before execution.

## Outputs

For `--outdir <outdir>`, the normal outputs are:

- `<outdir>/harmonized/<study>.h5ad`
- `<outdir>/prepared/<study>.cells.csv.gz`
- `<outdir>/prepared/<study>.source_obs.csv.gz` (retained source metadata audit sidecar)
- `<outdir>/prepared/<study>.prepare.json`
- `<outdir>/pseudobulk/<study>.pseudobulk.h5ad`
- `<outdir>/merged/pseudobulk_merged.h5ad`
- `<outdir>/merged/single_cell_merged.h5ad` (optional raw-count joined matrix)
- `<outdir>/reports/*.json`
- `<outdir>/qc/<study>/report.html`
- `<outdir>/qc/merged/report.html` (one report for both merge branches)
- `<outdir>/integration_benchmark/integration_benchmark.h5ad` (optional thin
  diagnostic artifact with labels and UMAP coordinates)
- `<outdir>/integration_benchmark/integration_benchmark.json` (optional method
  and input provenance)
- `<outdir>/integration_benchmark/{executed.ipynb,report.html}` (optional
  benchmark report)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/analysis.h5ad` (optional derived embeddings and clusters)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/report.html` (optional self-contained static, reader-facing downstream report without implementation-cell inputs)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/executed.ipynb` (optional executed technical/audit report)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/age_trajectory_{clusters,cluster_means}.csv` (optional per-gene labels, standardized profiles, UMAP coordinates, and cluster means)
- `<outdir>/cell_type_analysis/<aifi-l2-type>/fraction_model_diagnostics.tsv` (optional per-study fraction-model covariates, residual SD, and binomial-sampling reference)
- `<outdir>/cell_type_analysis/cell_type_manifest.json` (optional downstream provenance and completeness contract)
- `<outdir>/differential_expression/<cell-type-slug>/...csv` (per-study and combined covariate results)
- `<outdir>/differential_expression/<cell-type-slug>/combined/age_bin_pearson_residuals.npz` and `age_bin_pearson_residuals_samples.csv` (sample × gene residuals and sample metadata from the age-bin model)
- `<outdir>/differential_expression/<cell-type-slug>/age_model_diagnostics.csv` (sample support and raw-count inputs for the shared age-model diagnostic figure)
- `<outdir>/differential_expression/<cell-type-slug>/cell_type_result.json` (per-task fit record collected into the root manifest)
- `<outdir>/differential_expression/differential_expression.json` (DE status and cell-type result index)
- `<outdir>/differential_expression/trajectory_analysis/{report.html,executed.ipynb}` (combined trajectory report)
- `<outdir>/differential_expression/trajectory_analysis/*.csv` (shared-bin cross-type clusters and means, gene-level DE recurrence, and cross-type pattern concordance)
- `<outdir>/differential_expression/trajectory_analysis/expression_atlas_*.csv` (study-balanced common-gene expression matrix, depth support, gene-profile clusters, and independent 3′/5′ and intronic-read sensitivity summaries)
- `<outdir>/differential_expression/trajectory_analysis/*_residual_{clusters,pc_scores}.csv` (per-cell-type sample cluster assignments and PCA scores)
- `<outdir>/differential_expression/trajectory_analysis/cross_cell_type_pearson_residuals.npz` and residual metadata/cluster CSVs (wide matrix and analysis for shared samples across full-trajectory-supported cell types)
- `<outdir>/run_manifest.json` (requested, selected, and skipped studies;
  hashes of selected source inputs, metadata dependencies, and annotation models;
  selected studies are the complete requested set for explicit study lists)

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
`merged/single_cell_merged.h5ad`. It is disabled by default because the joined
matrix is large. An on-disk concat streams source matrices into a temporary
H5AD, then writes the raw joined counts and metadata. Per-cell QC fields
(`total_counts`, `n_genes_by_counts`, and `percent_mito`) are computed during
harmonization and carried through the join. The merge does not normalize the
full matrix, construct a global graph, or run another CellTypist prediction.

`make run-integration-benchmark` is a separate optional workflow over that
merged H5AD. It selects HVGs from the shared genes, runs PCA, applies Harmony
over `obs['study']`, and uses the resulting graph both for Harmony CellTypist
majority voting and UMAP. It then constructs one unintegrated PCA graph for the
comparison prediction. The published benchmark H5AD contains only metadata,
diagnostic labels, and UMAP coordinates, so it does not duplicate the raw
matrix. Its `benchmark_*_aifi_l2_majority` columns are diagnostic; the retained
per-study `obs['aifi_l2_majority']` remains the downstream label.
The accompanying integration-benchmark notebook renders the global UMAPs and
the AIFI-L2 concordance matrices formerly shown by merge QC.

### Cell-type analysis contract

`make run-cell-type-analysis` consumes the existing
`<outdir>/merged/single_cell_merged.h5ad`; it does not trigger the core merge.
The internal split task retains raw counts and canonical metadata. It computes
`n_cells_in_sample` and
`n_cells_in_sample_l1_parent` before splitting, so the reports need no auxiliary
inputs for their sample-level composition plots. The analysis task, rather than
the report notebook, computes the type-specific PCA/neighbours/UMAP/clusters.
It excludes configured V(D)J genes from local feature selection, uses
study-Harmony PCs for the local graph/UMAP/clusters, retains native PCs for the
descriptive PC diagnostics, and plots the inherited per-study L3 labels without
creating a merged L3 prediction. The renderer consumes the primitive split and
the completed analysis artifact, producing the executed notebook and HTML
without recomputing the analysis. Both the primitive split and the derived
analysis H5AD are published alongside the report tables.

The report visualizes native local PC1-versus-PC2 and PC3-versus-PC4 scores
with the standard study palette as colour and `technology` as marker shape.
These are descriptive, uncorrected-PC diagnostics; the Harmony-corrected PC
representation remains reserved for the local neighbour graph, UMAP, and
Leiden clustering.

The 50-cell local-analysis threshold does not guarantee enough cells per study
for reliable Harmony correction. The fixed 50-neighbour graph is also dense for
rare types. Local UMAP geometry and clusters must therefore be checked against
study colouring and composition, especially for small or imbalanced types.

The downstream manifest records checksums for its input, configuration,
reference models, and each published report. It is written only after validating
that every report emitted by the split fan-out has its H5AD, HTML, executed
notebook, and recognized completion status. `analysis_version` pins the
configuration-defined downstream analysis contract.

Cell-type splitting plus per-type analysis, analysis of an existing split, and
report rendering are named Nextflow subworkflows. The DE workflow remains an
independent producer: it lists labels, then each task loads the small merged
pseudobulk H5AD and subsets its assigned type in memory (without materializing
per-type pseudobulk inputs). A final collector writes the
`differential_expression.json` manifest that indexes per-study and combined
covariate results by cell-type slug. Each cell type also receives a separate
all-study age-bin trajectory model. It uses configured decade bins, a 20–30
reference, minimum sample support per bin, and a minimum of seven retained
decades; support decisions are recorded in cell-type metadata. The per-bin
coefficients are tested jointly with a Wald omnibus test and gene-level FDR
correction, then written with standardized shapes to `combined/age_bin.csv`.
Per-type reports cluster significant shapes with the configured linkage and
distance metric and plot cluster means and trajectory UMAPs. Gene assignments,
standardized profiles, and embedding coordinates are materialized beside each
per-type report, then checksummed in its downstream manifest. Types below the
seven-decade support threshold are omitted from trajectory figures, clustering, and cross-type
recurrence. A subsequent report process uses only bins shared by every
eligible type (at least five shared bins), clusters significant gene-by-type
trajectories, and summarizes per-gene recurrence and pattern agreement.
A single report-render process receives either the
matching result directory or a no-DE flag and publishes only beneath
`cell_type_analysis/`. DE results publish only beneath
`differential_expression/`.

Each study is fitted once with its maximal available design: age and sex plus
each optional covariate with at least two observed values, plus `study_site`
when it varies within the study, and `log10(total_counts)` as a technical
covariate. Its coefficients use one complete-case sample set.
Combined fits target one covariate at a time and include all complete-case
samples from studies that recorded it, adjusting for `study_site`, age, sex,
and log10 total counts. The combined shared model also reports the depth
coefficient as a volcano plot; age and depth have companion volcano panels
colored by log10 of the combined-fit `baseMean`.
The site term is omitted when only one site remains. Fit metadata and CSVs
record the formula, included studies and sites, complete-case sample count, and
exclusions for missing design values.

The explicit depth covariate is used in addition to PyDESeq2's size-factor
normalization. It was added after the production Naive CD8 T-cell age volcano
showed a strong age-associated abundance/depth shift, and sample-level residual
diagnostics showed patterns across pseudobulk depth after fitting. The observed
age–depth relationship could confound the age coefficient despite size-factor
normalization, so `log10(total_counts)` is included as a technical adjustment
in every per-study and combined fit, including the age-bin trajectory model.
This adjustment is motivated by those diagnostics; it is not evidence that
depth alone caused or fully explains the original pattern.

The configured `l2_parent_l1` mapping is an inferred taxonomy with recorded
provenance, rather than the independent `aifi_l1_majority` predictions. It is
applied only to retained per-study L2 labels; diagnostic benchmark L2 calls are
not used in this hierarchy. Validate the mapping against the AIFI atlas
before using within-L1 fractions for inference.

### Age-trajectory settings

The `differential_expression.age_trajectory` object in
`config/pipeline.json` is the source of truth for trajectory support,
significance, clustering, and display choices. Current settings are:

| Setting | Default | Meaning |
| --- | ---: | --- |
| `bin_width_years` | 10 | Width of the half-open age bins, such as `[20,30)`. |
| `minimum_samples_per_bin` | 10 | Minimum eligible sample-by-cell-type pseudobulks for a bin to be retained. |
| `minimum_bins` | 7 | Minimum retained bins for a cell type to enter trajectory displays and cross-type analyses. |
| `reference_bin_start_age` | 20 | Start age of the zero-reference bin; this is `[20,30)` with the current width. |
| `strict_age_cutoff_exclusive` | 90 | Samples with age ≥90 are excluded before trajectory binning and residual fitting. |
| `cross_cell_type_residual_minimum_sample_coverage` | 0.8 | Target fraction of the union of fitted samples present in every cell type retained for the combined residual embedding. |
| `cluster_fdr_threshold` | 0.001 | Omnibus adjusted-p-value cutoff for trajectories entering clustering. |
| `de_fdr_threshold` | 0.05 | Omnibus adjusted-p-value cutoff for cross-type significance and recurrence summaries. |
| `max_clusters` | 5 | Upper bound on hierarchical clusters. |
| `linkage_method` | `ward` | Hierarchical clustering linkage. Supported values are `single`, `complete`, `average`, `weighted`, `centroid`, `median`, and `ward`. |
| `distance_metric` | `euclidean` | Profile distance. Supported values are `euclidean`, `cityblock`, `cosine`, and `correlation`; Ward, centroid, and median require Euclidean distance. |
| `minimum_shared_bins` | 5 | Minimum common age bins required for the cross-cell-type trajectory comparison. |
| `residual_pca_components` | 100 | Maximum number of centered PCA components used before clustering bins-model Pearson residuals by sample. The count is capped by the available samples and genes. |
| `residual_umap_neighbors` | 50 | Neighbourhood size for residual-sample UMAPs; capped at one fewer than the available samples. |
| `umap_neighbors` | 15 | UMAP neighborhood size. |
| `minimum_umap_trajectories` | 4 | Minimum number of profiles required to compute UMAP coordinates. |
| `umap_min_dist` | 0.15 | UMAP minimum-distance parameter. |
| `random_state` | 413 | Seed for reproducible UMAP coordinates. |
| `report_top_n_genes` | 30 | Maximum gene rows shown in a per-cell-type trajectory table. |
| `cross_report_top_n_genes` | 100 | Maximum gene rows shown in cross-cell-type trajectory tables. |

The 90–100 age bin is explicitly removed before trajectory fitting. Its sample
count is measured before filtering and recorded per cell type in run_metadata.json;
per-cell-type reports and the combined trajectory support table display that count.
The strict age cutoff also excludes any eligible samples aged 100 or older.

The omnibus test asks whether any retained age-bin coefficient differs from
the reference bin; it tests any age-bin difference, including a monotone
pattern, rather than only curvature beyond a linear trend. Standardized shapes
are computed only when every retained bin is estimable and the across-bin
population SD is non-zero. The trajectory report also clusters samples from
each cell type's bins-model Pearson residual matrix. For observed counts `y`,
fitted means `mu`, and PyDESeq2 dispersions `alpha`, it stores Pearson residuals
`(y - mu) / sqrt(mu + alpha * mu^2)` as float32 values in compressed NPZ files;
sample identifiers and gene names are included, with sample covariates stored
beside each matrix. The age-bin fit excludes samples aged 90 years or older,
requires at least ten eligible pseudobulks in each retained bin, and adjusts for
age-bin, sex, log10(total counts), and study site when that term is estimable. Samples are clustered
using up to `residual_pca_components` centered PCA scores, then hierarchical
clustering and UMAP use Euclidean distances in PC space. Each per-type report
shows a 2 × 2 UMAP grid for sex, age, log10(total counts), and study, plus a
second 2 × 2 grid for residual clusters, BMI, and CMV; BMI and CMV panels show
recorded-value counts.

The cross-type matrix starts with every cell type that has residual artifacts,
then removes the type with the fewest fitted samples until samples present in
every remaining type reach the configured fraction of the union sample set, or
one type remains. Its report names included and removed types, the shared sample
count, and achieved coverage. It contains a 2 × 2 UMAP grid for model covariates
and a second grid for BMI and CMV coverage. Wide columns are labeled
`cell type::gene`. Per-type and cross-type cluster assignments, PC scores,
matrices, and sample metadata are materialized for downstream analysis.

### Expression-atlas settings

`differential_expression.expression_atlas` configures a descriptive,
cross-cell-type expression summary in the combined trajectory report. It uses
the gene intersection across all studies, pools raw counts to one study ×
AIFI-L2 cell-type profile, applies the depth filter to that exact intersection,
then calculates `log2(CPM + pseudocount)` and averages log values equally across
retained studies. This is a relative expression measure, not calibrated
absolute RNA abundance or a covariate-adjusted model.

The current depth threshold is 1,000,000 counts, `cpm_pseudocount` is 1, and a
cell type needs two retained studies. Genes with at least 1 CPM in one retained
cell type enter average-linkage/correlation profile clustering, capped at eight
clusters. The report separately contrasts the study-balanced 5′ and 3′
matrices (`5′ − 3′`) and the intronic-read and non-intronic matrices
(`intronic − non-intronic`). Both are calculated once each side has at least
one depth-qualified study in a cell type. `minimum_studies_per_technology` is
the two-study-per-side replication flag shown alongside each contrast, not an
estimability gate: one-versus-one contrasts remain visible but are explicitly
unreplicated. The two contrasts are non-interaction diagnostics; they do not
split into technology × intronic groups. Both labels can be study-confounded,
so neither contrast identifies a causal technical effect.

The positive synthetic DE fixture must retain at least seven supported bins
under these production settings. It therefore creates 28 independent samples
per study and cell type: four in each of seven decades below 90 across three
studies, giving 12 combined samples per bin before any covariate-specific exclusions.
The generator derives its sample count from `minimum_samples_per_bin` and the
number of synthetic studies, then records the active trajectory settings and
planted profiles in the fixture metadata. Keep the synthetic workflow report
checks enabled when changing trajectory support or filtering behavior. Its raw
count scale makes each study × cell-type group pass the expression atlas's 1M
gate; it labels two studies as 3′ and one as 5′, which exercises the atlas but
only supports an explicitly unreplicated 3′/5′ contrast. It also labels one
study as intronic and one as non-intronic for an independent, explicitly
unreplicated intronic-read contrast.

## Running and verification

Use Make targets rather than raw Nextflow or Docker commands:

```bash
make help
make lint test-unit
make workflow-lint
make run-test STUDIES=wang25
make run-all-test
make verify
make run-all STUDIES=all
make split-cell-types-test
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
make render-trajectory-report-test
make run-trajectory-test
make docs-check
```

Use the smallest check that covers the edit: `make lint test-unit` for
Python-only work, `make workflow-lint` for workflow/configuration edits, and a
focused `make run-test STUDIES=<study>` when core execution wiring needs
testing. For a single cell-type report, reuse an existing split with
`make run-cell-type-test CELL_TYPE=<slug>` rather than rebuilding the core
workflow. When only the report template changes, use
`make render-cell-type-test CELL_TYPE=<slug>` to consume the published split
and analysis artifact without recomputing local analysis. `make docs-check`
validates local Markdown links and documented Make targets without adding a
documentation-tool dependency.
`make render-trajectory-report[-test]` reruns the cross-cell-type analysis and
notebook from an existing DE manifest; it does not refit any model or rerender
the per-cell-type reports.
`make run-trajectory-test` reuses an existing test single-cell merge (run
`make run-test` if it is missing), fits the config-sized synthetic DE fixture,
then renders and checks per-cell-type and cross-cell-type trajectory reports.
It omits the separate integration benchmark stage from `make run-all-test`.
`make verify` runs linting, unit tests, Nextflow lint, cached image builds, and
non-resumed Docker workflows across all configured studies by default. After
the core merge it runs and checks the integration benchmark report, then
generates synthetic pseudobulk data for positive DE fits, checks the completed
models and planted markers, and renders all cell-type reports from the fresh
core test merge. Reserve it for cross-cutting or release-level validation;
pass `STUDIES=<list>` for a narrower verification. The synthetic fixture has
only two cell types, so other reports exercise the ordinary no-DE path. It is
workflow coverage rather than a realistic model of study balance or biological
effects: covariate availability and gene presence are deliberately simple, and
the planted effects have no study-specific heterogeneity.
`make run` is the resumable production entry point; `make run-no-qc` omits only
report rendering.
`MERGE_SINGLE_CELL=true` enables the optional single-cell branch.
The current Nextflow profiles leave CPU and time requests unspecified except
for one CPU per parallel cell-type DE task; they do not set per-process memory
limits. The full production pipeline and a fresh-clone input acquisition plus
`make run-all` completed successfully before the latest manual metadata
corrections. Resource reports/traces are not yet part of the tracked run
artifacts; record measured peak RAM before sizing a different runner or adding
resource directives.

## Per-cell-type reports

The main fraction model uses adults only and independently fits each
study with eligible age, sex, BMI, and CMV covariates. Age effects are reported
per decade. The modeled response is a 0–1 fraction; the forest plot scales
coefficients, confidence limits, and residual SDs to percentage points.
`cell_type_analysis.fraction_model` configures the deterministic
parametric binomial bootstrap (currently 1,000 replicates and seed 0) used as
the conditional cell-sampling reference in the report and diagnostics TSV.
Categorical forest-plot labels show each contrast and its reference level;
sex uses female as the reference, and CMV uses no/negative as the reference.
Each donor's OLS-fitted fraction is held fixed as its binomial probability
(clipped to `[0, 1]`, with the count and fraction reported in diagnostics), its
observed parent-cell count supplies the number of trials, and the same
covariate design is refit in each replicate.

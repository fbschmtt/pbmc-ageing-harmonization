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
                                      ├─ integration benchmark → method embeddings, graphs, labels, report
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
  optional global-integration task and its report-rendering task, both reading
  an existing raw merged single-cell H5AD.
- `src/pbmc_pipeline/integration_benchmark.py` and
  `reports/integration_benchmark_report.py`: benchmark fitting and its
  reader-facing notebook source.
- `reports/cell_type_report.py`: Jupytext source for the per-type notebook, rendering visualizations from the completed type-specific analysis artifact.
- `reports/differential_expression_report.py`: standalone Jupytext source for
  the per-type pseudobulk DE report. Rendering executes abundance models,
  trajectory and residual analyses, and residual WGCNA from saved DE outputs.
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
  shared-bin cross-type hierarchical clustering of gene trajectories,
  cross-type residual PCA/UMAP and sample clustering, recurrence summaries,
  and data for the combined cross-cell-type report.
- `reports/trajectory_analysis_report.py`: source notebook for the combined
  cross-cell-type analysis report rendered by `run-de`.
- `scripts/render_report_index.py`: builds the lightweight HTML index from
  harmonization QC, integration benchmark, the combined cross-cell-type report,
  and the residual and DE reports currently published per type.
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
  `scripts/profile_top_age_genes.py` applies the same report to the most
  significant merged age-DE genes per cell type.
- `scripts/check_gene_profile.py`: output-contract and notebook smoke checks
  for the focused gene report. `make run-gene-profile-synthetic-test` fits and
  renders a planted synthetic gene profile; `make profile-top-age-genes-test`
  renders two top merged age-DE genes per cell type from the synthetic DE
  fixture. Both report workflows are included in `make run-all-test`.
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

## Outputs

For `--outdir <outdir>`, the normal outputs are:

- `<outdir>/harmonization/harmonized/<study>.h5ad`
- `<outdir>/harmonization/prepared/<study>.cells.csv.gz`
- `<outdir>/harmonization/prepared/<study>.source_obs.csv.gz` (retained source metadata audit sidecar)
- `<outdir>/harmonization/prepared/<study>.prepare.json`
- `<outdir>/harmonization/pseudobulk/<study>.pseudobulk.h5ad`
- `<outdir>/harmonization/merged/pseudobulk_merged.h5ad`
- `<outdir>/harmonization/merged/single_cell_merged.h5ad` (optional raw-count joined matrix)
- `<outdir>/harmonization/provenance/*.json` (preparation, harmonization, and merge run summaries)
- `<outdir>/harmonization/harmonization_qc/<study>/report.html`
- `<outdir>/harmonization/harmonization_qc/merged/report.html` (one report for both merge branches)
- `<outdir>/harmonization/integration_benchmark/integration_benchmark.h5ad` (optional thin
  diagnostic artifact with labels, method embeddings, and neighbor graphs)
- `<outdir>/harmonization/integration_benchmark/integration_benchmark.json` (optional method
  and input provenance)
- `<outdir>/harmonization/integration_benchmark/{executed.ipynb,report.html}` (optional
  benchmark report)
- `<outdir>/analysis/per_type/<cell-type-slug>/split.h5ad` (reusable raw-count split)
- `<outdir>/analysis/per_type/<cell-type-slug>/analysis/{analysis.h5ad,report.html,executed.ipynb}` (optional cell-type analysis and reports)
- `<outdir>/index.html` (links to harmonization QC, integration benchmark, combined cross-type analysis, and available residual and DE reports per cell type)
- `<outdir>/analysis/per_type/cell_type_manifest.json` (optional downstream provenance and completeness contract)
- `<outdir>/analysis/combined/differential_expression/*.csv` (combined cross-cell-type DE results and manifest)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/...csv` (per-study and per-type combined covariate results)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/combined/age_bin_pearson_residuals.npz` and `age_bin_pearson_residuals_samples.csv` (sample × gene residuals and sample metadata from the age-bin model)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/age_model_diagnostics.csv` (sample support and raw-count inputs for the shared age-model diagnostic figure)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/developmental/age_support_samples.csv` (samples retained by the [0,20) developmental model, used for its study-site age-support plot)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/cell_type_result.json` (per-task fit record collected into the root manifest)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/{abundance_input.tsv,report.html,fraction_model_diagnostics.tsv}` (pseudobulk abundance input and per-type DE report outputs)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/residual_wgcna/` (residual WGCNA module assignments, eigengenes, soft-threshold and sample-connectivity diagnostics, covariate associations, hub genes, and leave-one-study-out module stability)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/residual_wgcna/{module_gsea.csv,cmv_associated_module_scores.csv,smoking_associated_module_scores.csv}` (preranked GO enrichment and sample scores for the strongest exposure-associated modules)
- `<outdir>/analysis/combined/differential_expression/differential_expression.json` (DE status and cell-type result index)
- `<outdir>/analysis/combined/trajectory/{report.html,executed.ipynb}` (combined trajectory report)
- `<outdir>/analysis/combined/trajectory/*.csv` (shared-bin cross-type clusters and means, gene-level DE recurrence, and cross-type pattern concordance)
- `<outdir>/analysis/combined/trajectory/developmental_adult_age_effects.csv` (per-cell-type paired raw continuous-age coefficients, mean normalized expression, and adjusted p-values for the [0,20) developmental and ≥20 adult fits; the combined report plots genes retained by both models' independent expression filters)
- `<outdir>/analysis/combined/trajectory/developmental_adult_age_effect_support.csv` (sample and gene support for that comparison)
- `<outdir>/analysis/combined/trajectory/expression_atlas_*.csv` (study-balanced common-gene expression matrix, gene UMAP, depth support, and independent 3′/5′ and intronic-read sensitivity summaries)
- `<outdir>/analysis/combined/trajectory/*_residual_{clusters,pc_scores}.csv` (per-cell-type sample cluster assignments and PCA scores)
- `<outdir>/analysis/combined/trajectory/cross_cell_type_pearson_residuals.npz` and residual metadata/cluster CSVs (wide matrix, Leiden-primary assignments, K-means/HDBSCAN alternatives, marker contrasts, and cell-type feature composition for shared samples across full-trajectory-supported cell types)
- `<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/` (focused gene-profile HTML and executed notebook, `report_render.log` for notebook-render diagnostics, model coefficient and residual tables, mean–dispersion diagnostics, normalized sample counts, and compressed VST matrix)
- `<outdir>/run_manifest.json` (requested, selected, and skipped studies;
  hashes of selected source inputs, metadata dependencies, and annotation models;
  selected studies are the complete requested set for explicit study lists)

### Harmonized expression contract

Each harmonized H5AD has exactly one expression representation: `.X` contains
the configured unnormalized integer count matrix. `.raw` is `None` and
`.layers` is empty. Normalized matrices may be created transiently for
annotation or embeddings, but are never retained in the output. This prevents
an alternate normalized matrix from being mistaken for analysis input.
For CellTypist harmonization, the validated count object is written to the final
H5AD before normalization. Analysis then uses the in-memory `.X`; after
annotation, that normalized matrix is released and the H5AD's metadata and
embeddings are updated in place while preserving the saved count matrix.

Merged pseudobulk `.var` records per-study gene availability. A gene absent
from an input study is represented by an outer-join zero only for that study;
use `available_in_<study>`, `n_studies_with_gene`,
`synthetic_zero_filled_studies`, and `has_synthetic_zeros` to distinguish those
synthetic zeros from observed zero counts.

The `--merge_single_cell` branch writes
`merged/single_cell_merged.h5ad`. The supported `make run` target enables it by
default; set `MERGE_SINGLE_CELL=false` to skip the merge. On-disk concat streams
the source matrices directly into the planned output H5AD. Temporary concat
inputs omit unused source embeddings and analysis metadata; non-CSC inputs use
external-link views, while CSC inputs are converted to CSR one at a time.
Per-cell QC fields
(`total_counts`, `n_genes_by_counts`, and `percent_mito`) are computed during
harmonization and carried through the join. The merge does not normalize the

`make run-integration-benchmark` is a separate optional workflow over that
merged H5AD. It opens the input in backed mode and defaults to processing all
cells up to `integration_benchmark.max_cells` (20,000,000). Above that cap, it
materializes a uniform random sample without replacement using
`processing.random_seed`, shared across all methods. Sample size, source size,
fraction, selection strategy, and seed are recorded in the JSON sidecar,
artifact metadata, and report. Benchmark UMAPs and comparisons therefore
describe all input cells unless the configured cap is exceeded. The active
benchmark configuration compares an unintegrated PCA baseline with Harmony
and BBKNN using one normalized HVG/PCA workspace; Scanorama is temporarily
disabled because its correction step can allocate large dense arrays. Its PCA preparation is shared
through `pbmc_pipeline.embedding.prepare_sparse_pca`: mean-dispersion HVGs are
selected using the configured minimum mean and dispersion cutoffs with no
maximum mean cutoff, and V(D)J genes are excluded case-insensitively,
and if fewer than two genes remain, all non-V(D)J genes are used. The selected
matrix is cast to float32 and scaled with centering disabled to preserve
sparsity; ARPACK PCA applies centering implicitly. PCA computes up to
`processing.neighbors.n_pcs`, capped by the available cells and selected genes.
Harmony adjusts embeddings, BBKNN balances the neighbor graph, and the
unintegrated PCA graph is the baseline. The benchmark implementation still
supports Scanorama if it is re-enabled in the configuration.
Every method's graph supplies both CellTypist majority voting
and its UMAP. The integration report compares graph edge overlap with
unintegrated PCA and AIFI-L2 label concordance without applying acceptance
thresholds. The published benchmark H5AD contains metadata, diagnostic
labels, method embeddings, and neighbor graphs, but no expression matrix. Its
`benchmark_*_aifi_l2_majority` columns are diagnostic; the retained per-study
`obs['aifi_l2_majority']` remains the downstream label. The notebook separates
embedding integration from type calling and reports graph edge overlap and
AIFI-L2 concordance.
The benchmark emits timestamped stage progress to the task's captured stdout
and stderr files in its Nextflow `work/` directory. Production and test
targets use the standard Python image. Use
`INTEGRATION_BENCHMARK_LOG_LEVEL=DEBUG make run-integration-benchmark` for
additional diagnostics.

### Cell-type analysis contract

`make run-cell-type-analysis` consumes the existing
`<outdir>/harmonization/merged/single_cell_merged.h5ad`; it does not trigger the core merge.
The split task reads the merged H5AD in backed read-only mode and materializes
one selected type at a time before writing its raw-count split. The complete
sparse expression matrix remains on disk; AnnData still loads the full
observation metadata into memory, and the largest individual type must fit in
memory during the write. The split retains raw counts and canonical metadata.
Per-sample abundance denominators are derived later from the merged pseudobulk
object.
The analysis task writes raw counts to the planned `analysis/analysis.h5ad`
before normalizing `.X` in memory. After PCA, Harmony neighbours/UMAP, Leiden
clustering, and marker ranking, it releases the normalized `.X` and updates
metadata and embeddings in that H5AD, preserving the on-disk counts without an
in-memory raw-count copy. PCA preparation uses the shared
`pbmc_pipeline.embedding.prepare_sparse_pca` function described above. Stage
orchestration, graph construction, Harmony, clustering, and markers are in
`pbmc_pipeline.cell_type.analyse_cell_type`; harmonization orchestration is in
`pbmc_pipeline.harmonize._compute_embedding`.
During report execution, the notebook also
computes a comparison UMAP from the saved native PCs; it is not persisted back
to the analysis H5AD.
It excludes configured V(D)J genes from local feature selection, uses
study-Harmony PCs for the local graph/UMAP/clusters, retains native PCs for the
descriptive PC diagnostics, and plots the inherited per-study L3 labels without
creating a merged L3 prediction. The renderer consumes the primitive split and
the completed analysis artifact, producing the executed notebook and HTML
without repeating the local PCA, Harmony clustering, or marker analysis. It
does compute the additional unintegrated UMAP diagnostic. Both the primitive
split and the derived analysis H5AD are published alongside the report tables.

The report visualizes native local PC1-versus-PC2 and PC3-versus-PC4 scores
with the standard study palette as colour and `technology` as marker shape.
These are descriptive, uncorrected-PC diagnostics; the Harmony-corrected PC
representation remains reserved for the local neighbour graph, UMAP, and
Leiden clustering. Per-type reports show one study-coloured UMAP from native
PCs before the existing Harmony-integrated UMAP plots.

The 50-cell local-analysis threshold does not guarantee enough cells per study
for reliable Harmony correction. The fixed 50-neighbour graph is also dense for
rare types. Local UMAP geometry and clusters must therefore be checked against
study colouring and composition, especially for small or imbalanced types.

The downstream manifest records checksums for its input, configuration,
reference models, and each published report. It is written only after validating
that every report emitted by the split fan-out has its H5AD, HTML, executed
notebook, and recognized completion status. `analysis_version` pins the
configuration-defined downstream analysis contract.
The configured `l2_parent_l1` taxonomy supplies the parent label in
`report.json` and the cell-type manifest; it is not copied into split or
analysis H5AD observations. Per-sample cell-type and parent denominators are
built from pseudobulk `n_cells` in the DE workflow.

Cell-type splitting plus per-type analysis, analysis of an existing split, and
report rendering are named Nextflow subworkflows. The DE workflow remains an
independent producer: it lists labels, then each task loads the small merged
pseudobulk H5AD and subsets its assigned type in memory (without materializing
per-type pseudobulk inputs). A final collector writes the
`differential_expression.json` manifest that indexes per-study and combined
covariate results by cell-type slug. It publishes only each type's
covariate results by cell-type slug. Each cell type also receives a separate
all-study age-bin trajectory model. It uses configured decade bins, a 20–30
reference, minimum sample support per bin, and a minimum of seven retained
decades; support decisions are recorded in cell-type metadata. The per-bin
coefficients are tested jointly with a Wald omnibus test and gene-level FDR
correction, then written with standardized shapes to each type's
`differential_expression/combined/age_bin.csv`. Types below the seven-decade
support threshold are omitted from trajectory figures, clustering, and
cross-type recurrence. The combined cross-cell-type report uses only bins
shared by every eligible type (at least five shared bins), clusters significant
gene-by-type trajectories, and summarizes per-gene recurrence and pattern
agreement. Its cluster tables and means are published under
`analysis/combined/trajectory/`.
The single cell analysis report process consumes only the split and single-cell analysis
artifacts. It has no DE manifest dependency. The DE workflow separately
renders each per-type DE report after the collector completes. Cell-type
analysis and reports publish under
`analysis/per_type/<cell-type-slug>/analysis/`; per-type DE results publish
under `analysis/per_type/<cell-type-slug>/differential_expression/`. The
collector writes the shared manifest under
`analysis/combined/differential_expression/`, and the cross-type report is
under `analysis/combined/trajectory/`. The single cell analysis report does
not require DE results, while `run-de` renders a separate DE report for each
cell type from the pseudobulk outputs.

Each study is fitted once with its maximal available design: age and sex plus
each optional covariate with at least two observed values, plus `study_site`
when it varies within the study, and `log10(total_counts)` as a technical
covariate. Its coefficients use one complete-case sample set.
Combined fits target one covariate at a time and include all complete-case
samples from studies that recorded it, adjusting for `study_site`, age, sex,
and log10 total counts. The per-type DE report presents fit status, covariate
volcano plots, age-bin trajectories, residual clustering, residual WGCNA,
abundance models, and result tables, including the depth coefficient where
estimable. Report rendering executes the abundance OLS/bootstrap, trajectory
and residual clustering, and WGCNA calculations from saved DE outputs.
Rerendering skips PyDESeq2 fitting but reruns those report-stage calculations.
The site term is omitted when only one site remains. Fit metadata and CSVs
record the formula, included studies and sites, complete-case sample count, and
exclusions for missing design values.

The explicit depth covariate is used in addition to PyDESeq2's size-factor
normalization. It adjusts for the age-associated pseudobulk-depth pattern
observed in Naive CD8 T cells and is included as a technical covariate in every
per-study and combined fit, including the age-bin trajectory model.

The configured `l2_parent_l1` mapping is an inferred taxonomy with recorded
provenance, rather than the independent `aifi_l1_majority` predictions. It is
applied only to retained per-study L2 labels; diagnostic benchmark L2 calls are
not used in this hierarchy. Validate the mapping against the AIFI atlas
before using within-L1 fractions for inference.

### Age-trajectory settings

The `differential_expression.age_trajectory` object in
`config/pipeline.json` is the source of truth for trajectory support,
significance, clustering, and display choices. The settings table below is
generated from that object with `make docs-update` and
checked for drift by `make docs-check`.

<!-- generated-age-trajectory-settings:start -->

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
| `max_clusters` | 7 | Upper bound on trajectory hierarchical clusters and candidate residual K-means groups. |
| `linkage_method` | `"ward"` | Hierarchical clustering linkage. Supported values are `single`, `complete`, `average`, `weighted`, `centroid`, `median`, and `ward`. |
| `distance_metric` | `"euclidean"` | Profile distance. Supported values are `euclidean`, `cityblock`, `cosine`, and `correlation`; Ward, centroid, and median require Euclidean distance. |
| `umap_neighbors` | 15 | UMAP neighborhood size. |
| `minimum_umap_trajectories` | 4 | Minimum number of profiles required to compute UMAP coordinates. |
| `umap_min_dist` | 0.15 | UMAP minimum-distance parameter. |
| `report_top_n_genes` | 30 | Maximum gene rows shown in a per-cell-type trajectory table. |
| `cross_report_top_n_genes` | 100 | Maximum gene rows shown in cross-cell-type trajectory tables. |
| `random_state` | 413 | Seed for reproducible UMAP coordinates. |
| `minimum_shared_bins` | 5 | Minimum common age bins required for the cross-cell-type trajectory comparison. |
| `residual_pca_components` | 100 | Maximum number of centered PCA components used before residual UMAP. The count is capped by the available samples and genes. |
| `residual_umap_neighbors` | 50 | Neighbourhood size for residual-sample UMAPs; capped at one fewer than the available samples. |
| `residual_hdbscan_min_cluster_size` | 15 | Smallest sample group retained by the HDBSCAN alternatives; points outside selected clusters are labeled noise. |
| `residual_hdbscan_min_samples` | 5 | Local-density conservativeness, used for both HDBSCAN comparison embeddings. |
| `residual_hdbscan_cluster_selection_method` | `"eom"` | HDBSCAN selection strategy, used for both UMAP and PC-space comparison labels. No HDBSCAN parameter sweep is shown. |
| `residual_pc_hdbscan_components` | 100 | Maximum number of leading residual PCs used for the HDBSCAN-in-PC-space alternative and primary Leiden graph. |
| `residual_leiden_neighbors` | 15 | Number of PC-space nearest neighbors used to construct the Leiden graph; capped at one fewer than the available samples. |
| `residual_leiden_resolutions` | `[0.25, 0.5, 1.0, 1.5]` | Candidate Leiden resolutions tested on the same PC-space graph. |
| `residual_leiden_min_clusters` | 10 | Target used to select the lowest candidate resolution yielding at least this many clusters; if none do, use the candidate yielding the most clusters. |
| `residual_marker_max_clusters` | 16 | Maximum number of largest clusters included in marker and feature-origin contrasts. |
<!-- generated-age-trajectory-settings:end -->

The 90–100 age bin is explicitly removed before trajectory fitting. Its sample
count is measured before filtering and recorded per cell type in
`run_metadata.json`; DE metadata and the combined trajectory support table
display that count. The strict age cutoff also excludes samples aged 100 or
older.

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
using up to `residual_pca_components` centered PCA scores. UMAP uses Euclidean
distance in PC space. Primary residual labels use Leiden on a PC-space
nearest-neighbor graph. The selected resolution is the lowest configured
candidate yielding at least `residual_leiden_min_clusters` groups; if no
candidate reaches the target, the candidate with the most groups is selected.
The combined cross-cell-type report shows candidate cluster counts and compares the selected result with
silhouette-selected K-means on UMAP and one configured HDBSCAN result on UMAP
and on up to the configured number of leading PCs. The HDBSCAN comparisons
both use the same `min_cluster_size`, `min_samples`, and selection method;
those parameters are not swept. Marker and feature-origin contrasts include at most the
`residual_marker_max_clusters` largest clusters; all labels remain in the
exported assignments. The combined report contains the sample-level residual
UMAP/PCA, clustering comparisons, and cross-cell-type cluster composition.
The single cell analysis reports instead show the local cell-level embedding,
markers, and fraction views described in the per-type report workflow.

The cross-type matrix starts with every cell type that has residual artifacts,
then removes the type with the fewest fitted samples until samples present in
every remaining type reach the configured fraction of the union sample set, or
one type remains. Its report names included and removed types, the shared sample
count, and achieved coverage. It contains a 2 × 2 UMAP grid for model covariates
and a second grid for BMI and CMV coverage; only the non-model grid is open by
default. It also plots unique sample counts for each technical batch within
each study site, with full counts in `technical_batch_sample_counts.csv`. A
collapsed PCA panel shows PC1 versus PC2 by cluster and the residual
UMAP colored by PC1 and PC2. It shows a clustering summary table and UMAP
overlays for the K-means and two HDBSCAN alternatives; the table includes the
tested Leiden resolution counts. Wide columns are labeled
`cell type::gene`. Per-type and cross-type cluster assignments, PC scores,
matrices, and sample metadata are materialized for downstream analysis.

Each per-type DE report runs residual WGCNA from the fitted age-bin residual
matrix. It publishes module assignments and eigengenes together with
soft-threshold, sample-connectivity, eigengene-covariate, hub-gene, CMV and
smoking checks, preranked GO Biological Process enrichment, and leave-one-study-out
stability diagnostics. CMV and smoking are checked because the residual model
already adjusts for age bins, sex, study site, and depth but omits those plausible
sources of remaining biological structure. These exploratory checks ask whether
signal remains; they do not assume an effect is present or adjust the residuals.
WGCNA color names are arbitrary module identifiers; grey is reserved for genes
not assigned to a module and is excluded from GSEA because it is not a coherent
module. The stability panel compares gene membership after omitting each study.

### Expression-atlas settings

`differential_expression.expression_atlas` configures a descriptive,
cross-cell-type expression summary in the combined trajectory report. It uses
the gene intersection across all studies, pools raw counts to one study ×
AIFI-L2 cell-type profile, applies the depth filter to that exact intersection,
then calculates `log2(CPM + pseudocount)` and averages log values equally across
retained studies. The report records the outer-union gene availability matrix
and overlap display next to this atlas, along with per-study counts outside the
intersection and each study's largest excluded gene, so the intersection is
interpreted where it is actually applied. This is a relative expression measure, not calibrated
absolute RNA abundance or a covariate-adjusted model.

The current depth threshold is 1,000,000 counts, `cpm_pseudocount` is 1, and a
cell type needs two retained studies. Genes with at least 1 CPM in one retained
cell type enter the gene-expression UMAP. Each gene is centered and scaled only
for its PCA (up to 100 components) and UMAP (50 neighbors); the report shows
the unscaled mean expression and the per-gene SD used for that scaling, then
uses the scaled values only for its requested AIFI-L1 summaries (mean across
each L1 class's retained L2 types). The report separately contrasts the study-balanced 5′ and 3′
matrices (`5′ − 3′`) and the intronic-read and non-intronic matrices
(`intronic − non-intronic`). Both are calculated once each side has at least
one depth-qualified study in a cell type. `minimum_studies_per_technology` is
the two-study-per-side replication flag shown alongside each contrast, not an
estimability gate: one-versus-one contrasts remain visible but are explicitly
unreplicated. The two contrasts are non-interaction diagnostics; they do not
split into technology × intronic groups. Both labels can be study-confounded,
so neither contrast identifies a causal technical effect.

The positive synthetic DE fixture must retain at least seven supported bins
under these production settings. It therefore creates 36 independent samples
per study and cell type: four in each of two developmental decades and seven
adult decades below 90 across three studies, giving 12 combined samples per
bin before any covariate-specific exclusions.
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
make docs-update
make docs-check
```

`make test-data`, `make run-test`, `make run-all-test`, and `make verify`
default to the representative `aida25,cima,terekhova26` set. It covers the
AIDA supplementary-metadata join, CIMA counts-layer input, and Terekhova26
raw-count input. Pass `STUDIES=all` to exercise every available fixture, or a
comma-separated list for another focused selection. Production `make run`
continues to default to all studies.

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
`run-de` fits pseudobulk models, writes the combined cross-cell-type analysis
report, then renders each per-type DE report. With `-resume`, completed PyDESeq2
fits remain cached when only a downstream report input changes. Use
`make render-de-reports` to rerun per-type report notebooks from existing model
outputs; this skips PyDESeq2 fitting but reruns analyses executed inside those
notebooks.
Core, integration-benchmark, DE, and cell-type report targets update
`<outdir>/index.html` after publishing. `make render-report-index OUTDIR=<outdir>`
rebuilds the index using the cell-type manifest when available, or scanning
published per-type reports otherwise, which keeps its links synchronized when
cell types change. Harmonization and benchmark reports are linked when their
published files are present.
`make verify` runs linting, unit tests, Nextflow lint, image builds, and fresh
Docker workflows for the representative test set by default. After the core
merge it runs and checks the integration benchmark report, then generates
synthetic pseudobulk data for positive DE fits, checks the completed models and
planted markers, profiles the top two merged age-DE genes per cell type,
exercises the focused synthetic gene-profile report, renders per-type DE
reports from synthetic data, and renders single cell analysis reports from the
fresh core test merge. Reserve it for cross-cutting or release-level validation;
pass `STUDIES=all` for every available fixture or `STUDIES=<list>` for another
focused verification. The synthetic fixture has
only two cell types, so other reports exercise the ordinary no-DE path. It is
workflow coverage rather than a realistic model of study balance or biological
effects: covariate availability and gene presence are deliberately simple, and
the planted effects have no study-specific heterogeneity.
Its synthetic WGCNA output currently has two modules. Thus verification executes
the ordinary WGCNA and leave-one-study-out report path, but does not exercise
the 25-module display cap or provide a production-scale runtime benchmark.
`make run` is the resumable production entry point and merges single-cell data
by default. Set `MERGE_SINGLE_CELL=false` to skip that optional branch.
The current Nextflow profiles leave CPU and time requests unspecified except
for one CPU per parallel cell-type DE task; they do not set per-process memory
limits. Resource reports and traces are not yet part of the tracked run
artifacts; record measured peak RAM before sizing a different runner or adding
resource directives.

## Per-type abundance in DE reports

Each per-type DE task writes an `abundance_input.tsv` with one row for every
sample observed in the merged pseudobulk input, including zero-count
sample/type combinations. Pseudobulk `n_cells` supplies the type numerator and
the retained-PBMC and parent-compartment denominators. These reports therefore
do not require per-cell abundance columns in split H5ADs. The per-type DE
inclusion filter does not exclude zero-count abundance rows.

The report fits adult-only, study-specific OLS models with the existing age,
sex, BMI, CMV, and smoking-status covariate rules; age effects are reported per
decade, and forest-plot effects and residual SDs are shown in percentage points. Each
fraction view independently requires at least the configured denominator
count (currently 10 cells). The deterministic parametric binomial bootstrap
uses 1,000 replicates and seed 0. For each sample, the OLS-fitted fraction is
held fixed as the binomial probability (clipped to `[0, 1]`), and the observed
parent-cell count supplies the number of trials. The same covariate design is
refit in each replicate. Counts, fractions, and model diagnostics are written
alongside each per-type DE report.

For an abundance model, a candidate term is selected within a study only when
at least 95% of adult samples have valid values and the term varies; the fit
then uses complete cases across all selected terms. Smoking status is valid
only as canonical `yes` or `no`; all other values are treated as missing, and
the report labels the contrast `yes vs no` with `no` as the reference. This
95% coverage rule applies only to the abundance OLS models. DE optional
covariates instead enter per-study designs when they vary and are fitted on
complete cases, with no minimum coverage percentage.

The PyDESeq2 design and contrast fitting is reusable across categorical
covariates, but covariates are not registered in one shared list: configuration
validation, metadata preparation, reference-level selection, and report ordering
have explicit field lists. The abundance OLS report has its own candidate list
and categorical handling. When adding another covariate, update each path and
add a planted signal to the synthetic DE fixture/check so both workflows
exercise it.

# Output contracts

[Documentation](../README.md)

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
  diagnostic artifact with labels and method UMAP coordinates; PCA matrices and
  neighbor graphs are not included)
- `<outdir>/harmonization/integration_benchmark/integration_benchmark.json` (optional method
  and input provenance)
- `<outdir>/harmonization/integration_benchmark/methods/method_<name>.json` (optional
  per-method completion status)
- `<outdir>/harmonization/integration_benchmark/{executed.ipynb,report.html}` (optional
  benchmark report)
- `<outdir>/analysis/per_type/<cell-type-slug>/split.h5ad` (reusable raw-count split)
- `<outdir>/analysis/per_type/<cell-type-slug>/analysis/{analysis.h5ad,report.html,executed.ipynb}` (optional cell-type analysis and reports)
- `<outdir>/analysis/per_type/<cell-type-slug>/analysis/site-assets.json` and `site-assets/*.png` (named, checksummed figures and fraction summaries consumed by the atlas)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/diagnostics.json` (input, configuration, code, dependency and artifact fingerprints for prepared DE diagnostics; diagnostic tables coexist with the original fitted-model outputs)
- `<outdir>/index.html` (links to harmonization QC, integration benchmark, combined cross-type analysis, and available per-type residual, DE, and age-clock reports)
- `<outdir>/analysis/per_type/<cell-type-slug>/ageing_clock/` (optional clock models, metadata, and reports)
- `<outdir>/analysis/per_type/age_clock_all_types_summary.json` (optional all-type clock run summary)
- `<outdir>/analysis/per_type/age_clock_cross_type/` (optional donor-matched mean-clock
  predictions, component and combined metrics, pairwise residual correlations, and report)
- `<outdir>/analysis/per_type/cell_type_manifest.json` (optional downstream provenance and completeness contract)
- `<outdir>/analysis/combined/differential_expression/*.csv` (combined cross-cell-type DE results and manifest)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/...csv` (per-study and per-type combined covariate results)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/combined/age_bin_pearson_residuals.npz` and `age_bin_pearson_residuals_samples.csv` (sample × gene residuals and sample metadata from the age-bin model)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/age_model_diagnostics.csv` (sample support and raw-count inputs for the shared age-model diagnostic figure)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/reference_expression{,_per_study}.csv` and `reference_expression_metadata.json` (reference-bin CPM summaries and their rules)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/developmental/age_support_samples.csv` (all samples retained by the shared juvenile/adult age-slope model; the per-type report filters this to juvenile samples for its age-support plot)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/cell_type_result.json` (per-task fit record collected into the root manifest)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/{abundance_input.tsv,report.html,fraction_model_diagnostics.tsv}` (pseudobulk abundance input and per-type DE report outputs)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/residual_wgcna/` (residual WGCNA module assignments, eigengenes, soft-threshold and sample-connectivity diagnostics, covariate associations, and hub genes)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/residual_wgcna/{module_gsea.csv,cmv_associated_module_scores.csv,smoking_associated_module_scores.csv}` (preranked GO enrichment and sample scores for the strongest exposure-associated modules)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/residual_wgcna/{module_gsea_specific_summary.csv,module_gsea_specificity_audit.csv,module_gsea_provenance.json}` (specific-term summary, collapse audit, and observed GSEA/GO input provenance)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/go_enrichment/{merged_age_gsea.csv,merged_age_gsea_specific_summary.csv,merged_age_gsea_specificity_audit.csv,merged_age_gsea_provenance.json}` (merged-age GO enrichment, leading-edge summary, collapse audit, and observed input provenance)
- `<outdir>/analysis/combined/differential_expression/differential_expression.json` (DE status and cell-type result index)
- `<outdir>/analysis/combined/trajectory/{report.html,executed.ipynb}` (combined trajectory report)
- `<outdir>/analysis/combined/trajectory/*.csv` (shared-bin cross-type clusters and means, gene-level DE recurrence, and cross-type pattern concordance)
- `<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/developmental/age_25_vs_5.csv` (per-gene direct Wald contrast for model-predicted age 25 vs age 5, with DESeq2 adjusted p-values; the contrast holds other covariates fixed and uses the fitted group-specific linear age trends)
- `<outdir>/analysis/combined/trajectory/developmental_adult_age_effects.csv` (per-cell-type paired raw continuous-age slopes, mean normalized expression, and adjusted p-values for juvenile and adult contrasts from one shared all-age model; the combined report plots genes retained by both contrasts' independent expression filters)
- `<outdir>/analysis/combined/trajectory/developmental_adult_age_effect_support.csv` (sample and gene support for that comparison)
- `<outdir>/analysis/combined/trajectory/expression_atlas_*.csv` (study-balanced common-gene expression matrix, gene UMAP, depth support, and independent 3′/5′ and intronic-read sensitivity summaries)
- `<outdir>/analysis/combined/trajectory/*_residual_{clusters,pc_scores}.csv` (per-cell-type sample cluster assignments and PCA scores)
- `<outdir>/analysis/combined/trajectory/cross_cell_type_pearson_residuals.npz` and residual metadata/cluster CSVs (wide matrix, Leiden-primary assignments, K-means/HDBSCAN alternatives, marker contrasts, and cell-type feature composition for shared samples across full-trajectory-supported cell types)
- `<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/` (focused gene-profile HTML and executed notebook, `report_render.log` for notebook-render diagnostics, model coefficient and residual tables, mean–dispersion diagnostics, normalized sample counts, compressed VST matrix, and cross-cell-type age-bin trajectory/support tables when DE results are available)
- `<outdir>/run_manifest.json` (requested, selected, and skipped studies;
  hashes of selected source inputs, metadata dependencies, and annotation models;
  selected studies are the complete requested set for explicit study lists)

### Harmonized expression contract

Each harmonized H5AD has exactly one expression representation: `.X` contains
the configured unnormalized integer count matrix. `.raw` is `None` and
`.layers` is empty. Normalized matrices may be created transiently for
annotation or embeddings, but are never retained in the output. This prevents
an alternate normalized matrix from being mistaken for analysis input.
For CellTypist harmonization, the validated count object is written to the
output H5AD before normalization. Analysis then uses the in-memory `.X`; after
annotation, that normalized matrix is released and metadata and embeddings are
updated in place, preserving the saved count matrix.

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
labels, and method UMAP coordinates, but no expression matrix, PCA matrices,
or neighbor graphs. Per-method checkpoint H5ADs retain graph data in the
Nextflow work cache so aggregation can calculate sparse edge overlap. Its
`benchmark_*_aifi_l2_majority` columns are diagnostic; the retained per-study
`obs['aifi_l2_majority']` remains the downstream label. The notebook reports
graph edge overlap and AIFI-L2 concordance.
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
metadata and embeddings in place. This preserves the on-disk counts without an
in-memory raw-count copy.
PCA preparation uses the shared
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
artifacts. It has no DE manifest dependency. Each per-type DE report renders
from its fit result as an independent task; the combined cross-type report
renders separately from the collected manifest. Cell-type
analysis and reports publish under
`analysis/per_type/<cell-type-slug>/analysis/`; per-type DE results publish
under `analysis/per_type/<cell-type-slug>/differential_expression/`. The
collector writes the shared manifest under
`analysis/combined/differential_expression/`, and the cross-type report is
under `analysis/combined/trajectory/`. The single cell analysis report does
not require DE results, while `run-de` renders a separate DE report for each
cell type from the pseudobulk outputs. When the per-type merged age result is
available beside the single-cell report, the notebook adds a 3 × 3 UMAP panel
for up to five positive and four negative age-DE genes. It selects FDR-significant
genes (adjusted p < 0.05) by adjusted p-value, overlays per-cell log-normalized
expression, and draws low-expression cells before high-expression cells. The
overlay is descriptive; significance comes from the merged sample-level
pseudobulk model. The site release extracts PNG outputs from the executed
notebook, without per-cell expression values. Single-cell figures, including
per-cluster fraction plots, are stored as separate PNG assets referenced by
summary JSON. Export applies lossless level-9 recompression only when it saves
bytes and deduplicates identical optimized PNGs by SHA-256. The original
notebooks and copied HTML reports remain unchanged. The single-cell page places
the age-DE image alongside the sample-level cluster-fraction plot. See the
[static atlas data contract](../site/data-contract.md) for publication formats and
the distinction between linear DE results, per-type age-bin previews, and
clustered UMAP trajectories.

Each study is fitted once with its maximal available design: age and sex plus
each optional covariate with at least two observed values, plus `study_site`
when it varies within the study, and `log10(total_counts)` as a technical
covariate. Its coefficients use one complete-case sample set.
Combined fits target one covariate at a time and include all complete-case
samples from studies that recorded it, adjusting for `study_site`, age, sex,
and log10 total counts. The per-type DE report presents fit status, covariate
volcano plots, age-bin trajectories, residual clustering, residual WGCNA,
abundance models, and result tables, including the depth coefficient where
estimable. A separate diagnostic task executes the abundance OLS/bootstrap,
trajectory and residual clustering, WGCNA and enrichment calculations from
saved DE outputs. Rendering reads the resulting tables. Rerendering skips
PyDESeq2 fitting and reuses diagnostics when their input, configuration, code,
and dependency fingerprints match.
The site term is omitted when only one site remains. Fit metadata and CSVs
record the formula, included studies and sites, complete-case sample count, and
exclusions for missing design values.
Continuous design columns are centered within each fit before PyDESeq2 runs.
This changes only the intercept parameterization and prevents an artificial
age-zero or zero-depth intercept from destabilizing the optimizer.

The explicit depth covariate is used in addition to PyDESeq2's size-factor
normalization. It adjusts for the age-associated pseudobulk-depth pattern
observed in Naive CD8 T cells and is included as a technical covariate in every
per-study and combined fit, including the age-bin trajectory model.

The configured `l2_parent_l1` mapping is an inferred taxonomy with recorded
provenance, rather than the independent `aifi_l1_majority` predictions. It is
applied only to retained per-study L2 labels; diagnostic benchmark L2 calls are
not used in this hierarchy. Validate the mapping against the AIFI atlas
before using within-L1 fractions for inference.

## Workflow execution records

Each Make workflow writes `workflow_summary.json` at the output root and keeps
its trace and summary under `provenance/workflows/<target>-<run-id>/`. A `partial`
summary identifies ignored analysis/report failures while preserving successful
outputs. See [orchestration](../development/orchestration.md#failure-isolation).

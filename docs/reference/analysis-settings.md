# Analysis settings

Analysis commands validate their consumed settings before reading H5AD inputs or
fitting models. DE, cell-type analysis, gene profiles, and report preparation use
the shared `load_analysis_configuration` loader in `pbmc_pipeline.config`.
Age clocks reuse the DE sample-inclusion validator without requiring DE model
or trajectory settings. Required configuration documents declare `schema_version: 1`.

The analysis loader accepts the standalone pipeline JSON staged by Nextflow;
it does not open study registries, observation schemas, or source files. Core
configuration loading also checks study references against the registry. Both
paths reuse the same analysis validators. Validation returns the document
unchanged, preserving scientific defaults and configuration digests.

[Documentation](../README.md)

## Age-trajectory settings

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
| `residual_leiden_min_clusters` | 10 | Target cluster count; select the candidate resolution whose cluster count is closest to this value, with ties favoring fewer clusters. |
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
nearest-neighbor graph. The selected resolution is the configured candidate
whose cluster count is closest to `residual_leiden_min_clusters`, with ties
favoring fewer clusters. This avoids selecting an extreme overshoot when the
candidate resolutions produce a large jump in cluster count.
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

Each per-type DE diagnostic task runs residual WGCNA from the fitted age-bin
residual matrix. It publishes module assignments and eigengenes together with
soft-threshold, sample-connectivity, eigengene-covariate, hub-gene, CMV and
smoking checks and preranked GO Biological Process enrichment. GSEA uses the GO Biological Process 2023 library downloaded to
`input_data/go/`
and its output is controlled by GSEA's native FDR q-value within each module.
The compact summary collapses a parent term only if a same-direction descendant
has a leading-edge Jaccard overlap of at least 0.5; it retains a visible
representative for every collapsed term. The complete GSEA table, specificity
audit, and a provenance manifest containing resource hashes and parameters are
published alongside the report. CMV and smoking are checked because the residual model
already adjusts for age bins, sex, study site, and depth but omits those plausible
sources of remaining biological structure. These exploratory checks ask whether
signal remains; they do not assume an effect is present or adjust the residuals.
WGCNA color names are arbitrary module identifiers; grey is reserved for genes
not assigned to a module and is excluded from GSEA because it is not a coherent
module.

The merged continuous-age DE report uses the same downloaded GO Biological Process
collection. It ranks all genes tested in the merged model by signed Wald
statistic when the statistic is complete, otherwise by signed `−log10(p-value)`.
Its collapsible summary lists FDR-significant specific terms, the number of
leading-edge genes, and a preview of the first twelve ranked leading-edge genes;
the complete list is retained in the published CSV. Positive and negative NES
values correspond to genes increasing and decreasing with age, respectively.
GO Molecular Function and Cellular Component are not currently tested.

## Expression-atlas settings

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

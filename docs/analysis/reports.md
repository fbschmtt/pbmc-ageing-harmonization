# Reading and rerendering reports

[Documentation](../README.md)

## QC reports

Passing `--qc` executes the tracked `reports/qc_report.py` template after each
harmonized H5AD is written. Each study gets a separate directory:

```text
output/harmonization/harmonization_qc/<study>/
├── executed.ipynb
└── report.html
```

The report includes provenance, metadata completeness, cohort
composition, count-depth checks, UMAPs, canonical PBMC marker expression, and
study-label/AIFI-label comparisons. All plots use the complete harmonized study.
When stored UMAP coordinates and AIFI L2 annotations are available, the report
includes a standalone AIFI L2 plot with a transparent background. When study
labels and AIFI L2 annotations are available, it also includes a comparison
heatmap with a transparent background. The optional published color table and
the deterministic fallback are described in the
[AIFI model setup](../../aifi_models/README).

To generate a report for an existing harmonized output without rerunning the
study:

```bash
pbmc-qc \
  --input output/harmonization/harmonized/onek1k.test.h5ad \
  --run-report output/harmonization/provenance/onek1k.test.json \
  --output-dir output/harmonization/harmonization_qc/onek1k.test \
  --study onek1k
```

The merge workflow renders one combined document at
`output/harmonization/harmonization_qc/merged/`. It
covers pseudobulk composition, technical covariates, and a covariate-completeness
heatmap showing the fraction of pseudobulk profiles missing each covariate per study. String
metadata uses the canonical `not_provided` value for both source missingness and
explicitly unreported fields, while numeric metadata uses `NaN`. The pseudobulk outer
union retains all input counts. Shared-gene availability, hypothetical count
loss under the all-study intersection, and the highest-count excluded gene are
reported with the cross-cell-type expression atlas, where that intersection is
actually used. Cross-study UMAP and benchmark AIFI-L2 comparisons belong to the
separate integration benchmark; they are computed on the cells selected for
that benchmark and do not provide a full-data embedding if the configured cap
is exceeded.
The merged artifact retains per-study labels in `aifi_l2_majority` for all
downstream grouping and splitting. Run `make run-integration-benchmark` after
the core workflow to create a diagnostic H5AD with UMAPs and CellTypist labels
majority-voted on each method's graph. Per-method checkpoint files retain the
graphs in Nextflow's work cache for sparse edge-overlap calculations; the final
H5AD omits the graphs and PCA matrices. Scanorama is temporarily disabled in the active
benchmark configuration because its correction step can require large dense
temporary arrays. By default, the benchmark processes the
full merged H5AD unless it exceeds `integration_benchmark.max_cells`, currently
20,000,000. Above the cap, it selects cells uniformly at random without
replacement using `processing.random_seed` (42 by default); that selection is
shared by all methods. Set the cap in `config/pipeline.json` to change this
behavior. The report and JSON sidecar record source and selected cell counts,
fraction, strategy, and seed. The input is opened in backed mode, and only the
selected rows are materialized. The published artifact does not duplicate the
merged expression matrix. Graph edge overlap and per-study versus benchmark
AIFI-L2 label concordance are diagnostic comparisons, not acceptance
thresholds.
The benchmark runs each configured method as its own sequential Nextflow task.
Its checkpoint H5AD remains in Nextflow's work cache for `-resume`, while the
small status JSON is published under
`harmonization/integration_benchmark/methods/` as soon as the task completes.
`-resume` reuses completed method tasks after an interruption. The final
artifact and report combine successful methods; ordinary method errors are
recorded and do not prevent later methods from running. A partial report lists
failed or interrupted methods, and no HTML report is rendered if every method
fails. Add an implementation in `METHOD_ADAPTERS` and its configuration
validation in `pbmc_pipeline.config` when introducing another method.
Timestamped benchmark stage progress is captured in the Nextflow task log
(`.command.out` and `.command.err` in that task's `work/` directory), like
other process output. The terminal shows Nextflow task progress. Set
`INTEGRATION_BENCHMARK_LOG_LEVEL=DEBUG` for additional diagnostics. Production
and test targets use the standard Python image.
Its `report.html` has separate embedding-integration and type-calling sections,
including graph edge overlap and per-study-versus-benchmark AIFI-L2 concordance.
`output/run_manifest.json` records requested, selected, and skipped studies,
image references, configuration checksum, and checksums of the merged
deliverables. For Make-driven runs, `pipeline_revision` defaults to
`git describe --always --dirty`: a tag-based description when available or an
abbreviated commit ID, with `-dirty` for tracked changes. Untracked files are
not identified, and a dirty marker does not preserve the changes themselves.
`PIPELINE_REVISION` can override the value; direct Nextflow invocation defaults
to `unknown` unless `--pipeline_revision` is supplied. The core run manifest
identifies the core run; it does not prove that separately run DE analyses or
reports used the same revision. DE metadata records its own configuration
digest and PyDESeq2 version, while a site export has a separate export ID and
exporter revision.

### Report artifacts and web delivery

Each report run publishes a self-contained static `report.html`, which can be
opened directly in a browser or served from a static file host, plus an
`executed.ipynb` that retains code and outputs for technical reproduction. The
HTML report hides implementation cells and execution prompts, adds a compact
summary/navigation layer, and is the reader-facing analysis artifact; the
notebook is the audit artifact. The separate [static presentation site](../site/README.md#static-presentation-site)
uses exported structured results rather than parsing notebook HTML.

#### Report navigation

Rendered report navigation is generated from the final HTML headings. When
editing a report, change the headings in the template rather than maintaining a
separate table of contents. Give a heading an explicit ID when its link must
remain stable if its text changes.

## Single cell analysis reports

The optional downstream workflow opens the completed merged single-cell H5AD in
read-only backed mode and splits it by the retained per-study
`aifi_l2_majority` call. It materializes and writes one type at a time, keeping
the full sparse expression matrix on disk during splitting. AnnData still loads
the observation metadata into memory, and the largest individual type must fit
in memory while its split is written. Experimental benchmark labels are never
used for this split. For each type, `ANALYSE_CELL_TYPE`
creates the type-specific PCA/UMAP/Leiden embedding, marker tables, PC--age
correlation table, and `analysis.h5ad`. The PCA section shows each principal
component's individual share of variance explained, not cumulative variance.
Analysis writes the raw-count `analysis.h5ad` before normalization, then updates
its metadata and embeddings in place after analysis finishes. This avoids
holding a second full counts matrix in memory. Shared feature selection and PCA preparation live in
`pbmc_pipeline.embedding.prepare_sparse_pca`; per-type graph construction,
Harmony correction, clustering, and marker ranking are in
`pbmc_pipeline.cell_type.analyse_cell_type`.
`RENDER_CELL_TYPE_REPORT` then loads
that analysis artifact plus the primitive raw-count H5AD to create the executed
notebook and HTML. It also plots the retained per-study `aifi_l3_majority`
labels; no merged L3 prediction is created. Global neighbours are deliberately
not copied into the split artifacts.
During report execution, the notebook also computes an unintegrated UMAP from
the retained native PC scores for comparison; this embedding is not written
back to `analysis.h5ad`.

To analyze a subset, pass comma-separated slugs to either the production or test
target, for example:

```bash
make run-cell-type-analysis CELL_TYPES=cd14-monocyte,naive-cd4-t-cell
make run-cell-type-analysis-test CELL_TYPES=cd14-monocyte,naive-cd4-t-cell
```

Every requested slug must occur in the selected input. Unknown or absent types
stop the workflow before per-type analysis starts, with the available slugs in
the error message; a mixed valid/invalid request also fails. Omitting `CELL_TYPES`
analyzes all available types. Direct Nextflow `--cell_types` selection also works
with existing splits or split/analysis pairs for rerendering. Rerender inputs
must have matching, unique slugs; unmatched paths cannot silently disappear in
a join. Selection cannot be combined with `--split_only`.

Sample-level abundance and fraction views are generated with the pseudobulk
DE reports, where cell counts are summed from the merged pseudobulk object.
The cell-type report retains local-cluster composition views, includes zero
fractions for eligible samples with no cells in a cluster, and excludes
samples below `cell_type_analysis.min_cluster_fraction_denominator_cells`.
Cluster-fraction plots show a pooled study-unadjusted linear trend and separate
per-study linear trends when each has at least two distinct ages.
The configured `cell_type_analysis.l2_parent_l1` map supplies the parent label
in `report.json` and the cell-type manifest, and defines the parent denominator
for DE abundance. The split and analysis H5ADs do not carry a derived parent
column or per-sample abundance totals.

```bash
make run-cell-type-analysis
```

This downstream workflow requires an existing merged single-cell H5AD at
`<outdir>/harmonization/merged/single_cell_merged.h5ad`. Override that location with
`MERGED_INPUT=/path/to/single_cell_merged.h5ad`.
When a per-type analysis or report fails, the merged-input workflow continues
with other types. The report index lists successful types and names failed
types; the cell-type manifest records whether each failure occurred during
analysis or report rendering. Nextflow logs contain the underlying error.

For a complete production run, `make run-all` executes the core workflow with
single-cell merging enabled, then runs differential expression, detailed
profiles for the ten lowest-adjusted-p-value age-DE genes in each cell type,
and cell-type analysis/reporting in order. The profiles are written under
`<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/`. Run
`make profile-top-age-genes TOP_AGE_GENE_COUNT=5` to change the number or to
rerun profiling from existing merged DE results and pseudobulk data. The
workflow is resumable. It uses the configured `STUDIES`, `OUTDIR`, and
`WORK_DIR` values. For test output, run `make run-test`
first, then `make run-cell-type-analysis-test` to consume the test merge under
`output/test/`. All test targets publish only beneath `output/test/`.

The central HTML report index is written to `<outdir>/index.html`. It links
to the merged and per-study harmonization QC reports, integration benchmark,
combined cross-cell-type analysis, and per-type residual, DE, and age-clock
reports. When present, the run manifest limits study QC links to its selected
studies, and the cell-type manifest defines the per-type list. Without a
cell-type manifest, the index discovers types from available report files.
Age-clock entries include the merged method comparison and each method's report,
based on the clock manifest. Core, integration-benchmark, DE, cell-type, and
age-clock Make targets refresh it after publishing. To rebuild it manually
after moving or removing reports, run
`make render-report-index OUTDIR=<outdir>`.

The raw-count splits are published under `<outdir>/analysis/per_type/<cell-type-slug>/split.h5ad` so one
report can be rerun without re-splitting the merged H5AD:

```bash
make run-cell-type CELL_TYPE=cd14-monocyte
make run-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type CELL_TYPE=cd14-monocyte
make render-cell-type-test CELL_TYPE=cd14-monocyte
make render-cell-type-reports
make render-cell-type-reports-test
```

Published outputs and the Nextflow work cache are managed directly in the
filesystem. Remove only the directories you intend to discard; production and
test artifacts live under separate output roots, while Nextflow's cache is
stored separately under the configured work directory.

The targeted test command deliberately does not build a core test merge or all
reports; it requires the corresponding split to already exist. It runs the
analysis and rendering stages for that one type. Create only the reusable splits
from an existing merge with `make split-cell-types-test`, or run `make run-test`
followed by `make run-cell-type-analysis-test` for all types. A targeted rerun
refreshes that type directory but not `cell_type_manifest.json`; rerun the
all-type workflow before treating the published report set as a new
manifest-backed release. `cd14-monocyte` is the default targeted test type
because it is well represented in the bundled multi-study fixture.

For a report-only edit, `make render-cell-type[-test]` rerenders one type, while
`make render-cell-type-reports[-test]` rerenders every type. These targets
require the published splits and `analysis.h5ad` artifacts. They schedule only
report execution and overwrite the executed notebooks and HTML without
rerunning local PCA, Harmony integration, Leiden clustering, or marker ranking.
The report execution does calculate its comparison UMAP from saved native PCs.

Use `make check-cell-type-test-prerequisites CELL_TYPE=<slug>` before a
targeted test rerun when unsure whether its input exists. It prints the minimal
next command: create splits from the existing test merge, or create a fresh core
test merge and reports when that merge is absent. The production equivalent is
`make check-cell-type-prerequisites CELL_TYPE=<slug>`.

Each published type directory under `output/analysis/per_type/` contains its
raw-count `split.h5ad`, `analysis/analysis.h5ad`, HTML report, executed notebook,
and report tables. `cell_type_manifest.json` records the exact merged input,
configuration, reference-model checksums, analysis specification, completed
type reports, and any failed type and stage. A failed type does not prevent
other reports from being published: the manifest records partial status and
the workflow exits with an error. Each type marked complete must have all
required deliverables. Site export includes only types marked complete, so a
stale report for a failed type is not presented as part of the current release.
PC--age correlations in these reports are explicitly descriptive at the cell
level; use sample/subject-level pseudobulks or mixed models for inference.
The local neighbour graph, UMAP, and Leiden clustering use Harmony correction
by study; native (uncorrected) PCs remain the basis of the current PC--age
diagnostic. This separation is deliberate but should be revisited before any
inferential downstream use, since integration can also alter age-associated
structure. Each report now shows one study-coloured UMAP from native PCs before
the existing Harmony-integrated UMAP plots. V(D)J genes matching
`processing.exclude_vdj_regex` are excluded from the local HVG/PCA input; the
report records the count excluded. HVGs use Scanpy's configured mean-dispersion
selection on normalized, log-transformed expression.

Harmonization and per-cell-type analysis use the same PCA preparation function,
`pbmc_pipeline.embedding.prepare_sparse_pca`; the selected HVG thresholds and
maximum PC count are set in `config/pipeline.json` under `processing.hvg` and
`processing.neighbors`. It selects mean-dispersion HVGs using the configured
minimum mean and dispersion cutoffs, without a maximum mean cutoff, and excludes
genes matching `processing.exclude_vdj_regex`
(case-insensitively), falls back to all non-V(D)J genes if fewer than two HVGs
remain, casts the selected matrix to float32, scales without centering to retain
sparsity, and runs ARPACK PCA with implicit centering. The PCA component count
is capped at `processing.neighbors.n_pcs` and at the dimensions available in the
selected matrix, so analysis does not compute PCs that the configured neighbor
graph will not use. Stage-specific embedding orchestration is in
`pbmc_pipeline.harmonize._compute_embedding` and
`pbmc_pipeline.cell_type.analyse_cell_type`; Harmony correction is shared via
`pbmc_pipeline.harmonize._harmony_integrate`.

Local embeddings are exploratory. The current eligibility threshold is 50 total
cells and does not require balanced study representation; small or highly
imbalanced types can yield visually structured UMAPs even after Harmony. In
particular, the configured 50-neighbour graph is relatively dense for rare
types. Inspect study-coloured UMAPs and per-study cluster composition before
interpreting local clusters as biological states.

Nextflow's `work/` directory is the cache for task intermediates. Published
outputs are copied to `output/`, including final reports, reusable cell-type
split H5ADs, and the derived per-type analysis H5ADs; do not use paths inside
`work/` as report inputs.

Workflow targets build their required content-tagged images automatically; use
`make images` only to build them without a run. For a production run whose
repository and task filesystem are separate, keep the repository readable by
Docker and choose external work/output locations:

```bash
make run NF_PROFILE=docker STUDIES=all \
  WORK_DIR=/path/to/run/work OUTDIR=/path/to/run/output
```

Nextflow stages declared inputs into the task directory. Keep the repository
readable to Nextflow and use a `WORK_DIR` that supports symbolic links.

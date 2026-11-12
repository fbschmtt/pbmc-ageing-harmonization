# Pseudobulk differential expression

[Documentation](../README.md)

Each cell-type DE run also writes `reference_expression_per_study.csv`,
`reference_expression.csv`, and `reference_expression_metadata.json` under its
`differential_expression/` directory. For each eligible pseudobulk in the configured
reference age bin, CPM uses the sum of measured gene counts as its denominator.
Within each study, mean `log(CPM + 1)` is calculated for each measured gene; the
combined summary averages those study means equally and transforms back to CPM.
The CSVs report sample counts and, for the combined summary, contributing study
counts. Genes unavailable in a study are excluded from that study's summary.
These are descriptive, unadjusted summaries of the available cohort, not fitted
or sex-standardized expression. The reference cohort applies the minimum-cell and
configured subject-exclusion rules but does not require sex or study-site metadata.
Test mode bypasses the minimum-cell threshold as recorded in the metadata.
For older DE outputs without these CSVs, run `make run-de` using the existing
pseudobulk merge; rerendering reports or exporting the site does not calculate
them. The site's trajectory details show the combined CPM and support counts
when the CSV is available.

### Gene universes

The workflow intentionally uses different gene universes for different tasks:

| Context | Gene universe | Purpose |
|---|---|---|
| Merged pseudobulk | Outer union across studies | Retains each study's observed genes for per-study models. |
| Per-study DE | Genes available in that study | Fits only observed counts for that study. |
| Combined DE | Intersection across studies included in that fit | Prevents study-absent outer-join zeros from entering the model. |
| Study-balanced expression atlas | Intersection across all studies | Makes its equally study-weighted expression values comparable. |

### Covariate policy across analyses

The general rule is to use the broadest adjustment set that the analysis
population can support. For a covariate-focused estimate, use the largest
study set with usable observations and enough variation to estimate the target
contrast. Missing or unclear values are excluded from that model; they are not
treated as another biological level. Each model records its design, included
studies, and sample counts so the actual support can be checked.

| Analysis | Covariate use and sample selection |
|---|---|
| Per-study DE | Uses age, sex, library depth, and study site when estimable, plus each optional covariate (BMI, CMV, smoking) that varies in that study. Samples missing a term in that model are excluded. This is the maximal varying optional-covariate set within each study. |
| Combined primary age DE | Uses the configured core terms: study site when estimable, age, sex, and library depth. It does not add every optional covariate, which would unnecessarily restrict the shared age-analysis cohort when those fields are sparsely recorded. |
| Combined covariate-focused DE | Fits one optional covariate at a time on top of the core terms. It uses complete cases and retains studies with usable observations; smoking is explicitly `yes vs no` with `no` as the reference. Unlike abundance models, DE has no 95% coverage rule. The current fit checks variation and estimability across the combined design, but does not impose a minimum sample count in each study or exposure group; inspect the recorded study counts before interpreting a focused result. |
| Juvenile versus adult age effects | Fits one model across all eligible ages with shared study-site intercepts, an age-group intercept at the juvenile/adult boundary, and separate within-group age slopes. It exports the adult age slope, the juvenile age slope, and a direct Wald contrast for predicted age 25 vs age 5, with gene-wise adjusted p-values for each contrast. The predicted contrast holds the other model covariates fixed and uses the fitted piecewise-linear trend; no study-specific slopes are fitted. Every otherwise eligible juvenile sample is retained. A rank-deficient shared design is recorded as non-estimable. The report shows juvenile-only age support, the adult-versus-juvenile slope comparison, and the age-25-versus-age-5 contrast against the juvenile slope, followed by the juvenile volcano plot and gene table. |
| Age-trajectory DE | Replaces continuous age with age bins and retains the other core terms. It uses the configured age and bin-support filters; optional covariates are not added. |
| Abundance changes | Fits adult-only, study-specific OLS models. Within each study, an optional covariate is included only if at least 95% of adult samples have a valid value and it varies. The model then uses complete cases across its selected terms. Smoking accepts only `yes` and `no`; other values are missing. These are study-specific estimates, not a pooled abundance effect. |
| Age-bin abundance | Fits separate OLS fraction models per study and denominator, with each study's 20–30 reference contrast set to zero. The Composition tab anchors those contrasts at each study's adjusted reference fraction to show absolute estimated fractions. Uses ten-year bins below age 90, the configured per-study minimum support, and available varying study site, sex, BMI, CMV and smoking covariates with the same 95% rule and complete cases. |
| Gene profiles | The PyDESeq2 age model uses the combined primary design. The exploratory age-bin model uses the same core adjustments; the per-study-slope and mixed models adjust for sex and library depth. These profile models do not add BMI, CMV, or smoking. |
| Expression atlas and residual WGCNA | Atlas technology contrasts are descriptive, not covariate-adjusted DE. WGCNA's metadata associations are post-fit associations with module eigengenes, not covariate adjustment in the gene-level fit. |

For a future covariate-focused analysis with a within-study support rule, define
the minimum sample and group counts explicitly before fitting. The current
combined DE implementation uses complete cases and global design
estimability, so it can retain a study even when that study does not contain
both levels of a categorical contrast.

The separate `differential_expression.nf` workflow consumes the merged
sample-level pseudobulk H5AD. It first lists retained per-study
`aifi_l2_majority` labels, then runs one independent task per label. Each task
loads the same small merged pseudobulk H5AD and subsets its label in memory; it
does not materialize a per-type pseudobulk H5AD. Per-study models adjust for
`study_site` when multiple sites remain in that study, along with age, sex, and
available optional covariates. The merged age-and-sex model uses
`~ study_site + age + sex + log10_total_counts` across studies. `log10(total_counts)`
from each pseudobulk is included as a technical covariate in all DE models and
its coefficient is shown in the combined-model result tables. When
present and variable, BMI, CMV, and smoking status each receive separate
combined fits adjusted for study site, age, sex, and log10 total counts.
Before PyDESeq2 fitting, continuous design columns are centered within each
fit. This changes the intercept parameterization, while preserving the
interpretation of reported slopes and contrasts and improving optimizer stability.
For each merged age contrast, the report also runs preranked enrichment against
the downloaded GO Biological Process 2023 collection. It ranks all tested genes by
their signed Wald statistic when complete; otherwise it uses signed
`−log10(p-value)`. Positive normalized enrichment scores indicate terms among
genes increasing with age and negative scores terms among genes decreasing with
age. The collapsible summary uses the native GSEA FDR q-value and shows each
term's leading-edge gene count plus the first twelve ranked leading-edge genes;
the complete untruncated list remains in the CSV.
DE fits have no minimum covariate-coverage percentage: an optional term enters
a per-study design when it varies, and each fit uses complete cases for its
own design. Smoking status is accepted only as `yes` or `no`; other labels are
treated as missing. Its categorical contrast is `yes vs no`, with `no` as the
reference.
DESeq2 estimates size factors to normalize library depth, and this additional
depth term adjusts for the age-associated pseudobulk-depth pattern observed in
Naive CD8 T cells, in addition to size-factor normalization. It is included in
the per-study and combined DE models.
Each per-type DE task writes
`analysis/per_type/<type>/differential_expression/abundance_input.tsv` and its
report includes sample-level abundance plots and adult-only fraction models.
The separate `ANALYSE_CELL_TYPE_DE_DIAGNOSTICS` task prepares abundance models,
trajectory/residual clustering, WGCNA and enrichment before report rendering.
It publishes `diagnostics.json` plus the existing diagnostic tables. The manifest
records input and artifact checksums, configuration, analysis code and dependency
versions; unchanged diagnostics can be reused when rerendering. The notebook
reads the prepared tables and does not fit these models or rewrite their results.
Direct `pbmc-differential-expression-report --result-dir ...` invocations require
`--analysis-dir ...` from `pbmc-de-diagnostics --result-dir ... --output-dir ... --config ...`.

Pseudobulk `n_cells` values are summed within `study × sample` for the
retained-PBMC denominator and within `study × sample × configured L1 parent`
for the parent denominator. The abundance input includes every observed sample
for each reported type, with zero cells for sample/type combinations absent
from the pseudobulk table. Each fraction view applies its own denominator
cutoff (currently 10 cells); the per-type DE pseudobulk cell filter does not
exclude zero-count abundance observations. The configured L2-to-L1 map remains
in `cell_type_analysis.l2_parent_l1`; it is applied to retained L2 labels, not
the independent `aifi_l1_majority` classifier output. The existing OLS,
bootstrap, and plot implementation moved from the cell-type report section.
Denominator and bootstrap settings are configured under
`differential_expression.abundance`; diagnostics are published alongside each
per-type DE report.

The report also fits separate study-specific age-bin OLS fraction models for
both retained-PBMC and L1-parent denominators. Each supported bin is contrasted
with the same study's 20–30 bin; the fitted coefficients are aligned at zero
and are not absolute fractions. In the Composition tab, these contrasts are
offset by each study's model-adjusted reference fraction so the trajectory
points show absolute estimated fractions. Its point-based y-axis clips longer
contrast intervals and marks truncated ends with arrows. The minimum per bin is the configured
`study_comparison_minimum_samples_per_bin`; a supported reference and at least
one other bin are required. Intervals are pointwise 95% OLS intervals and
reported p-values are unadjusted across bins. Missing reference bins, insufficient
support, rank-deficient models, and legacy inputs without study-site metadata
are recorded as skipped in `fraction_age_bins_metadata.json` rather than
imputed. `fraction_age_bins.tsv` contains the estimates. Older DE outputs must
be regenerated to include study-site metadata before rerendering these models.
Fraction models retain zero-cell samples and therefore do not use cell-type
pseudobulk sequencing depth, which does not exist for zero-cell samples.
The sample selection follows the existing abundance analysis rather than the
DE-specific subject and cell-type pseudobulk filters.

The diagnostics stage fits adult-only (age ≥20), study-specific OLS models. A candidate
covariate is included only when at least 95% of that study's adult samples have
valid values and it varies; the fit then uses complete cases across all
included covariates. Smoking status is accepted only as `yes` or `no`; other
labels are treated as missing, and its contrast is `yes vs no` with `no` as
the reference. This 95% rule applies to the abundance models, not the DE fits.
The site term is dropped when only one site remains. A final collector writes
the result directories and manifest.
Per-type DE reports show model status and result tables, with up to 25 rows per
fit sorted by adjusted p-value. The gene-universe table above is the canonical
definition of which genes enter each fit. The parallel type-fitting task is
configured for one CPU, which is passed directly to PyDESeq2; all other CPU,
time, and memory requests use the executor defaults.
In addition to the continuous-age fit, each cell type gets an exploratory,
site- and covariate-adjusted age-trajectory model. Per-type DE reports show
model status and result tables for available per-study and combined fits.
`run-de` also renders the combined cross-cell-type analysis report. It
summarizes shared age trajectories, cross-type recurrence, and Pearson-residual
structure, and includes a descriptive, study-balanced expression atlas with
gene availability, count-loss diagnostics, expression maps, and technical-
label contrasts. It also reports unique sample counts per technical batch
within each study site to show batch-size imbalance and batch coverage. These
are exploratory outputs, not covariate-adjusted
estimates of absolute RNA abundance. Per-type DE results are stored under each type's
`analysis/per_type/<cell-type-slug>/differential_expression/`; combined DE and
cross-type artifacts are under `analysis/combined/`. See
[analysis settings](../reference/analysis-settings.md#age-trajectory-settings) for the model,
support rules, clustering, output tables, and current configurable settings.

The single cell analysis reports are generated separately by
`run-cell-type-analysis` and contain local embeddings, markers, coverage, and
fraction views. They do not consume the DE manifest. The residual WGCNA helper
is run by each per-type DE report from the fitted age-bin residual matrix. The
report includes soft-threshold and sample-connectivity diagnostics, eigengene
structure and covariate associations, hub genes, CMV evidence, and
Module GSEA ranks all selected network
genes by correlation with each module eigengene and tests GO Biological Process
2023 terms from the downloaded `input_data/go/` library. It publishes complete results as
`module_gsea.csv`, a hierarchy-aware specific-term summary, and a provenance
manifest recording the library and ontology files actually used. The summary uses the
native GSEA FDR q-value within each module and collapses a GO parent only when
a same-direction descendant has at least 0.5 leading-edge Jaccard overlap.
The strongest CMV- and smoking-associated modules are
also plotted on the residual UMAP and by observed yes/no/not-provided status.
Binary covariate heatmap values are
signed point-biserial correlations; multi-level categorical values are η², the
fraction of eigengene variance explained by category membership. These are
descriptive post-fit associations, not adjusted effects. CMV and smoking are
checked because the residual model adjusts for age bins, sex, study site, and
depth but omits these two plausible sources of biological residual structure.
The checks ask whether omitted traits remain associated with modules; they do
not assume that such effects are present. WGCNA color names are module IDs,
not biological classes. Grey is the special unassigned-gene label, so it is
excluded from module GSEA because those genes do not form a coherent module.
Module colors are assigned independently when a network is fit, so color names
are not stable across fits.
The tables are published beside the report under `residual_wgcna/`. To keep large networks readable, module-indexed
figures show at most the 25 largest modules and state when that limit applies;
the complete module-level tables remain available in that directory.

The current enrichment scope is GO Biological Process only. Molecular Function
and Cellular Component are not tested or combined into the reported FDR family.
The reusable runner accepts a caller-provided gene-set collection when a future
analysis needs one, and records that collection's content hash in its output
provenance.

Samples must be age 20 or older and have at least 10 cells in that
sample × cell-type pseudobulk. Samples missing required metadata or covariates,
or with zero counts over the fit's gene universe, are also omitted. Optional
covariate missingness receives additional handling within each study: a
per-study model includes an optional covariate only when it varies in that
study, and samples missing a covariate are omitted from fits that include it.
Combined optional-covariate models use complete cases from studies that recorded
the covariate. All sample exclusions and model fits that cannot estimate the
configured covariates are recorded in the JSON manifest. Combined results and
the manifest are published beneath
`<outdir>/analysis/combined/differential_expression/`; per-cell-type results
and task fit records are under
`<outdir>/analysis/per_type/<cell-type-slug>/differential_expression/`.

The DE manifest's `results_by_cell_type` index lists the per-study and combined
result paths for each covariate (or an empty result set) for each cell-type
slug. `run-de` uses it to render separate per-type DE reports. The cell-type
workflow does not read this manifest: `make render-cell-type-reports` rerenders
single cell analysis reports from published splits and analysis artifacts
without repeating analysis. That renderer writes single-cell artifacts beneath
`<outdir>/analysis/per_type/<cell-type-slug>/analysis/`; it never publishes
report files into the DE directory. Named Nextflow subworkflows group splitting
plus analysis, analysis of an existing split, and report rendering while
retaining the separate single-cell and pseudobulk entry points.

## Model status and partial results

Each cell type's `differential_expression/run_metadata.json` records the outcome
of every requested fit. Model statuses differ from Nextflow task statuses:

| Model status | Meaning |
|---|---|
| `complete` | The fit produced result tables for its requested contrasts. |
| `skipped` | The available data cannot support the model, for example after inclusion filters or because the design is rank deficient. The record gives a `reason`. |
| `failed` | Fitting raised an error, such as `LinAlgError: Singular matrix`. The record gives an `error`; it is not converted into an eligibility skip. |

A per-type task can complete successfully after recording one failed fit and
continuing its other models. Consequently, Nextflow's failed-task count is not
the number of failed DE models. Inspect both levels with:

```bash
make workflow-failures OUTDIR=output
```

The collector writes a manifest with model outcomes in its work directory.
By default, it then fails if any model failed or no model completed. Sibling
analyses and available per-type reports continue, and the workflow records
partial status; combined reporting cannot run without a successful collector.
Independent report-rendering failures also make the workflow partial.

To allow combined collection and reporting when at least one fit completes:

```bash
make run-de ALLOW_PARTIAL_RESULTS=true
```

The same option applies to `make run-all` and `make run-analysis`. Failed models
retain their errors, successful results are published, and both the DE manifest
and workflow summary remain `partial`. This option does not change inclusion
filters, model designs, or fitting rules. Collection still fails when no model
completes. Strict verification rejects partial results even when partial
collection is allowed. See [failure inspection](../development/orchestration.md#failure-isolation)
for archived summaries and Nextflow task-log commands.

## Per-type abundance in DE reports

Each per-type DE task writes an `abundance_input.tsv` with one row for every
sample observed in the merged pseudobulk input, including zero-count
sample/type combinations. Pseudobulk `n_cells` supplies the type numerator and
the retained-PBMC and parent-compartment denominators. These reports therefore
do not require per-cell abundance columns in split H5ADs. The per-type DE
inclusion filter does not exclude zero-count abundance rows.

The diagnostics stage fits adult-only, study-specific OLS models with the existing age,
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

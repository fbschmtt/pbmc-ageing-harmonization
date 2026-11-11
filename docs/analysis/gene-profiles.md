# Focused gene profiling

[Documentation](../README.md)

`make profile-gene CELL_TYPE=<label-or-slug> GENE=<gene-id>` fits the configured
merged DESeq2 model and a naive normalization linear model for one cell type
from the pseudobulk H5AD, then renders a focused report for one gene inside
the Python Docker image through Nextflow. Model preparation, gene calculations,
and notebook rendering are separate tasks. For example:

```bash
make profile-gene CELL_TYPE=cd14-monocyte GENE=IL7R
```

The executed notebook and HTML report give each model its own top-level
heading. The PyDESeq2 section includes the gene's position in the fitted
mean–dispersion plot, its age coefficient and adjusted p-value, sample-level
size-factor normalized counts, and Pearson residuals against age and
`log10(total_counts)`, each faceted by study and colored by study site. The
mean–dispersion plot marks the selected gene's gene-wise dispersion before
shrinkage alongside the final dispersion used in the model. Correlated-gene
scatter plots color points by study. The
residual section shows a table with the selected gene followed by its ten
strongest residual correlates, including each gene's mean normalized count and
final dispersion, and plots the four strongest correlates. A separate Naive
normalization section documents geometric-mean total-count scaling and the
`log(normalized count + 1)` transform once, shows the normalized-count
histogram, and compares naive expression with PyDESeq2 VST. Continuous-age
models use a geometric mean from their included samples; the age-bin model
calculates its geometric mean across all samples in that cell type before
applying model filters, so lower-age samples contribute to its scale. The naive
linear age model uses the PyDESeq2 design covariates and reports its coefficient
table and residual plots against age and library depth. The age-bin model uses
the configured decade width, includes samples below age 20, adjusts for the
other varying design terms, and plots standardized fitted levels alongside the
sample expression values, with residual plots against age and library depth.

A further model allows a separate linear age slope and intercept per study,
adjusting for sex and library depth. `age * study` estimates study-specific
intercept contrasts and age-by-study slope differences, without pooling; its
residuals are also plotted against age and library depth. The mixed model keeps
study intercepts fixed and estimates random age slopes around a shared slope.
It reports fixed effects, random-slope and residual standard deviations, fit
convergence, and partially pooled conditional slopes for each study, plus
residual plots against age and library depth. Treat slope-variation estimates
as exploratory when only a few studies are available. Sample-level expression
points use 0.3 opacity so overlapping samples remain visible. The displayed
geometric means use scientific notation. The HTML table
of contents includes each implemented model heading, deferred `limma-voom` and
R DESeq2 headings marked `TODO: implement R execution`, and rendered
subheadings. The
target also writes a compressed, full fit-gene-set
VST matrix and the included sample metadata under
`<outdir>/analysis/per_gene/<cell-type-slug>/<gene-slug>/`. VST values are
transformed data for downstream exploration and modeling; differential
expression inference should continue to use raw integer counts.

The focused profile refits the combined primary age DESeq2 model from the
pseudobulk input. Nextflow prepares that fit once per cell type and stages it as
an immutable input to the selected gene
profiles of the same cell type. It does not currently reuse a fitted model
from `run-de`; that workflow publishes coefficient tables and fit metadata,
not a reusable fitted PyDESeq2 object. A shared fit artifact would be cleaner
for the matching combined primary design, provided it records and matches the
exact included samples, gene universe, design, and software version. Per-study
fits and covariate-focused fits have different samples or designs and must not
be substituted for that age model.

When completed `run-de` results are available, each focused profile also shows
the selected gene's existing age-bin trajectories across cell types. These are
the adjusted log2 fold changes relative to each cell type's own reference age
bin, so they describe within-type trajectories and do not compare expression
abundance between cell types. The report publishes the displayed rows and a
per-cell-type availability table beside the profile.

To smoke-test the complete fit and report with planted synthetic expression,
run `make run-gene-profile-synthetic-test`. It writes the pseudobulk fixture to
`output/test/analysis/combined/gene_profile_synthetic/` and the checked report artifacts to
`output/test/analysis/per_gene/cd14-monocyte/synth-age-marker-1/`. This standalone check is separate from pytest because it fits DESeq2 and executes
a notebook. The composed `make run-all-test` checks the same planted gene using
the synthetic input shared with DE and clocks.

`make profile-top-age-genes` applies the same detailed report to the top age-DE
genes from each cell type's merged model results. Ranking uses finite adjusted
p-values in `<outdir>/analysis/per_type/<cell-type>/differential_expression/merged.csv`; it
profiles up to ten genes per type by default, using the merged pseudobulk H5AD.
The report paths follow the same `analysis/per_gene/<cell-type>/<gene>/` layout as
the single-gene target.

See [workflow testing](../development/testing.md#de-test-fixtures) for synthetic DE coverage and test-mode inclusion rules.

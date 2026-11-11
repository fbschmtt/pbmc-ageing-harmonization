# Age clocks

[Documentation](../README.md)

`make run-age-clock` fits sample-weighted, study-balanced, within-sample rank, within-study
gene z-score, age-aware study-offset, and age-and-sex-aware study-offset elastic-net clocks
for one type selected with `CELL_TYPE` (default `cd14-monocyte`; label or slug). Set
`AGE_CLOCK_METHODS=<comma-separated IDs>` to fit a subset; omit it to fit all eligible
methods.

The z-score clock is the label-free naïve study normalization: it centers and scales every
gene within the available target-study cohort, so it cannot score a single sample.

The age-aware clock fits per-gene age trends and additive training-study offsets before
regression, then uses `AGE_CLOCK_CALIBRATION_SAMPLES` known-age target donors (default 10)
to estimate a target-study offset before predicting the remaining donors. Its calibration
donors are excluded from CV and external scores; it is a labeled target-site protocol, not
calibration-free deployment.

The `age_sex_aware_study_offset_elastic_net` method additionally uses sex to estimate the
training and target study offsets; sex is not an ElasticNet predictor, and sex effects
remain in its expression features. Its reference fit includes a sex main effect but assumes
a shared age slope across sexes. It requires canonical female/male labels in the external
samples and both sexes in training and every leave-one-study-out training fold. If all
methods are requested and this metadata or overlap is unavailable, the method is recorded
under `skipped_methods` in the clock manifest while other methods continue; explicitly
requesting this method raises an error.

The age-aware method assumes stable age effects, additive site shifts, and reasonably
comparable calibration and evaluation populations. The standard
`age_aware_study_offset_elastic_net` report also shows an exploratory post-hoc linear age
mapping fit from leave-one-calibration-donor-out calibration predictions; it is reported
separately from the original external predictions and metrics. Both age-aware methods use
pooled training-derived elastic-net feature scaling, not target-study SD scaling.

Training and held-out external evaluation apply the shared adult inclusion rule (age ≥20);
juvenile pseudobulks are excluded and counted in the inclusion diagnostics. The low-level
raw-pseudobulk prediction helper can technically score juveniles, but those predictions are
outside this adult clock's validated domain.

CV candidates run with `AGE_CLOCK_N_JOBS=20` by default; lower it to limit CPU use. Use
`make run-age-clocks` to fit every cell type present in the pseudobulk. It logs and skips
types with no eligible samples in any configured training or external study; an unexpected
type-specific fit failure is recorded in
`analysis/per_type/age_clock_all_types_summary.json` and the workflow trace, while
successful types are still reported. The Nextflow run continues and its summary
records `partial`; strict verification rejects partial runs after completion.
Fitting, per-type reporting, and cross-type aggregation are separate tasks.

Defaults to holding out `nehar_belaid26,perez22`; these relatively small studies retain as
much training data as possible while spanning a broad validation age range. Override with
any comma-separated list via `AGE_CLOCK_EXTERNAL_STUDIES=study_a,study_b`. Per-method
reports and a merged comparison go to `analysis/per_type/<slug>/ageing_clock/`. The
all-cell-type run also writes an exploratory cross-cell-type clock analysis to
`analysis/per_type/age_clock_cross_type/`. For each clock method shared across at least two
cell types, it averages the per-type predictions for each donor and compares the combined
performance with its component clocks on matched donors. It also records pairwise
correlations between per-type prediction residuals for training and held-out donors in
`residual_correlations.csv`, and renders a report at
`analysis/per_type/age_clock_cross_type/report.html`. Donors need predictions from at least
two cell types, with matching study and age metadata. These results assess whether combining
cell-type-specific expression improves age prediction; they do not combine different clock
methods. The site export also derives the matched-donor mean from available per-type
predictions when an older analysis output does not yet contain the aggregate summary and
the component clocks share the same input, configuration, and study splits.
To build the aggregate and diagnostics from already fitted per-type clocks without
refitting them, run `make aggregate-age-clocks`; then run `make site-export` to publish
the selector option. The aggregate report remains in the analysis outputs;
it is not copied into the site's HTML report bundle.

Smoke-test with `make run-age-clock-synthetic-test` and check that expected output files
exist with `make check-age-clock-synthetic-test`. The composed `make run-all-test`
and `make verify` exercise the same clock workflow against their shared synthetic input.

## Synthetic age-clock experiments

The synthetic pseudobulk generator is deliberately compact: its single
implementation, [`src/pbmc_pipeline/synthetic_de.py`](../../src/pbmc_pipeline/synthetic_de.py),
defines the studies, assay-availability flags, planted effects, and count
sampling in `create_synthetic_pseudobulk`. For one-off clock diagnostics,
modify that source directly rather than adding a scenario framework; record
the command and seed, and revert temporary source edits afterwards.

The default fixture is a software smoke test, not evidence that a clock
transfers between real sites. In particular, when `synthetic_study_c` is the
external study, development contains only A and B. Their study-only genes are
unavailable in the opposite leave-one-study-out fold and deliberately have
opposite age slopes. This is a useful assay-coverage stress case, but can
produce large opposite fold biases and should not be interpreted as clock
performance. To isolate that mismatch, temporarily set the two
`study_specific_age_effects` slopes to zero. The sample count is allocated
equally over nine age bins, so `--samples-per-study` must be divisible by nine;
for example, 360 creates the tenfold diagnostic fixture.

The compact unit-test fixture in `tests/test_age_clock.py` counterbalances
female and male labels across studies at every age while retaining both sexes
within each study. Its default-fit test checks that the sex-aware clock is
actually fitted when the required metadata are available. The larger synthetic
age-clock smoke test uses the separate planted-expression pseudobulk generated
by `synthetic_de.py`.

Keep exploratory inputs and reports in a dedicated ignored workspace directory
such as `output/experiments/<experiment-name>/`, not only under system-wide
`/tmp`: files outside the repository may be inaccessible to browsers or the
shared workspace UI. Inspect gene-intersection loss, missing-gene coverage,
fold bias, and top coefficients before interpreting feature selection. Any
retained fixture change must preserve the default smoke-test contract, update
`adata.uns["synthetic_test_data"]`, and have focused unit-test coverage.

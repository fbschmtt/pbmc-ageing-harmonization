# Testing the workflows

[Documentation](../README.md)

## Running and verification

Use Make targets rather than raw Nextflow or Docker commands. For routine
changes, choose the smallest focused validation target; reserve `make verify`
for cross-cutting or release-level checks:

The root `Makefile` owns common defaults and includes `make/checks.mk`,
`make/workflows.mk`, and `make/site.mk`. The public target names are preserved;
`scripts/image_context_hash.py` defines and hashes each Docker build context.
Default unit tests require no external study inputs or AIFI binaries. The
separate `make test-models` target runs compatibility checks against locally
acquired models. CI runs Python static/unit/documentation checks and synthetic
frontend checks without acquiring scientific inputs. A separate Nextflow CI job runs
workflow lint, publication, selection, and failure-isolation checks.

`make workflow-test-publication` exercises the production DE process definitions
with small Nextflow stubs. It checks directory staging and verifies that model
tables, diagnostic files, nested WGCNA artifacts and reports remain available
after publication. It uses no scientific inputs or Docker containers, keeps its
artifacts under `output/tests/`, and is included in `make verify`.

`make workflow-test-selection` checks valid selection and pairing, then invokes
the production cell-type workflow with empty, unknown, mixed, duplicate, and
unmatched requests. It asserts failure before any task starts or output is
published. These checks use placeholder files, require no Docker or scientific
inputs, and keep their logs under `output/tests/`.

`make workflow-test-failures` checks that failed fits and renders do not stop
successful siblings or replace their artifacts. It also exercises an entirely
failed cell-type branch, which must still publish a partial manifest. See
[orchestration](orchestration.md) for the failure policy and per-run summaries.

`make run-cell-type-analysis-test` checks the newly published manifest after
Nextflow finishes. It rejects stale manifests, ignored task failures, missing
requested types, changed reports, and missing or corrupt site images. An
`insufficient_cells` report remains a valid outcome; a missing report does not.

The documentation check discovers pages recursively under `docs/`, resolves
links relative to each page, validates heading fragments, and checks every Make
target in documented commands, including targets declared in `make/*.mk`.

```bash
make help
make lint test-unit
make workflow-lint
make run-test STUDIES=wang25
make run-cell-type-test CELL_TYPE=cd14-monocyte
make run-all STUDIES=all
make split-cell-types-test
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
`run-de` fits pseudobulk models, fans out per-type DE report tasks from the fit
results, and renders the combined cross-cell-type analysis report separately.
With `-resume`, completed PyDESeq2 fits remain cached when only a downstream
report input changes. Use `make render-de-reports` to rerun per-type report
notebooks from existing model outputs; this also runs one render task per type
and skips PyDESeq2 fitting. The preceding diagnostics task reuses validated
published calculations or regenerates them when their scientific inputs or
implementation changed; notebook rendering itself only reads their results.
Core, integration-benchmark, DE, and cell-type report targets update
`<outdir>/index.html` after publishing. `make render-report-index OUTDIR=<outdir>`
rebuilds the index using the cell-type manifest when available, or scanning
published per-type reports otherwise, which keeps its links synchronized when
cell types change. Harmonization and benchmark reports are linked when their
published files are present.
`make run-all-test` invokes one composed Nextflow graph. Its synthetic pseudobulk
is generated once for DE, clocks, and gene profiles. `make run-analysis-test`
exercises the same downstream graph using existing core test merges. Both finish
with strict checks on the current run summary and published outputs.

`make verify` is a broad, non-resumed check for cross-cutting changes and
release validation; routine changes should use the focused targets above.
It runs linting, unit tests, Nextflow lint, documentation checks, image builds,
and fresh Docker workflows for the representative test set by default. After
the core merge it checks the integration benchmark, generates synthetic
pseudobulk data for positive DE fits, checks the completed models and planted
markers, runs the default one-method synthetic age-clock smoke test, profiles
the top two merged age-DE genes per cell type, exercises the synthetic
gene-profile report, renders per-type DE reports, and analyzes three
single-cell types by default. Set `VERIFY_CELL_TYPES` to choose those types.
Pass `STUDIES=all` for every available fixture. The synthetic fixture has
only two cell types, so other reports exercise the ordinary no-DE path. It is
workflow coverage rather than a realistic model of study balance or biological
effects: covariate availability and gene presence are deliberately simple, and
the planted effects have no study-specific heterogeneity.
Its synthetic WGCNA output currently has two modules. Thus verification executes
the ordinary WGCNA report path, but does not exercise
the 25-module display cap or provide a production-scale runtime benchmark.
`make run` is the resumable production entry point and merges single-cell data
by default. Set `MERGE_SINGLE_CELL=false` to skip that optional branch.
The current Nextflow profiles leave CPU and time requests unspecified except
for one CPU per parallel cell-type DE task; they do not set per-process memory
limits. Task traces and summaries are published under `<outdir>/provenance/workflows/`.
Use measured peak RAM before sizing another runner or adding resource directives.
Clock tasks run one type at a time while retaining their configured internal
CV parallelism; gene model preparation and profile calculations are capped at two
concurrent tasks each.

## DE test fixtures

The Nextflow `test` profile used by `make run-de-test` enables a test-only mode
that bypasses only the minimum-cell cutoff. It keeps the adult-age, metadata, and model-estimability
checks, so the small fixture can exercise PyDESeq2 without changing the
production inclusion rules. DE test-mode runs are labeled `test_only`; their
DE metadata, result CSVs, and DE reports carry a warning that the results are
for software validation and must not be treated as biological evidence. Use
`make run-de` for the configured analysis rules.

The small core fixture may produce no DE fits because bypassing a cell cutoff
does not add independent samples or fix a rank-deficient design.
`make run-de-synthetic-test` instead generates a deterministic, test-only
pseudobulk fixture with balanced age support, planted covariate and
study-specific age signals, gene-availability gaps, and technical-label
contrasts. It validates fit and report wiring—not biological power or expected
real-data results. The fixture design and marker contract are documented in
[fixture coverage above](#running-and-verification). Do not
upsample real test samples to simulate replication; more cells within a sample
do not add independent samples.

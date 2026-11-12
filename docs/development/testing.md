# Testing the workflows

[Documentation](../README.md) · [Code-to-test map](test-map.md)

## Choosing validation scope

Run the smallest checks that exercise the changed behavior. Start with the
[code-to-test map](test-map.md), which names the relevant test files, workflow
checks, and their limits. Choose tests by the contract they exercise, not just
by filename or the number of files in the diff.

1. Identify what changed: documentation, a Python function, CLI/routing,
   scientific calculation, report execution, publication, or browser behavior.
2. Select the affected tests or test nodes. Include callers of a changed shared
   helper. For Python code, add `make lint`; for docs or Make-target changes,
   add `make docs-check`.
3. Add a workflow or browser check only for behavior the selected unit tests
   cannot exercise. Say what gap the more expensive check will cover before
   starting it. Reuse existing fixtures and prepared artifacts.
4. Stop when those checks pass. Broaden only after a new change, a failure,
   or an unresolved coverage gap. Do not repeat a passed workflow merely to
   collect reassurance, or run every target that could touch the edited code.

Documentation-only changes normally need **only `make docs-check`**. Editing
commands in a guide does not require launching the production analyses they
illustrate. If the documentation checker itself changes, also run its unit tests.

Mode names, branch defaults and Make argument routing normally need focused
orchestration unit tests, `make workflow-test-modes`, and `make docs-check`.
Changing a command's routing does not by itself justify rerunning DE, clocks,
gene fitting, or notebooks. Dataflow, publication and failure-policy changes
need the corresponding stub contract check; add an actual scientific run when
changed runtime behavior lies outside the stubs.

Reserve `make verify` for an explicit request, release validation, or changes
across scientific branches whose combined behavior cannot be covered with
focused checks. A large diff is not sufficient reason on its own. It builds
images and executes a fresh composed scientific graph plus frontend suites.
`make run-all-test` and `make run-analysis-test` also run multiple scientific
branches; they are not lightweight routing checks.

## Focused unit tests

Use Make entry points to retain the locked repository-local environment:

```bash
make test-unit PYTEST_ARGS=tests/test_workflow_orchestration.py
make test-unit PYTEST_ARGS="tests/test_config.py tests/test_cell_type_manifest.py"
make test-unit PYTEST_ARGS=tests/test_config.py::test_staged_analysis_config_needs_no_registries_and_preserves_digest
make test-unit PYTEST_ARGS="tests/test_workflow_orchestration.py -k make_routes"
```

`PYTEST_ARGS` selects test files, nodes, or pytest options. Omit it to run the
full deterministic unit suite. Default unit tests need no external study inputs
or AIFI model binaries; some execute small fits while others mock fitting or
notebook tools. `make test-models` separately checks locally acquired AIFI
binaries; use it when model compatibility is affected.

For frontend logic, select its adjacent Vitest file:

```bash
make site-test-unit SITE_TEST_FILE=src/app/navigation.test.ts
```

For changed interactions, add a matching browser journey using
`make site-test-browser-focus SITE_BROWSER_FILTER=<journey-fragment>`, or a
related navigation/release group. See the [frontend test map](test-map.md#frontend-and-export).
Use `make site-test site-test-browser` when shared release contracts or changes
across views need full coverage.

## Workflow contract checks

These use Nextflow and placeholder inputs or stubs, without Docker or scientific
data. The real preflight check additionally uses the locked repository `.venv`. Select the row matching the changed contract; do not run every row by default.
Their logs/artifacts remain under `output/tests/`.

| Target | What it exercises | What it cannot establish |
|---|---|---|
| `make workflow-lint` | Syntax/structure of production and test entry points. | Runtime inputs, calculations or output publication. |
| `make workflow-test-modes` | Production `main.nf` graph previews: harmonization/analysis/full branches, boolean overrides, benchmark defaults, required merges and preflight rejection. | Task execution, scientific results or notebook rendering. |
| `make workflow-test-preflight` | Real Python validators in the locked `.venv`: unchanged valid config reaches a tiny consumer; invalid analysis/clock settings or core registries are fatal before scientific tasks, including conversion. | Scientific calculations or data-dependent eligibility. Requires the local dev environment. |
| `make workflow-test-selection` | Valid cell-type subsets and production preflight rejection of empty, unknown, duplicate or unmatched inputs. | Cell-type fitting or rendering. |
| `make workflow-test-publication` | Production DE processes with stubs: raw tables, nested diagnostics and report artifacts coexist after publication. | The correctness of calculations inside those tasks. |
| `make workflow-test-failures` | Failed fits/renders and an all-failed cell-type branch: siblings continue, scientific artifacts survive, partial summaries remain visible. | Real scientific/task input compatibility beyond stub outputs. |
| `make workflow-test-main-failures` | Actual `main.nf` analysis-mode connections: a single DE fit, DE report or first benchmark method fails; siblings and unrelated branches finish. Cell reports receive current DE tables independently of DE report success. Strict validation rejects the partial run. | Scientific calculations; fixtures use production process stubs. |
| `make workflow-test-resume` | Production gene-profile graph with stubs: shared preparation, identical-run reuse and render-only invalidation. | Scientific fitting or automatic cache reuse across different entry points. |

CI and `make verify` run the contract checks together. Local iteration can use
just the affected check. See [orchestration](orchestration.md) for failure policy,
publication ownership and run summaries.

Some checks deliberately fail stub tasks or reject invalid configuration.
They print `[EXPECTED FAILURE TEST]` before the case and `[PASS]` only when
its assertions succeed. Raw error output is retained at the printed
`console.log` path. An unexpected outcome prints `[FAIL]` with diagnostics
and causes Make to fail. A recorded partial/failed workflow inside one of these
cases is the expected fixture outcome, not the verdict for `make verify`.

## Focused scientific checks

Use these when a changed calculation, task input contract, or report template
requires actual execution. Check existing artifacts first. Missing inputs
justify creating only the missing prerequisite, not rebuilding the whole graph.
The default single cell-type selection is `cd14-monocyte`.

| Changed behavior | Focused check | Inputs and selection |
|---|---|---|
| Study adapter, harmonization or core task contract | `make run-test STUDIES=<affected-study>` | Existing small fixture and required annotation models. Create a missing fixture with `make test-data STUDIES=<affected-study>`. |
| Local cell-type calculations | `make run-cell-type-test CELL_TYPE=cd14-monocyte` | Reuse its existing `split.h5ad`. If splitting changed, use `make run-cell-type-analysis-test CELL_TYPES=cd14-monocyte` against the existing merge. |
| Cell-type report template | `make render-cell-type-test CELL_TYPE=cd14-monocyte` | Reuse its existing split and analysis; no need to refit local analysis or core. |
| DE report template or prepared diagnostic consumption | `make render-de-reports-test CELL_TYPE=cd14-monocyte` | Existing DE tables and abundance input. Diagnostics reuse valid prepared calculations; missing/invalid diagnostics are regenerated. |
| DE/trajectory fitting or inclusion rules | `make run-de-synthetic-test check-synthetic-de-results check-synthetic-de-reports` | Generates its synthetic input independently of core. If the DE-to-cell-report contract changes, also rerender the ordinary cell-type report to check that those artifacts reach it. |
| Gene-profile calculation/report | `make run-gene-profile-synthetic-test` | One planted gene with output checks. For ranking/fan-out changes, use `make profile-top-age-genes-test` with existing synthetic DE outputs. |
| Clock calculation/report | `make run-age-clock-synthetic-test check-age-clock-synthetic-test AGE_CLOCK_METHODS=elastic_net AGE_CLOCK_ALPHAS=0.1 AGE_CLOCK_L1_RATIOS=0.5 AGE_CLOCK_N_JOBS=1` | One type and small search grid. Select the affected method when changing a different clock. |
| Integration benchmark calculation/report | `make run-integration-benchmark-test check-integration-benchmark-test` | Existing test single-cell merge; executes configured integration methods without other scientific branches. |
| Scientific artifact-to-site contract | `make site-test-workflow` | Existing test workflow outputs. Builds an isolated export; no scientific fitting or replacement of local site data. |

`make run-cell-type-analysis-test` validates its newly published manifest after
Nextflow returns. It rejects stale success, ignored task failures, missing
requested types, and missing/corrupt report images. Other focused report targets
still require inspection of their task exit codes and expected deliverables;
starting Nextflow successfully is not proof of a passing run.

Notebook wrapper unit tests mock the toolchain. For a notebook-template change,
execute the affected render and inspect its HTML/prepared artifacts. For
scientific changes, a syntactically valid notebook or stub output is not enough.
Use synthetic DE to obtain positive fitting coverage with production inclusion
filters; the ordinary small study merge may legitimately yield no DE fits.

## Running and verification

The root `Makefile` owns common defaults and includes `make/checks.mk`,
`make/workflows.mk`, and `make/site.mk`. Docker images use build-context-derived
tags. CI runs Python/static/docs, synthetic frontend, and Nextflow contract
checks without downloading scientific inputs.

`make test-data`, `make run-test`, `make run-all-test`, and `make verify`
default to `aida25,cima,terekhova26`, covering the AIDA supplementary-metadata
join, CIMA counts-layer input, and Terekhova26 raw-count input. Unavailable
fixtures are skipped. Select one affected study for focused checks; use
`STUDIES=all` only when all available adapters/inputs need coverage. Production
`make run` defaults to all available studies and remains resumable. Fresh test
runs and `make verify` must not acquire `-resume`.

Use existing splits, analysis artifacts, and DE tables for report work.
`run-cell-type-analysis-test` must not gain an implicit core-workflow prerequisite.
When working in another output root, pass `OUTDIR=<path>` and reuse its inputs;
keep exploratory inputs, logs and outputs under `output/experiments/<name>/`.

When several changed scientific branches need to be exercised together,
`make run-analysis-test` runs the downstream graph using existing core merges.
It still refits DE, clocks, gene profiles and cell analyses, renders reports,
and runs benchmarking. Use `make run-all-test` when the fresh core-to-analysis
composition also needs validation. Both end with strict artifact checks and
use synthetic pseudobulk generated once for positive DE, clock and gene checks.

`make verify` adds full unit/model/static/documentation checks, workflow
contracts, frontend unit/browser checks, a fresh composed Docker run, and the
workflow-backed atlas check. `VERIFY_CELL_TYPES` selects its three default
single-cell types. The synthetic fixture has only two types; other reports
exercise the ordinary no-DE path. Its synthetic WGCNA currently has two modules,
so this run does not cover the 25-module display cap or production-scale runtime.
The fixture has balanced age support and planted age, sex, BMI, CMV, smoking,
and study effects, but no realistic study-specific heterogeneity.

This synthetic fixture provides positive fitting coverage; a passing run does
not establish that sparsely supported study subsets can fit every model.
Gene-profile unit tests separately cover nullable text with missing values,
rank-deficient and saturated optional designs, and retention of independent
per-study age-bin means. The workflow failure checks above use injected stub
failures to test orchestration, not numerical fitting. Tests of report wrappers
mock notebook execution; changes to unavailable-model report sections still
need a focused real render using the affected prepared artifacts.

Before the existing artifact validators, composed tests also run
`make check-workflow-inventory`'s check. It requires the expected current
processes (including core only for full mode), successful task exit codes, the
requested cell/clock selections, and one preparation per profiled type. It reads
this run's gene-selection files and requires matching fit/render tasks and
published outputs for every selected gene plus the planted verification gene.
It compares those outputs with their current producers, so old published files
cannot substitute for missing tasks. It does not pin the graph's total task count.

Rerun this read-only check against an existing completed composed test without
any scientific or Docker execution:

```bash
make check-workflow-inventory OUTDIR=output/test
```

Use `PIPELINE_MODE=analysis` when the latest summary came from
`run-analysis-test`. Keep its task work directories: the check reads their exit
codes, selections and producer artifacts. The ordinary published outputs remain
usable independently of this test's work-directory checks.

Every Make workflow uses the shared runner to record task trace/exit codes and
`complete`, `partial`, or `failed` summaries under
`<outdir>/provenance/workflows/`. Strict validators reject partial current runs
even when older published artifacts exist. Confirm task `.exitcode` and expected
outputs. The runner refreshes `<outdir>/index.html` once, including after failure;
`make render-report-index` refreshes the inventory without running notebooks.

For diagnosis, `make workflow-failures OUTDIR=output/test` reads the latest
saved test run and lists task errors and failed DE models, including fitting
errors recorded by tasks that exited zero. This is a read-only inspection
command, not a substitute for strict artifact validation. See
[failure inspection](orchestration.md#failure-isolation) for archived runs and
native Nextflow commands.

Frontend tests use isolated synthetic releases. `make site-test-workflow` uses
existing workflow outputs without changing the local export. `make site-build`
validates/builds the current export; `make site-export` replaces it and must not
be used merely for testing.

Report the checks run and their outcomes, plus any relevant coverage gap or
broader validation deliberately left out. Do not imply that a unit, preview, or
stub check validated scientific execution.

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

# Repository guidance

## Environments and dependencies

- Nextflow is expected to be available on `PATH`. Agent-only local fallback:
  it is also available in the `scanpy` Conda environment. Do not encode that
  fallback in repository scripts, Make targets, or user-facing documentation.
- Keep dependency installs local to this repository. Use the existing `.venv`
  for Python development dependencies (for example,
  `.venv/bin/python -m pip install -e '.[dev,qc]'`); do not install Python or
  system packages globally.
- Docker image tags are derived from a hash of the current build context. Build
  or run them through the Make entry points when Docker execution is requested.
- `requirements.lock` is the exact runtime/report environment for the Python
  image. `requirements.dev.lock` extends it with local test and lint tools; use
  the Make install target rather than resolving an ad-hoc environment.

## Workflow checks

- Start with the [Choosing a Make target](README.md#choosing-a-make-target)
  table in `README.md` to select the smallest validation target; the rules
  below add mandatory constraints and take precedence if they conflict.
- Use Make targets as the entry points for linting, testing, Docker builds, and
  workflows. Choose the smallest target that covers the change; do not run
  `make verify` by default. `make verify` is the broad, non-resumed Docker
  verification target and uses Docker layer caching plus all configured test
  studies unless `STUDIES` is explicitly supplied.
- Use `make run` only for resumable production runs. Do not add `-resume`
  to `make verify` or `make run-test`.
- Run `make lint test-unit` for Python-only changes and `make workflow-lint`
  for Nextflow/configuration-only changes. Use a focused workflow target (for
  example, `make run-test STUDIES=wang25` or a targeted cell-type report) when
  wiring must be exercised. Reserve `make verify` for broad cross-cutting,
  release-level, or deliberately full-suite validation; report the focused
  checks run and any relevant broader check not run.
- For a change isolated to one cell-type report, reuse an existing split with
  `make run-cell-type-test CELL_TYPE=<slug>` rather than rebuilding the core
  workflow or all reports. If no type is specified, use `cd14-monocyte`: it is
  well represented in the bundled multi-study fixture. If the split is absent
  but the test merge exists, run `make split-cell-types-test`; if the merge is
  absent, use `make run-cell-type-analysis-test`. Reserve the full all-type
  test workflow for changes that affect splitting, shared analysis, or publication.
- Use `make run-de-synthetic-test` for positive DE model coverage with production
  inclusion filters. Then use `make run-cell-type-analysis-test-existing` to
  check that synthetic DE artifacts reach the ordinary cell-type reports.
- `make verify` runs the fresh core test once, generates and fits the synthetic
  DE fixture, then renders all cell-type reports from that merge with matching
  synthetic DE results. Do not reintroduce a core-workflow prerequisite on
  `run-cell-type-analysis-test-existing`.
- When adding a required CLI argument or published output, find every invocation
  with `rg` (Make, Nextflow, tests, and docs) and run the smallest public Make
  target that exercises the changed contract.
- Treat Nextflow submission as incomplete validation. Confirm task `.exitcode`
  and the published deliverables before reporting a focused workflow as passed.
- Run `make docs-check` after editing tracked Markdown documentation or Make
  target names.
- Development test targets omit `nehar_belaid26` by default because its larger
  fixture dominates runtime. This is a Make-level selection only; use
  `STUDIES=all` for all-study coverage. `make verify` is the exception: it
  restores all configured studies by default, while an explicit `STUDIES` list
  still permits a narrower verification run.

## Working style

- Preserve existing user changes in a dirty worktree.
- Ask a clarifying question before making an assumption that would materially
  change scientific analysis, workflow behavior, or output semantics.

## Study adapter evidence

- Do not comment every canonical field mapping. Source-backed mappings are
  self-evident from their source column or supplement.
- Clearly mark inferred values or transformations with a short adjacent
  `# inferred:` comment, especially assumptions such as genome version or
  thawed/frozen status.
- Keep the corresponding assumption in the study configuration provenance as
  well; do not silently present an inferred value as reported metadata.

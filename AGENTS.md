# Repository guidance

See the [documentation index](docs/README.md); these rules take precedence over guide examples.

## Working style

- Preserve existing user changes in a dirty worktree.
- Ask before making an assumption that would materially change scientific
  analysis, workflow behavior, or output semantics.
- Keep exploratory inputs, diagnostics, and reports under
  `output/experiments/<name>/` so they remain accessible in the shared workspace.

## Environments and dependencies

- Use `make install` for Python development dependencies in the repository-local
  `.venv`; do not install Python or system packages globally.
- `requirements.lock` defines the exact runtime/report environment;
  `requirements.dev.lock` extends it for development. Do not resolve an ad-hoc environment.
- Build and run Docker through Make; preserve its build-context-derived image tags.

## Architectural boundaries

- Put scientific calculations in analysis modules and explicit workflow stages.
  Report templates consume prepared artifacts through the shared notebook runner.
- Export site data from named, versioned artifacts and manifests; do not infer
  content from notebook headings, cell positions, HTML, or source-code comments.
- Keep independent analysis/report failures isolated: allow siblings to finish,
  retain successful outputs, and record partial runs. Keep essential input/configuration
  and core-processing failures fatal. Verification may reject partial results after completion.
- Preserve output contracts and provenance during refactors. When changing a CLI
  argument or published artifact, find all invocations with `rg`, including Make,
  Nextflow, tests, and documentation, and exercise the changed contract.

## Validation

- Run the smallest tests that cover the changed behavior. Use the
  [code-to-test map](docs/development/test-map.md) to select files and contracts,
  and the [testing guide](docs/development/testing.md#choosing-validation-scope)
  to decide whether integration coverage is needed. Python changes normally use
  `make lint` and `make test-unit PYTEST_ARGS="tests/test_<area>.py"`; narrow to a
  test node or add related files when appropriate. Documentation-only edits need
  `make docs-check`, not scientific or frontend runs.
- Entry-point modes and Make routing use focused orchestration unit tests,
  `make workflow-test-modes`, and `make docs-check`. Do not rerun calculations
  or reports unless the change also affects task inputs, dependencies, or outputs.
- Broaden validation only for an identified coverage gap, a failure, or a new
  change. Stop once the relevant checks pass. Before starting a costly workflow,
  state which changed contract it will exercise and why narrower checks cannot.
  Reserve `make verify` for an explicit request, release validation, or changes
  spanning multiple scientific branches that focused checks cannot cover.
- Keep production `make run` resumable. Do not add `-resume` to `make verify` or `make run-test`.
- Reuse generated fixtures, splits, and analysis artifacts for focused report checks;
  default to `cd14-monocyte`. Do not add a core-workflow prerequisite to
  `run-cell-type-analysis-test`.
- Use synthetic DE checks for positive fitting coverage with production inclusion
  filters, and check that their artifacts reach the ordinary cell-type reports.
- For producer/consumer wiring, publication ownership, or failure-policy changes,
  run `make workflow-test-failures`. Add only the affected scientific workflow
  when changed runtime inputs/calculations need coverage beyond its stubs;
  entry-point routing alone does not require scientific execution. For changes to
  main-entry failure propagation or DE handoffs, also use
  `make workflow-test-main-failures` (stubs, no scientific fitting).
- Keep configuration preflight fatal and ahead of scientific fan-out, including
  conversion. For changes to validation or its gates, use focused preflight/config
  unit tests and `make workflow-test-preflight`; use the main stub failure check
  when the dependency barrier changes.
- For DE publication or staging changes, run `make workflow-test-publication`
  and the focused scientific workflow when calculations or inputs also change.
- Nextflow submission alone is not validation: confirm task `.exitcode` and
  published deliverables. Use the artifact validators to detect ignored failures,
  stale manifests, missing selections, and missing outputs.
- For frontend changes, choose affected Vitest files with `make site-test-unit
  SITE_TEST_FILE=<path>` and, for browser behavior, a matching
  `make site-test-browser-focus SITE_BROWSER_FILTER=<journey>` or browser group.
  Use `make site-test site-test-browser` for shared release/data contracts or
  cross-view changes needing full coverage. Use `make site-test-workflow` for
  existing workflow outputs; do not replace the local site export just to test.
  `make site-build` validates it; `make site-export` replaces it.
- Run `make docs-check` after documentation or Make-target changes. Report checks
  run, their outcomes, and any relevant broader validation not performed.

## Study adapter evidence

- Keep source-backed field mappings concise; do not comment every mapping.
- Mark inferred values or transformations with an adjacent `# inferred:` comment,
  especially genome version and frozen/thawed status, and record the same
  assumption in study-configuration provenance. Do not present inferred metadata
  as reported facts.

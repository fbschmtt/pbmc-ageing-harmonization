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

- Use Make entry points and the [target table](docs/development/README.md#choosing-a-make-target)
  to choose the smallest relevant check. Follow the [testing guide](docs/development/testing.md)
  for fixture setup and broader coverage; reserve `make verify` for broad changes,
  release checks, or an explicit request.
- Keep production `make run` resumable. Do not add `-resume` to `make verify` or `make run-test`.
- Reuse generated fixtures, splits, and analysis artifacts for focused report checks;
  default to `cd14-monocyte`. Do not add a core-workflow prerequisite to
  `run-cell-type-analysis-test`.
- Use synthetic DE checks for positive fitting coverage with production inclusion
  filters, and check that their artifacts reach the ordinary cell-type reports.
- For orchestration, publication ownership, or failure-policy changes, run
  `make workflow-test-failures` and the relevant focused scientific workflow.
- For DE publication or staging changes, run `make workflow-test-publication`
  and the focused scientific workflow when calculations or inputs also change.
- Nextflow submission alone is not validation: confirm task `.exitcode` and
  published deliverables. Use the artifact validators to detect ignored failures,
  stale manifests, missing selections, and missing outputs.
- For frontend changes, run `make site-test site-test-browser`. Use
  `make site-test-workflow` for existing workflow outputs; do not replace the local
  site export just to test. `make site-build` validates it; `make site-export` replaces it.
- Run `make docs-check` after documentation or Make-target changes. Report checks
  run, their outcomes, and any relevant broader validation not performed.

## Study adapter evidence

- Keep source-backed field mappings concise; do not comment every mapping.
- Mark inferred values or transformations with an adjacent `# inferred:` comment,
  especially genome version and frozen/thawed status, and record the same
  assumption in study-configuration provenance. Do not present inferred metadata
  as reported facts.

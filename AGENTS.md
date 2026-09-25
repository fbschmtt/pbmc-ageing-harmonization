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

## Workflow checks

- Use Make targets as the entry points for linting, testing, Docker builds, and
  workflows. `make verify` is the required full verification target for
  workflow, container, or notebook changes; it uses Docker layer caching and
  runs the all-study Docker test workflow without `-resume`.
- Use `make run` only for resumable production runs. Do not add `-resume`
  to `make verify` or `make run-test`.
- After a large feature or architectural change, run the smallest end-to-end
  Make target that exercises the affected container, Make, and Nextflow wiring.
  Consider the full verification suite (`make verify`) for broad cross-cutting
  or release-level changes, balancing its runtime and resource cost; report any
  relevant verification that was not run.

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

# Development and Make targets

[Documentation](../README.md)

## Local Python development

Create the repository-local environment and run the data-independent checks:

```bash
make install
make lint test-unit docs-check
```

Dependencies are pinned in `requirements.lock` and `requirements.dev.lock`.
Python 3.12 matches the locked development environment; `.venv` stays local to
the repository.
`make test-models` separately checks locally acquired AIFI model binaries.
See [workflow testing](testing.md) for fixture-based checks and full verification.

## Nextflow and Docker

Nextflow runs the configured studies in parallel and caches completed tasks.
The Docker profile is the supported workflow environment and is selected by
default. The test profile uses isolated 200-cell fixtures from `test_data/`;
RDS fixtures still use the normal conversion process. `make test-data` creates
only missing fixtures and never changes `input_data/`; use
`TEST_DATA_OVERWRITE=true` to refresh them. Tests default to
`aida25,cima,terekhova26`; pass `STUDIES=<list>` (including `STUDIES=all`) to
change that selection. Production `make run` still defaults to all studies.

## Choosing a Make target

Start with the [code-to-test map](test-map.md) and the
[validation scope guide](testing.md#choosing-validation-scope). Select checks
for the changed contract; do not use a whole scientific run to validate a mode
name, Make argument, or documentation edit. Broaden only for a stated coverage
gap, failure, or additional change, and stop once the relevant checks pass.

| Changed behavior | Start with | Add broader coverage only when |
|---|---|---|
| Documentation, including AGENTS.md | `make docs-check` | The checker itself changed: add `make test-unit PYTEST_ARGS=tests/test_docs.py`. Example production commands do not need execution just because their explanation changed. |
| Python helper, configuration, or CLI | `make lint` plus `make test-unit PYTEST_ARGS="<affected test files or nodes>"` from the [test map](test-map.md) | Callers or runtime contracts also changed; add those tests or one affected scientific task. Full unit tests are useful when a shared change affects many consumers. |
| Entry-point modes or Make routing | `make test-unit PYTEST_ARGS=tests/test_workflow_orchestration.py workflow-test-modes docs-check` | Task inputs, dependencies or output contracts changed beyond routing; use their checks below. |
| Nextflow definitions or container profiles | `make workflow-lint`, then the affected contract check from the [workflow table](testing.md#workflow-contract-checks) | Stubs/previews cannot cover the changed runtime command, environment or data contract; execute that branch only. |
| Cell-type selection or split/analysis pairing | `make workflow-test-selection` | Runtime splitting/analysis also changed: add one selected cell type against the existing merge. |
| DE staging or publication ownership | `make workflow-test-publication`; add `make workflow-test-failures` for ownership/failure changes | Calculations or actual task inputs also changed; add focused synthetic DE. |
| Producer/consumer wiring or failure policy | `make workflow-test-failures` and affected unit tests | A changed runtime input/calculation is not covered by the stubs; add the corresponding focused scientific check. |
| Failure propagation through the main entry point | `make workflow-test-main-failures` | Actual calculations or CLI inputs changed beyond the stubs; add just that affected scientific task. |
| Completed composed test task/selection coverage | `make check-workflow-inventory` using the existing output root | Scientific results changed; reuse the corresponding artifact validators before deciding whether a new calculation is needed. |
| Gene-profile caching/template invalidation | `make workflow-test-resume` | Scientific model/cache contents also changed; add profile unit tests and one synthetic gene run. |
| Scientific calculation or report template | Affected unit tests and one [focused scientific check](testing.md#focused-scientific-checks) | A shared producer-to-consumer contract affects more branches; test those consumers. Reuse existing inputs and prepared artifacts. |
| Frontend logic or interaction | Matching `make site-test-unit SITE_TEST_FILE=<path>`; add the affected browser journey/group for interactions | Shared data/release contracts or changes across views need `make site-test site-test-browser`. |
| Existing scientific outputs consumed by the atlas | `make site-test-workflow` | Required workflow outputs are absent or changed; create only those outputs. Never replace the local export merely to test. |
| Multiple interacting scientific branches | `make run-analysis-test` when focused checks cannot cover their interaction | Fresh core-to-downstream composition is also affected: use `make run-all-test`. Both execute substantial calculations and rendering. |
| Release or explicitly requested full verification | `make verify` | Full verification is deliberate, not the default after an edit. Use `STUDIES=all` only when every available fixture needs coverage. |

Production commands are documented in [workflow modes](orchestration.md#entry-points).
They are not substitutes for selecting focused validation.

Nextflow submission alone is not a passing workflow check: confirm each task's
`.exitcode` and the published deliverables. Development test targets use the
representative AIDA25/CIMA/Terekhova26 fixtures and skip unavailable inputs;
pass an explicit comma-separated `STUDIES` list to change a run.

Harmony integration uses the pinned Harmony2 C++ backend. It uses the backend's
native thread policy; the legacy OpenMP/OpenBLAS thread overrides that crashed
on large production runs have been removed.

Single cell analysis reports provide local embeddings, marker contrasts, study
coverage, and local-cluster composition views. Separate per-type DE reports and a
combined cross-cell-type analysis report are generated by `run-de`; those DE
reports need only the merged pseudobulk input.

`make run-test` and `make run` create the single-cell merge by default. To skip
the production merge when it is not needed, set `MERGE_SINGLE_CELL=false`:

```bash
make run STUDIES=all MERGE_SINGLE_CELL=false
```

To benchmark global embedding and type-calling methods without rerunning that
merge:

```bash
make run-integration-benchmark
```

Study-QC, merge-QC, integration-benchmark, cell-type, DE, and gene-profile
report templates are Jupytext Python sources under `reports/`. Their runners
materialize a temporary notebook immediately before execution; generated
templates are not tracked. The published `executed.ipynb` is the audit artifact.

The merge QC report summarizes sample-level technology, intronic-read
inclusion, and frozen-status distributions with provenance status. OneK1K's
sample counts use `donor_id`; its `pool_number` is a multiplexed technical
processing group, not the biological sample.

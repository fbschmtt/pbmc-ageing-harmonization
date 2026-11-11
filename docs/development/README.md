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

Choose the smallest target that directly exercises the changed contract. The
repository instructions in [AGENTS.md](../../AGENTS.md) use this table when deciding
which validation to run.

| Change or goal | Preferred target | When to use a broader target |
|---|---|---|
| Python library, CLI, or deterministic test change | `make lint test-unit` | Run a focused workflow too if a task contract or published artifact changed. |
| Nextflow workflow or configuration change | `make workflow-lint` | Add `make run-test STUDIES=<study>` when task wiring must execute. |
| Shared scientific orchestration | `make workflow-lint workflow-test-failures`, then `make run-analysis-test` | Uses existing test merges; `make run-all-test` also exercises fresh core composition. |
| Cell-type selection or split/analysis pairing | `make workflow-test-selection` | Checks valid subsets and production preflight errors without scientific inputs. Add `make run-cell-type-analysis-test CELL_TYPES=<slug,slug>` when analysis wiring changes. |
| DE artifact publication or task staging | `make workflow-test-publication` | Uses production process definitions with tiny test stubs, then asserts that model files, nested diagnostics and rendered reports coexist. Run synthetic DE too when calculations or report inputs change. |
| One single cell analysis report | `make run-cell-type-test CELL_TYPE=<slug>` (requires an existing split; create test splits with `make split-cell-types-test`) | Use `make render-cell-type-test CELL_TYPE=<slug>` to rerender one report from existing analysis artifacts, or `make render-cell-type-reports-test` to rerender all test reports. |
| Rerender single cell analysis reports from existing analysis artifacts | `make render-cell-type-reports[-test]` | Use `make render-cell-type[-test] CELL_TYPE=<slug>` to refresh one report. |
| Rerender per-type pseudobulk DE reports from existing model results | `make render-de-reports[-test] [CELL_TYPE=<slug>]` | Renders all types when `CELL_TYPE` is omitted, or one type when set. Requires current per-type model outputs and `abundance_input.tsv`. A separate diagnostics stage prepares abundance, clustering, WGCNA and enrichment artifacts; unchanged calculations are reused. Rendering reads those artifacts and does not fit models. Rerun `make run-de[-test]` to regenerate missing abundance inputs. |
| Rebuild the central report index | `make render-report-index` | It is also refreshed automatically after core, integration-benchmark, DE, cell-type, and age-clock workflow targets. |
| Cell-type splitting or shared analysis | `make run-test`, then `make run-cell-type-analysis-test CELL_TYPES=<slug,slug>` | Omit `CELL_TYPES` to analyze every type. Use `make run-all-test` only when the composed integration workflow is needed. |
| Positive pseudobulk DE fitting | `make run-de-synthetic-test` and `make check-synthetic-de-results` | This fixture is independent of the core merge and renders per-type DE plus cross-cell-type reports. To also render single-cell reports, run `make run-test` followed by `make run-cell-type-analysis-test CELL_TYPES=<slug,slug>`; omit `CELL_TYPES` to analyze all types. |
| Cell-type transcriptomic age clocks | `make run-age-clock` | Fit all types with `make run-age-clocks`; see [age clocks](../analysis/age-clocks.md#age-clocks) for selection and validation. |
| Focused gene-by-cell-type profile | `make profile-gene CELL_TYPE=<label-or-slug> GENE=<gene-id>` | Use `make run-gene-profile-synthetic-test` to fit and render a planted synthetic gene and check the report end to end. |
| Top merged age-DE gene profiles | `make profile-top-age-genes-test` | Uses the synthetic DE fixture and profiles two genes per cell type; included in `make run-all-test`. |
| DE from an existing test pseudobulk merge | `make run-de-test` | Requires `make run-test` first; the small core fixture may not have enough independent samples for successful fits. Use the synthetic target above for positive fitting coverage. |
| Focused profiles of top merged age-DE genes | `make profile-top-age-genes` | Runs ten profiles per cell type by adjusted p-value; set `TOP_AGE_GENE_COUNT` to change the number. Requires completed merged DE results and pseudobulk input. |
| Core workflow for one or a few studies | `make run-test STUDIES=<list>` | Omit `STUDIES` for the representative AIDA25/CIMA/Terekhova26 set; use `STUDIES=all` for every available fixture. |
| Global integration/annotation diagnostics | `make run-integration-benchmark[-test]` | Run after an existing single-cell merge; `make run-all-test` also exercises its test form. |
| Composed Docker integration test | `make run-all-test` | Includes the synthetic DE, clock, and profile checks plus the default three single-cell types. |
| Static atlas export, frontend, or data contract | `make site-test site-test-browser` | These run in the tracked frontend Docker images and need no host Node installation. Browser tests build an isolated synthetic release. Use `make site-test-workflow` to test existing workflow outputs; `make site-build` validates the existing export. |
| Markdown documentation or Make target names | `make docs-check` | Also run the target described by changed commands. |
| Broad cross-cutting or release validation | `make verify` | Reserve this non-resumed, full workflow for cross-cutting or release checks. Prefer the focused targets above for routine changes. Pass `STUDIES=all` to use every available fixture. |
| Production analysis | `make run` | It is resumable; use only after the relevant focused checks pass. |

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

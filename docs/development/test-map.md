# Code-to-test map

[Documentation](../README.md) · [Choosing validation scope](testing.md#choosing-validation-scope)

Select the row covering the changed behavior, then run its test files through
`make test-unit PYTEST_ARGS="<files or test nodes>"`. Paths in the Python tables
are relative to `src/pbmc_pipeline/`, except paths starting with `scripts/` or
`config/`. Include callers when a shared helper's behavior or contract changes.
These are coverage pointers, not a requirement to run every test in a row.

## Python configuration and core processing

| Changed code | Start with | What these tests check |
|---|---|---|
| `config.py`, `config/*.json` and schemas | [test_config.py](../../tests/test_config.py) | Registry/schema rules, staged analysis configuration without registries, shared inclusion/processing validation, and rejection before data reads. AIFI binary compatibility is separately marked `external_models`. |
| `preparation.py`, `studies/`, study mappings | [test_preparation.py](../../tests/test_preparation.py), affected cases in [test_metadata.py](../../tests/test_metadata.py) | Canonical per-cell rows, adapter joins/labels, cohort exclusions, and exact cell coverage or declared subsets. |
| `metadata.py`, feature repair in `harmonize.py`, count validation | [test_metadata.py](../../tests/test_metadata.py), [test_validation.py](../../tests/test_validation.py) | Missing metadata, feature-label repair/aggregation, and integer-count validation. |
| `harmonize.py`, `embedding.py`, feature/count output | [test_harmonize.py](../../tests/test_harmonize.py), relevant [test_merge.py](../../tests/test_merge.py) cases | Raw-count output, feature-index repair, Harmony adapter orientation, sparse matrix compatibility, and merge embedding inputs. This is not a production-scale integration benchmark. |
| `merge.py`, pseudobulk grouping and gene joins | [test_merge.py](../../tests/test_merge.py) | Sample/type pseudobulks, outer gene joins, gene order, metadata serialization, and preservation of original annotations. |
| `manifest.py` | [test_manifest.py](../../tests/test_manifest.py) | Published artifact/input provenance and manifest contents. |

For a changed adapter or core CLI/task contract, add
`make run-test STUDIES=<affected-study>` using its existing fixture. Run
`make test-data STUDIES=<affected-study>` only if the fixture is missing. These
checks do not establish that inferred study metadata are source-backed; review
the mapping evidence separately.

## Scientific analysis and reporting

| Changed code | Start with | What these tests check |
|---|---|---|
| `cell_type.py` | [test_cell_type.py](../../tests/test_cell_type.py) | Raw-count splits, retained type partitions, and local analysis behavior on small inputs. |
| `cell_type_figures.py`, `cell_type_report.py` | [test_cell_type_report.py](../../tests/test_cell_type_report.py) | Report command/environment wiring and estimable cluster-fraction trends. |
| `cell_type_assets.py`, `cell_type_manifest.py`, shared colors | [test_cell_type_assets.py](../../tests/test_cell_type_assets.py), [test_cell_type_manifest.py](../../tests/test_cell_type_manifest.py), [test_cell_type_colors.py](../../tests/test_cell_type_colors.py) | Named images/checksums, selections and partial manifests, and shared label colors. |
| `differential_expression.py`, `age_trajectory.py`, `synthetic_de.py` | [test_differential_expression.py](../../tests/test_differential_expression.py), [test_config.py](../../tests/test_config.py) when settings change | Inclusion/design rules, age-bin support, covariate/result separation, synthetic fixture properties, and per-type result collection. |
| `de_diagnostics.py`, `abundance.py`, prepared DE report inputs | [test_de_diagnostics.py](../../tests/test_de_diagnostics.py) | Fraction-age effects, prepared tables, source/config/resource invalidation, corruption detection and repair, and rendering from prepared results. |
| `fraction_age.py`, `reference_expression.py` | [test_fraction_age.py](../../tests/test_fraction_age.py), [test_reference_expression.py](../../tests/test_reference_expression.py) | Zero fractions, reference-bin alignment, support/rank-deficiency handling, and study-balanced reference expression with gene availability. |
| `trajectory_analysis.py`, `expression_atlas.py` | [test_trajectory_analysis.py](../../tests/test_trajectory_analysis.py) | Supported/shared bins, clustering inputs, recurrence, residual cluster selection, and expression-atlas gene/depth filters. |
| `residual_wgcna.py`, `go_enrichment.py` | [test_residual_wgcna.py](../../tests/test_residual_wgcna.py), [test_go_enrichment.py](../../tests/test_go_enrichment.py) | Residual-network/enrichment transformations and output contracts on small or mocked inputs. |
| `age_clock.py`, `age_clock_features.py`, `age_clock_cross_type.py` | Affected cases in [test_age_clock.py](../../tests/test_age_clock.py) | Training/external-study separation, normalization/features, inclusion, model selection, prediction coverage, aggregation, and failure handling. |
| `scripts/clock_task.py` | Clock cases in [test_workflow_orchestration.py](../../tests/test_workflow_orchestration.py), relevant [test_age_clock.py](../../tests/test_age_clock.py) cases | Task status, eligibility and cohort provenance, plus the underlying clock behavior. |
| `scripts/profile_gene.py`, `gene_profile_report.py` | [test_gene_profile.py](../../tests/test_gene_profile.py) | Gene availability, prepared-model mismatch rejection, trajectory/reference support, coefficient labels, and report runner failures. |
| `integration_benchmark.py` | [test_integration_benchmark.py](../../tests/test_integration_benchmark.py), relevant core embedding tests | Thin benchmark artifacts retain labels/UMAP while dropping large matrices. This test does not execute every integration method. |
| `notebook.py`, `report_toc.py`, report wrappers | [test_notebook.py](../../tests/test_notebook.py), the affected wrapper's [test_qc.py](../../tests/test_qc.py), [test_merge_qc.py](../../tests/test_merge_qc.py), [test_cell_type_report.py](../../tests/test_cell_type_report.py), [test_age_clock_report.py](../../tests/test_age_clock_report.py), or [test_gene_profile.py](../../tests/test_gene_profile.py) | Tool command order/environment, isolated caches, titles/TOCs, cleanup and failure behavior. Notebook execution is mocked; template changes need an actual focused render. |
| `reporting.py`, `covariates.py`, `report_theme.py` | [test_reporting.py](../../tests/test_reporting.py), [test_report_theme.py](../../tests/test_report_theme.py) | Fraction denominators/zeros, label/reference handling, reporting transformations and presentation helpers. |

Use the [focused scientific checks](testing.md#focused-scientific-checks) when
calculations, runtime inputs, or notebook templates change. Small or mocked unit
tests do not prove that a scientific fit or notebook execution succeeds.
Some wrappers, including `integration_benchmark_qc.py` and `trajectory_report.py`,
have no dedicated runner test file. For changes there, use the shared notebook
tests plus the affected actual report execution; add a focused regression for
new wrapper behavior rather than claiming existing tests cover it.

## Orchestration, publication, and tools

| Changed code or contract | Smallest relevant check | Coverage and limits |
|---|---|---|
| `make/workflows.mk`, `main.nf` mode/routing arguments | [test_workflow_orchestration.py](../../tests/test_workflow_orchestration.py), `make workflow-test-modes` | Make dry runs check paths, one invocation, mode selection, and fresh/resumed flags. Nextflow previews check branches/defaults and input preflight. Neither executes science or notebook reports. |
| `preflight.py`, `modules/configuration.nf`, configuration gates | [test_preflight.py](../../tests/test_preflight.py), `make workflow-test-preflight`; main stub failure checks for the dependency barrier | Consumed settings, clock options and full registry validation; preserves configuration bytes. Invalid config blocks fan-out and conversion before data access. |
| `scripts/run_workflow.py`, `scripts/workflow_diagnostics.py`, `scripts/show_workflow_failures.py` | [test_workflow_orchestration.py](../../tests/test_workflow_orchestration.py), [test_workflow_diagnostics.py](../../tests/test_workflow_diagnostics.py) | Task/application-level partial status, exit codes, one index finalization, and saved/read-only task and model error inspection. |
| `scripts/workflow_test_output.py`, failure-test console capture | [test_workflow_test_output.py](../../tests/test_workflow_test_output.py), affected workflow failure/preflight check | Expected-failure banners, pass-after-assertions behavior, retained raw output and unchanged exit status, plus visible diagnostics and propagation for unexpected errors. |
| `scripts/render_report_index.py` | [test_report_index.py](../../tests/test_report_index.py) | Index discovery, manifests, grouping and report links without rerunning reports. |
| `modules/cell_type_selection.nf`, split/analysis pairing | `make workflow-test-selection`; [test_cell_type_manifest.py](../../tests/test_cell_type_manifest.py) if manifest behavior changes | Valid subsets and rejection of empty/unknown/duplicate selections or unmatched inputs before dispatch. |
| DE staging or `publishDir` ownership | `make workflow-test-publication`; [test_workflow_checks.py](../../tests/test_workflow_checks.py) for validator changes | Production DE process definitions with stubs check coexistence of raw results, nested diagnostics, and rendered artifacts. |
| Task dependencies, ignored failures, aggregation or publication ownership | `make workflow-test-failures` | Stubbed production analysis processes check failed fits, failed renders, all-failed cell analyses, retained siblings and partial summaries. Add the affected scientific branch only for runtime behavior the stubs cannot check. |
| Real `main.nf` failure propagation and DE-to-report handoffs | `make workflow-test-main-failures` | Single-task failures in DE fitting, DE rendering and benchmarking; successful same-branch siblings, independent gene/clock/cell reports, current raw DE consumption and partial summaries. Benchmark siblings queued after the first failure must complete. |
| Composed test task/selection inventory | [test_workflow_inventory.py](../../tests/test_workflow_inventory.py), `make check-workflow-inventory` against existing outputs | Rejects missing branches, unsuccessful tasks, dropped/duplicated gene requests or model preparations, slug collisions and stale publication. Reads current task work directories without fitting or rendering. |
| Gene-profile caching or template invalidation | `make workflow-test-resume` | Real profile graph with stubs checks one preparation per type, identical-run caching, and report-only invalidation. It does not establish cross-entry-point cache reuse. |
| `scripts/check_cell_type_reports.py`, `scripts/check_de_publication.py` | [test_workflow_checks.py](../../tests/test_workflow_checks.py) | Rejection of stale success, ignored failures, missing/corrupt site assets, and missing published stages. |
| `scripts/image_context_hash.py` | [test_image_context_hash.py](../../tests/test_image_context_hash.py) | Build-context-derived image hashes and excluded generated inputs. Build only the affected image when the actual container environment changes. |
| `scripts/check_docs.py` | [test_docs.py](../../tests/test_docs.py), `make docs-check` | Recursive Markdown discovery, relative links/anchors, and Make target discovery including shared recipes. |
| `scripts/download_site_pathways.py` | [test_download_site_pathways.py](../../tests/test_download_site_pathways.py) | Pathway input parsing and generated metadata without scientific workflows. |

For other scripts, look for imports and CLI invocations with `rg` rather than
assuming a similarly named test covers them. For example,
`rg -n 'check_synthetic_de|profile_top_age_genes' tests scripts make workflows`
finds both tests and the relevant artifact validators. Add a focused regression
when existing coverage misses the changed contract.

## Frontend and export

Paths below are relative to `site/`. Use
`make site-test-unit SITE_TEST_FILE=<test-path>` to select one Vitest file.
For changed browser behavior, use
`make site-test-browser-focus SITE_BROWSER_FILTER=<journey-name-fragment>`;
journey names are in [journeys.mjs](../../site/tests/browser/journeys.mjs).
A focused journey still builds an isolated synthetic release and needs Docker.

| Changed code | Start with | What these tests check |
|---|---|---|
| `src/app/` navigation, session and labels | Adjacent `navigation.test.ts`, `session.test.ts`, `releaseLabel.test.ts` | Navigation/selection state, session behavior and release labels; add the affected navigation journey for interactions. |
| `src/data/` loaders and contracts | Adjacent `atlas.test.ts`, `contracts.test.ts`, `cellTypeHierarchy.test.ts`, `differentialDownload.test.ts`, `geneSet.test.ts`, `pointSelection.test.ts` | Asset parsing, selection, downloads, shared gene sets and data contracts. |
| Cross-type summaries and view | `scripts/export-cross-type.test.mjs`, `src/features/crossType/crossTypeView.test.ts`, cross-type browser journeys | Named artifacts, signed effects, tested/significant denominators, optional legacy data, selection/navigation and shared patterns. Workflow browser checks the exported real test results too. |
| Differential, composition and explorer features | Adjacent `*.test.ts` under `src/features/` | Model/plot transformations, filters, composition summaries and scale rules. A browser journey checks actual rendering/interactions. |
| `src/theme/` | Adjacent `format.test.ts`, `axisScale.test.ts`, `studyColors.test.ts` | Numeric formatting, plot scales and stable study colors. |
| `scripts/export-single-cell.mjs`, named image inputs | `scripts/export-single-cell.test.mjs`; Python asset tests if producers change | Versioned manifest/image consumption, metadata and checksum failures. |
| Other exporter, PNG, CSV, clock or provenance helpers | Matching `scripts/*.test.mjs` | The corresponding asset transformation, download values, compression, or provenance contract. |
| Browser/network failure helpers | `tests/browser/network.test.mjs`, affected release journey | Network failure handling plus visible loading/fallback behavior. |
| Shared release schemas, loaders or changes spanning views | `make site-test site-test-browser` | Export/release contracts, all Vitest tests, and browser navigation/release journeys. Uses synthetic inputs, not workflow calculations. |
| Scientific artifact-to-atlas integration | `make site-test-workflow` | Exports and browses existing test workflow results in isolation. Does not run science or replace the local site export. |

Browser groups `make site-test-browser-navigation` and
`make site-test-browser-release` are useful when several related journeys are
affected. `make site-build` validates/builds the current export; it is not a
substitute for interaction tests. Never run `make site-export` just to test:
it replaces the local export.

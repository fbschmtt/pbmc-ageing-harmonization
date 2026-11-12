# PBMC Ageing Harmonization

Combine PBMC single-cell RNA-seq studies into a common dataset and explore how
cell composition and gene expression vary with age.

The pipeline produces harmonized raw-count H5AD files, sample-level pseudobulk
matrices, and HTML reports for quality control and cell-type analyses. An
optional browser-based atlas presents the exported results.

## Before you start

You need **Docker, Nextflow 25.04 or newer, Make, and Python 3** on the host.
The analysis Python and R environments run in Docker; a local Node installation
is optional.

Expression data and trained annotation models are **not bundled**. Choose the
studies you want from the [input catalog](docs/inputs/files.md), then follow
[the acquisition guide](docs/inputs/README.md). Public files can be downloaded
through Make; some study inputs and the AIFI annotation models require manual
acquisition. Several source datasets are tens of gigabytes or larger.

For a step-by-step walkthrough, see [Getting started](docs/getting-started.md).

## Check your setup

Once the models and at least one study's inputs are available, create small
fixtures and run the test workflow:

```bash
make test-data
make run-all-test
```

These commands use the representative `aida25,cima,terekhova26` selection and
skip unavailable fixtures. Test outputs go to `output/test/`. They are software
checks, including synthetic data, and are not biological results.

## Run an analysis

Run the studies whose full inputs you have acquired. For example:

```bash
make run-all STUDIES=cima
```

This runs the full pipeline, including integration benchmarking, in one Nextflow
run. Benchmark methods run one at a time. On a constrained machine, use
`RUN_BENCHMARK=false` and run `make run-integration-benchmark` later. Age clocks
are optional: enable them with `RUN_CLOCKS=true`.

Replace `cima` with your study ID or a comma-separated list. This runs
harmonization, pseudobulk differential expression, focused gene profiles, and
cell-type reports, plus the integration benchmark. Completed production tasks can
be resumed on later runs. Use `make run` for harmonization only, or
`make run-analysis` for downstream analyses from existing merged files. See
[workflow modes](docs/development/orchestration.md#entry-points) for details.
The [input guide](docs/inputs/README.md) explains how to check that a local file
is a full study input rather than a development subset.

## Find your results

Open `output/index.html` in a browser to navigate the reports.

| Result | Location |
| --- | --- |
| Harmonized studies and merged count matrices | `output/harmonization/` |
| Cell-type reports and differential-expression results | `output/analysis/per_type/` |
| Cross-cell-type summaries | `output/analysis/combined/` |
| Focused gene reports | `output/analysis/per_gene/` |
| Selected and skipped inputs, software and provenance | `output/run_manifest.json` |

See [the report guide](docs/analysis/reports.md) for reading and rerendering
reports, and [the DE guide](docs/analysis/differential-expression.md) for model
assumptions and interpretation.

If a run reports failed or ignored tasks, use `make workflow-failures` to see
task errors, affected DE models, and log paths. See
[inspecting an incomplete run](docs/getting-started.md#inspect-an-incomplete-run)
for partial-result options and verification behavior.

## Explore the interactive atlas

After an analysis has finished:

```bash
make site-export
make serve-frontend
```

Visit `http://localhost:8000/`. Exporting replaces the local site data with a
validated release of the available results; serving the site does not run
analysis. See [the site guide](docs/site/README.md) for test data, static builds,
and publication.

## Learn more

- [Documentation index](docs/README.md) — guides and reference material by task.
- [Age clocks](docs/analysis/age-clocks.md) and [gene profiles](docs/analysis/gene-profiles.md) — optional analyses.
- [Adding a study](docs/development/studies.md) — adapters, configuration, and metadata evidence.
- [Development and validation](docs/development/README.md) — local setup and choosing a Make target.
- [Architecture](docs/reference/architecture.md), [output contracts](docs/reference/outputs.md), and [analysis settings](docs/reference/analysis-settings.md) — implementation reference.

`make help` lists the available commands.

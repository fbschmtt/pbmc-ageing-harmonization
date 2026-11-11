# Getting started

[Documentation](README.md)

This walkthrough takes you from an input-ready checkout to a small test run,
a production analysis, and the resulting reports. Run commands from the
repository root.

## 1. Prepare the environment

Install Docker, Make, Python 3, and Nextflow 25.04 or newer on the host. Make
builds the locked analysis images when they are needed. Analysis libraries,
RDS conversion, and report execution run in Docker.

A local scientific Python environment is only needed for development; see
[development setup](development/README.md). A local Node installation is
optional because the site also has Docker-backed Make targets.

## 2. Acquire inputs for your studies

The [input catalog](inputs/files.md) lists the required expression files,
auxiliary metadata, and AIFI CellTypist models. Follow the [acquisition
guide](inputs/README.md) to download public inputs and identify manual sources.
Downloads are explicit and can be large. Installation, tests, and analysis do
not download scientific inputs automatically.

The AIFI models require interactive Allen authentication and must be acquired
separately. Public and manual study inputs must be placed at their configured
paths. The optional AIFI color table supplies published label colors; reports
use stable fallback colors when it is absent.

Development workspaces can contain small subsets at production-looking paths.
Do not use them for scientific analysis. Check the acquisition guide and source
provenance before starting a production run.

## 3. Run a small test

After acquiring models and at least one study's inputs:

```bash
make test-data
make run-all-test
```

`test-data` creates missing small fixtures under `test_data/` without changing
`input_data/`. The default selection is `aida25,cima,terekhova26`; unavailable
fixtures are skipped. To select another set, pass the same `STUDIES` list to
both commands, for example `STUDIES=cima,onek1k`.

The ordered test covers the core workflow, integration benchmark, synthetic
DE, age clocks, gene profiles, and selected cell-type reports. Results appear
under `output/test/`; open `output/test/index.html` to inspect the reports.
Synthetic and test-mode results are software validation only. See
[workflow testing](development/testing.md) for coverage and limitations.

## 4. Run your analysis

Use an explicit list of studies with complete inputs:

```bash
make run-all STUDIES=cima
```

Replace `cima` with your selected IDs. `run-all` runs core harmonization,
pseudobulk DE, profiles of top age-associated genes, and cell-type reports.
Production stages resume completed tasks. For `make run`, an explicit list is
strict about missing inputs; `STUDIES=all` instead selects available studies and
records omissions. Review `output/run_manifest.json` to see what was included.

Open `output/index.html` for the report index. Count matrices live under
`output/harmonization/`, per-type results under `output/analysis/per_type/`,
and cross-type summaries under `output/analysis/combined/`. The
[report guide](analysis/reports.md) explains the views; the
[output reference](reference/outputs.md) lists individual artifacts.

## 5. Explore existing results in the atlas

```bash
make site-export
make serve-frontend
```

Open `http://localhost:8000/`. Exporting builds and validates a static release
from the current analysis outputs and replaces `site/public/data`. It does not
rerun analysis. The [site guide](site/README.md) covers alternate inputs,
static builds, and publication.

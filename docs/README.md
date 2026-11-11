# Documentation

Start with the [project README](../README.md) for prerequisites, a first run,
and where to find results. All commands below are run from the repository root;
paths shown in code refer to that root unless stated otherwise.

## Run the pipeline

- [Getting started](getting-started.md): acquire inputs, check your setup, run an analysis, and open the results.
- [Acquiring inputs](inputs/README.md): public downloads, manual sources, and integrity checks.
- [Input catalog](inputs/files.md): exact study paths, source links, sizes, checksums, and annotation models.

## Understand and explore results

- [Harmonization and annotation](analysis/harmonization.md): data flow and cell-label semantics.
- [Reports](analysis/reports.md): QC, cell-type views, published artifacts, and rerendering.
- [Differential expression](analysis/differential-expression.md): sample inclusion, covariates, gene universes, trajectories, abundance, and residual analyses.
- [Gene profiles](analysis/gene-profiles.md): inspect an individual gene or the top age-associated genes.
- [Age clocks](analysis/age-clocks.md): training, held-out evaluation, calibration, and synthetic experiments.
- [Interactive atlas](site/README.md): export, run, build, and publish the site.

## Develop or extend the pipeline

- [Development and Make targets](development/README.md): local setup and the smallest useful validation command.
- [Orchestration](development/orchestration.md): composed workflows, caching, failure isolation, and finalization.
- [Workflow testing](development/testing.md): test fixtures, synthetic coverage, and full verification.
- [Adding a study](development/studies.md): adapters, metadata joins, and provenance.

## Reference

- [Architecture](reference/architecture.md): components and their responsibilities.
- [Output contracts](reference/outputs.md): directory layout, H5AD contracts, and per-type analysis artifacts.
- [Metadata semantics](reference/metadata.md): samples, donors, batches, and missing values.
- [Analysis settings](reference/analysis-settings.md): configuration defaults and their meaning.
- [Atlas data contract](site/data-contract.md): asset schemas, provenance, and release validation.

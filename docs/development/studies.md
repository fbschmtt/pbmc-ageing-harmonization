# Adding a study

[Documentation](../README.md)

Start with [templates/study_adapter.py](../../templates/study_adapter.py) and
[templates/study_config.json](../../templates/study_config.json). Copy the adapter to
`src/pbmc_pipeline/studies/<study_id>.py`, replace its rules with the source's
cell and sample metadata, then add the configuration object under `studies` in
`config/studies.json`. Add every expression object and supplementary metadata
file to `config/input_sources.json` and `docs/inputs/files.md`; paths must be relative
to the repository root. Both registry documents must satisfy their tracked JSON
Schemas as well as the pipeline's cross-reference checks. Declare supplementary files in
`preparation.dependencies` so Nextflow stages them and configuration changes are
tracked in task identity.

Each study requires a `publication_doi`. When available, also provide
`publication_title` with the paper's full title and `publication_citation` with
its formal reference (authors, journal or preprint server, year, and available
volume, issue, and page details). These optional fields appear with the DOI in
the site's study cards; they credit the source authors and help readers identify
the underlying study. Keep the DOI as the canonical identifier even when title
or citation metadata is unavailable.

An adapter receives the expression object's `.obs` table and must return the
canonical metadata fields, with exactly one row for each retained cell. Use
`safe_left_join` for sample-level supplements: it rejects ambiguous joins and
preserves the expression-cell order. The usual route is to retain every input
cell. When an object deliberately mixes cohorts, filter it explicitly in the
adapter, explain the condition in `provenance.selection`, and set
`preparation.allow_cell_subset` to `true`. The preparation and harmonization
steps then verify that retained IDs are a strict subset of the source IDs and
drop the excluded cells from the count matrix before annotation. This prevents
metadata-only filtering from leaking excluded cells into downstream outputs.

Choose `counts_source` as `X`, `raw`, or `layer:<layer_name>` after checking
which matrix contains unnormalized integer counts. The harmonized output always
stores only those selected counts in `.X`; source `.raw` and every source layer
are removed. For CellTypist studies, harmonization writes counts to the output
H5AD before normalizing `.X` in memory. After annotation, it updates metadata
and embeddings in place while preserving the saved count matrix, then writes
the run report.

When an H5AD uses Ensembl IDs as `.var_names` but carries gene symbols in a
feature metadata column, configure `features.source_column` to that symbol
column (usually `feature_name`). The harmonizer repairs `.var_names` from this
column and sums duplicate symbols according to `duplicate_policy`. Named genes
use the configured source column; only missing or blank names fall back to the
original `.var_names` identifiers; unresolved identifiers cause harmonization
to fail rather than creating a shared placeholder gene. Fachrul26 and Perez22
follow this pattern. Harmonization reports record the source column, symbol-like
feature count, identifier-like feature count, missing names, and duplicate-name
counts so source annotation quality is visible rather than implicit.

Nehar-Belaid Treg per-study results contain both `C2ORF15` and `C2orf15` in the
sex model. These spellings are kept distinct: feature aggregation matches exact
symbol text, and the site exporter now matches exact capitalization when
counting study support. `C2orf15` is the [HGNC-approved symbol](https://www.ncbi.nlm.nih.gov/gene/150590).
The current handling prevents an export collision but does not combine the
features. If canonicalization is added later, resolve symbols through stable
gene IDs before aggregating counts; do not change capitalization globally.

Configuration loading also checks adapter module existence, count-source and
feature-repair policies, conversion declarations, and that every study input
and preparation dependency is represented in `config/input_sources.json`.

Set `technology` to the most specific reported assay information. For example,
use `10X3'v3` rather than `10X3'` when the reagent-kit revision is available;
use the assay-family value only when no more specific source information is
reported. Record the source field in `provenance`, and use the lowest-level
run, pool, library, or processing-group identifier that assigns one value per
sample and groups multiple samples where possible. Record
every inferred metadata value in `provenance` with its rationale; do not make
it look like a reported source field.

Before enabling the study in an all-study run, verify that the adapter's cell
selection, sample-level joins, count source, technology/provenance fields, and
input download declarations match the publication and supplements. Generate a
test fixture with
`make test-data`, then run `make run-test STUDIES=<study_id>` and inspect its
prepared-table report, harmonized H5AD, and QC report. Confirm that the output
has integer counts in `.X`, `raw is None`, and no layers before adding it to an
all-study run.

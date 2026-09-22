# Input data

This directory is the single entry point for immutable pipeline inputs. Data
files are ignored by Git; expected relative paths are declared in
`config/studies.json`. Each study owns a directory containing its expression
object and every study-specific supplementary file.

Downsampled smoke-test fixtures are kept separately under `test_data/`, so full
and test inputs can coexist. Create or refresh them with `make test-data`.
Wang25's test fixture is an RDS, so test runs exercise the same RDS conversion
process as production.

The current study directories are `aida25/`, `aifi/`, `fachrul26/`,
`nehar_belaid26/`, `onek1k/`, `perez22/`, `terekhova23/`, and `wang25/`. Keep
source filenames intact so they
remain recognizable against their accessions and publications. Generated
conversions belong in `cache/`, and smoke-test fixtures belong in `test_data/`;
neither is an immutable study input.

Fachrul26 and Perez22 carry the metadata needed by their adapters in the H5AD
`.obs` tables, so they have no supplementary runtime dependencies. Their source
`.var_names` are Ensembl IDs; the configuration uses each object's
`.var["feature_name"]` column as the gene-symbol source and sums duplicate
symbols during harmonization.

Nehar-Belaid26 uses Supplementary Data 1, sheet `1a`, as the authoritative
metadata for the in-house scRNA-seq cohort. Its all-PBMC H5AD also includes
reused public studies; the configured adapter retains only cells whose `Study`
value is `Nehar-Belaid_et_al`.

The Terekhova compact lookup is the canonical runtime input mapping each cell ID
to its demultiplexed `Tube_id`. It was generated from the published
`all_pbmcs_metadata.csv`, which is not reconstructible from `mmc2.xlsx` because
the workbook contains no cell barcodes. If the large source CSV is downloaded
again, recreate the lookup with:

```bash
python scripts/prepare_auxiliary_metadata.py
```

`config/input_sources.json` is the tracked input manifest. It records direct
URLs for the AIFI expression object and the Nehar-Belaid26 H5AD plus
Supplementary Data 1, checksums for local supplementary files, and manual
acquisition locations for the remaining inputs. Populate verified direct URLs
and expression checksums before treating the pipeline as release-ready.

For convenience, `make download-inputs STUDIES=<study>` reads the tracked
`config/input_sources.json` manifest and makes a best-effort download of direct
public files. It is opt-in and deliberately does not authenticate or scrape
portals. Manual-only sources are reported with their acquisition instructions;
the ignored `download_manifest.json` records what was obtained locally.

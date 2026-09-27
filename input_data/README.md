# Input data

This directory is the single entry point for immutable pipeline inputs. Data
files are ignored by Git; expected relative paths are declared in
`config/studies.json` and recorded in `config/input_sources.json`. Their tracked
Schemas and cross-reference checks require every configured expression input
and preparation dependency to have a manifest entry. Each study owns a
directory containing its expression object and every study-specific
supplementary file.

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
metadata for the in-house scRNA-seq cohort. The configured input is the GEO RAW tar, whose compressed 10x Matrix Market members are read directly without extraction. All complete GEX 10x triplets are discovered from archive member names; their `JB` subject IDs join to Supplementary Data 1a `IDs`, and its `Names` values become `sample_id`. In the published label H5AD, `obs["sample_id"]` is the sample identity; its cell index is not used for that purpose. Labels transfer only when `sample_id` and the A/C/G/T barcode sequence uniquely agree. This intentionally ignores the numeric 10x suffix because the published H5AD reassigned it during concatenation. Among the retained H5AD Nehar-Belaid cells, all 88 sample IDs occur in the supplement; the supplement has seven additional samples (`HC17`, `HI3`, `HI26`, `HI27`, `HO13`, `HO21`, `HY19`) with no retained label-H5AD cells, so their labels remain `not_provided`.

The Terekhova compact lookup is the canonical runtime input mapping each cell ID
to its demultiplexed `Tube_id`. It was generated from the published
`all_pbmcs_metadata.csv`, which is not reconstructible from `mmc2.xlsx` because
the workbook contains no cell barcodes. If the large source CSV is downloaded
again, recreate the lookup with:

```bash
python scripts/prepare_auxiliary_metadata.py
```

`config/input_sources.json` is the tracked input manifest. It records direct
URLs for the AIFI expression object and the Nehar-Belaid26 RAW archive, label H5AD, plus
Supplementary Data 1, checksums for local supplementary files, and manual
acquisition locations for the remaining inputs. Populate verified direct URLs
and expression checksums before treating the pipeline as release-ready.

For convenience, `make download-inputs STUDIES=<study>` reads the tracked
`config/input_sources.json` manifest and makes a best-effort download of direct
public files. It is opt-in and deliberately does not authenticate or scrape
portals. Manual-only sources are reported with their acquisition instructions;
the ignored `download_manifest.json` records what was obtained locally.

# Input data

This directory is the single entry point for immutable pipeline inputs. Data
files are ignored by Git; expected relative paths are declared in
`config/studies.json`.

Downsampled smoke-test fixtures are kept separately under `test_data/`, so full
and test inputs can coexist. Create or refresh them with `make test-data`.
Wang25's test fixture is an RDS, so test runs exercise the same RDS conversion
process as production.

Expression inputs currently use the top-level paths established by the original
conversion notebook. Supplementary inputs are grouped by study:

- `metadata/aida25/mmc1.xlsx`
- `metadata/terekhova23/mmc2.xlsx`
- `metadata/terekhova23/cell_to_tube.csv.gz`
- `metadata/wang25/41590_2024_2059_MOESM3_ESM.xlsx`

The Terekhova compact lookup is the canonical runtime input mapping each cell ID
to its demultiplexed `Tube_id`. It was generated from the published
`all_pbmcs_metadata.csv`, which is not reconstructible from `mmc2.xlsx` because
the workbook contains no cell barcodes. If the large source CSV is downloaded
again, recreate the lookup with:

```bash
python scripts/prepare_auxiliary_metadata.py
```

`config/input_sources.json` is the tracked input manifest. It currently records
one direct AIFI expression URL, checksums for local supplementary files, and
manual acquisition locations for the remaining inputs. Populate verified direct
URLs and expression checksums before treating the pipeline as release-ready.

For convenience, `make download-inputs STUDIES=<study>` reads the tracked
`config/input_sources.json` manifest and makes a best-effort download of direct
public files. It is opt-in and deliberately does not authenticate or scrape
portals. Manual-only sources are reported with their acquisition instructions;
the ignored `download_manifest.json` records what was obtained locally.

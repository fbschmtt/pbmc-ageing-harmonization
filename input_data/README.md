# Input data

This directory is the single entry point for immutable pipeline inputs. Data
files are ignored by Git; expected relative paths are declared in
`config/studies.json`.

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

Source URLs, versions, and checksums should be added to a tracked input manifest
before the pipeline is considered release-ready.

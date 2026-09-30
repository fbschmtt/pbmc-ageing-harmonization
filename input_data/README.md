# Input data

This directory is the single entry point for immutable pipeline inputs. Data
files are ignored by Git; expected relative paths are declared in
`config/studies.json` and recorded in `config/input_sources.json`. Their tracked
Schemas and cross-reference checks require every configured expression input
and preparation dependency to have a manifest entry. Each study owns a
directory containing its expression object and every study-specific
supplementary file.

Downsampled smoke-test fixtures are kept separately under `test_data/` and are
created or refreshed with `make test-data`. The current ignored files at several
production input paths in this checkout are also smoke-test-scale subsets
(200 or 500 cells), not the complete published sources. The four CELLxGENE H5AD
objects have authoritative expected byte sizes in the input registry, so the
optional downloader reports these local subsets as `size_mismatch`. Do not use
them as full production inputs. To replace a selected local file with its
published download, opt in explicitly, for example
`make download-inputs STUDIES=onek1k DOWNLOAD_FORCE=true`; this downloads the
multi-gigabyte source object. Wang25's test fixture is an RDS, so test runs
exercise the same RDS conversion process as production.
This describes the ignored files in the development checkout; it does not
describe or invalidate the full inputs used by the successful production run.

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

## Acquisition

`config/input_sources.json` is the tracked source registry. It records direct
URLs where programmatic retrieval is appropriate, integrity information where
available, and manual-source instructions where access conditions require it.
Acquisition is strictly opt-in and is documented once in the top-level
[Download inputs guide](../README.md#download-inputs); neither installation nor
workflow execution downloads data. Consult [INPUT_FILES.md](../INPUT_FILES.md)
for the exact paths, links, checksums, and reasons that a source is manual.

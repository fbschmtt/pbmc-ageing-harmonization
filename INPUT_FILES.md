# External input files

The pipeline source does not redistribute expression data, supplementary study
metadata, or trained CellTypist models. Place each external file at the exact
path below. Paths are relative to the repository root.

Checksums are recorded where the file is currently available locally. A future
download script should verify these values before starting the pipeline.

## Expression data

| Study | Required path | Role | SHA-256 |
|---|---|---|---|
| AIDA25 | `input_data/3a8a77a6-f069-479c-a2c1-252feafd0106.h5ad` | Raw-count H5AD | Not yet recorded |
| AIFI | `input_data/human_immune_health_atlas_full.h5ad` | Immune Health Atlas H5AD | Not yet recorded |
| OneK1K | `input_data/81d84489-bff9-4fb6-b0ee-78348126eada.h5ad` | Raw-count H5AD | Not yet recorded |
| Terekhova23 | `input_data/pbmc_gex_raw_with_var_obs.h5ad` | Raw-count H5AD | Not yet recorded |
| Wang25 | `input_data/scRNA-seqProcessedLabelledObject.rds` | Seurat RDS converted by `scripts/convert_rds.R` | Not yet recorded |

The Wang RDS conversion produces
`cache/converted/scRNA-seqProcessedLabelledObject.h5ad`. This is a generated
intermediate, not an external input.

## Supplementary metadata

| Study | Required path | Role | SHA-256 |
|---|---|---|---|
| AIDA25 | `input_data/metadata/aida25/mmc1.xlsx` | Donor and single-cell batch metadata | `7610edd1105181e7fecb9e61e25e2668d57433524cc801838ed9c76f92bb4eef` |
| Terekhova23 | `input_data/metadata/terekhova23/mmc2.xlsx` | Donor, visit, BMI, and ethnicity metadata | `af8da63efbaf6800da79f4660412ea3100aac88b18748b01679593e29b4444a2` |
| Terekhova23 | `input_data/metadata/terekhova23/cell_to_tube.csv.gz` | Cell ID to demultiplexed tube mapping | `bf0dede2434136b6f5515210c6c48fb2a7105c60e836ad007b7407856705df5a` |
| Wang25 | `input_data/metadata/wang25/41590_2024_2059_MOESM3_ESM.xlsx` | Published sample-level metadata | `a46e3aa3f877890fa6100e42a6ffc641ec481e7f122b00b5a68fa89689b9fd2f` |

`cell_to_tube.csv.gz` is the compact canonical Terekhova runtime input. It was
derived from the much larger published `all_pbmcs_metadata.csv` by
`scripts/prepare_auxiliary_metadata.py`. It cannot be reconstructed from
`mmc2.xlsx`, because that workbook contains no cell barcodes. To recreate the
compact lookup, first place the original file at
`input_data/metadata/terekhova23/all_pbmcs_metadata.csv` and run the preparation
script.

## CellTypist models

The models were downloaded from the Allen Institute Immune Health Atlas model
downloads page on 2025-06-27. Confirm redistribution terms before publishing
them; the binary files are ignored by Git.

| Level | Required path | SHA-256 |
|---|---|---|
| AIFI L1 | `aifi_models/ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl` | `c63ec36ec90e195706dbab0e5596082d3bc06d980fcb212781c60c5ddb63d973` |
| AIFI L2 | `aifi_models/ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl` | `1f5d49bb85bf1fdff4936a6b250a4a24362ec5b58472a73d891f77f706e0eb10` |
| AIFI L3 | `aifi_models/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl` | `df18d419fc3b73468b5d0e3fd940871c3930f207bddab91f1fba69e3905ff5cc` |

Source page:
<https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/>

## Test data

Files under `cache/test_inputs/` are deterministic 200-cell subsets produced from the
full expression inputs by `scripts/create_test_data.py`. They are generated test
artifacts rather than external dependencies and are ignored by Git.

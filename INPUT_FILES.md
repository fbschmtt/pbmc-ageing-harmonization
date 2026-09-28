# External input files

The pipeline source does not redistribute expression data, supplementary study
metadata, or trained CellTypist models. Place each external file at the exact
path below. Paths are relative to the repository root.

Checksums are recorded where the file is currently available locally. The
best-effort `make download-inputs` command reads `config/input_sources.json`;
it checks a SHA-256 where recorded, checks expected byte sizes for CELLxGENE
objects, and reports manual-only sources. Full published CELLxGENE H5ADs are
large; local files at these paths may be test-scale subsets and must be checked
before production use. To deliberately replace existing files, pass
`DOWNLOAD_FORCE=true` to the Make target.

## Expression data

| Study | Required path | Role and source | SHA-256 / expected size |
|---|---|---|---|
| AIDA25 | `input_data/aida25/3a8a77a6-f069-479c-a2c1-252feafd0106.h5ad` | [CELLxGENE AIDA Phase 1 Data Freeze v2](https://cellxgene.cziscience.com/collections/ced320a1-29f3-47c1-a735-513c7084d508), raw-count H5AD | 14,266,895,920 bytes expected |
| AIFI | `input_data/aifi/human_immune_health_atlas_full.h5ad` | [Allen Institute Immune Health Atlas download](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/scrna/), H5AD | Not yet recorded |
| Fachrul26 | `input_data/fachrul26/6a322de5-4cc6-43f5-b35d-6f0c246fb297.h5ad` | [CELLxGENE dataset](https://cellxgene.cziscience.com/collections/d1e0e64d-6d2a-4a3e-b7f4-43ed909a9d9c), raw-count H5AD with embedded observation metadata | 3,855,510,491 bytes expected |
| Nehar-Belaid26 | `input_data/nehar_belaid26/GSE233321_RAW.tar` | [Direct GEO RAW archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233321/suppl/GSE233321_RAW.tar); 10x members are read directly from tar | Not yet recorded |
| Nehar-Belaid26 | `input_data/nehar_belaid26/GSE233321_all_PBMCs.h5ad` | [GEO label H5AD](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233321); labels transfer only on unique `obs.sample_id` plus normalized barcode-sequence matches | Not yet recorded |
| OneK1K | `input_data/onek1k/81d84489-bff9-4fb6-b0ee-78348126eada.h5ad` | [CELLxGENE OneK1K dataset](https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1), raw-count H5AD | 4,436,424,475 bytes expected |
| Perez22 | `input_data/perez22/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad` | [CELLxGENE dataset](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2), raw-count H5AD with embedded observation metadata | 12,218,105,530 bytes expected |
| Terekhova23 | `input_data/terekhova23/pbmc_gex_raw_with_var_obs.h5ad` | [Synapse study data](https://www.synapse.org/Synapse:syn49637038); access conditions may apply | Not yet recorded |
| Wang25 | `input_data/wang25/scRNA-seqProcessedLabelledObject.rds` | [Synapse entity syn61609846](https://www.synapse.org/Synapse:syn61609846), Seurat RDS converted by `scripts/convert_rds.R` | Not yet recorded |

The Wang RDS conversion produces
`cache/converted/scRNA-seqProcessedLabelledObject.h5ad`. This is a generated
intermediate, not an external input.

## Supplementary metadata

| Study | Required path | Role | SHA-256 |
|---|---|---|---|
| AIDA25 | `input_data/aida25/mmc1.xlsx` | Donor and single-cell batch metadata | `7610edd1105181e7fecb9e61e25e2668d57433524cc801838ed9c76f92bb4eef` |
| Terekhova23 | `input_data/terekhova23/mmc2.xlsx` | Donor, visit, BMI, and ethnicity metadata | `af8da63efbaf6800da79f4660412ea3100aac88b18748b01679593e29b4444a2` |
| Terekhova23 | `input_data/terekhova23/cell_to_tube.csv.gz` | Cell ID to demultiplexed tube mapping | `bf0dede2434136b6f5515210c6c48fb2a7105c60e836ad007b7407856705df5a` |
| Nehar-Belaid26 | `input_data/nehar_belaid26/41467_2026_73729_MOESM3_ESM.xls` | Nehar-Belaid et al. (2026) Supplementary Data 1; sheet 1a in-house donor metadata | `c55555c2765b0fb4cbe872ffdd37271dde4d0aefa5a92d3343d2368778f8b0b5` |
| Wang25 | `input_data/wang25/41590_2024_2059_MOESM3_ESM.xlsx` | Published sample-level metadata | `a46e3aa3f877890fa6100e42a6ffc641ec481e7f122b00b5a68fa89689b9fd2f` |
| Terekhova23 | `input_data/terekhova23/all_pbmcs_metadata.csv` | [Synapse study source](https://www.synapse.org/Synapse:syn49637038); large source metadata used to generate the compact lookup below, not read by the normal workflow | Not yet recorded |

`cell_to_tube.csv.gz` is the compact canonical Terekhova runtime input. It was
derived from the much larger published `all_pbmcs_metadata.csv` by
`scripts/prepare_auxiliary_metadata.py`. It cannot be reconstructed from
`mmc2.xlsx`, because that workbook contains no cell barcodes. To recreate the
compact lookup, first place the original file at
`input_data/terekhova23/all_pbmcs_metadata.csv` and run the preparation script.

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

`make test-data` creates deterministic 200-cell fixtures under `test_data/` from
the production expression inputs. Python-source studies are first read through
their configured input reader and then sampled by `scripts/create_test_data.py`;
thus Nehar-Belaid26 fixtures exercise direct tar reading during fixture creation,
but test workflow runs read the resulting H5AD fixture. Each configured RDS
conversion study is sampled by `scripts/create_test_rds.R`. These generated
fixtures are ignored by Git and are used only by the Nextflow `test` profile.

# Sample and metadata semantics

[Documentation](../README.md)

`subject` is a biological individual. `sample` is one biological specimen;
it is the unit aggregated into sample × AIFI-L2
pseudobulks. `batch_single_cell` records the lowest available technical
grouping above samples that can be assigned consistently: each sample maps to
one value, and a useful grouping is shared by multiple samples. Examples
include multiplexed 10x processing runs and processing batches. A sample may
not map to multiple technical batches. If a source value is missing or
ambiguous for a sample, that sample receives a sample-specific fallback. If a
source only provides one study-wide value or a unique value per sample,
samples are used as their own batches so the field does not simply mirror the
study. The preparation report records how many batch groups and sample-specific
fallbacks were written. Study prefixes keep source batch IDs distinct in merged
data. This field describes technical correlation groups; it never defines a
pseudobulk. Missing string metadata in other fields is written as
`not_provided`; missing numeric metadata is `NaN`.

`study_celltype` preserves the highest-resolution published or author-provided
cell-type annotation available per cell. An unsupervised cluster is not used
merely because it has more values. The selected fields below are source labels
for audit and reporting; cross-study analyses use the separately harmonized
AIFI labels.

Every harmonized and merged H5AD also contains `obs.study_doi`. It is a
categorical per-cell DOI drawn from the required `publication_doi` study
configuration field; it is intentionally not duplicated in `uns`.

| Study | `study_celltype` source | Selection basis |
|---|---|---|
| AIDA25 | `Annotation_Level4` | Highest level of the published AIDA hierarchy. |
| AIFI | `AIFI_L3` | Highest AIFI level supplied in the source H5AD. |
| CIMA | `cell_type_l4` | Highest CIMA hierarchy level. |
| Fachrul26 | `cell_type` | One-to-one with `author_cell_type`; no finer author annotation is supplied. |
| Nehar-Belaid26 | `LS_L4` | Highest published LS hierarchy level. |
| OneK1K | `cell_type` | One-to-one with the supplied `predicted.celltype.l2`; no finer level is supplied. |
| Perez22 | `cell_type` | One-to-one with `author_cell_type`; no finer author annotation is supplied. |
| Terekhova23 | `Cluster_names` | Only supplied per-cell type-label field. |
| Terekhova26 | `author_celltype` | Highest author label; identical to `annotation`. |
| Wang25 | `secondary_type` | Finer of the supplied `primary_type`/`secondary_type` fields. |

| Study          | Reference DOI                                                                | Subject                     | Sample                        | Technical batch                  | Notes                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| -------------- | ---------------------------------------------------------------------------- | --------------------------- | ----------------------------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| AIDA25         | [10.1016/j.cell.2025.02.017](https://doi.org/10.1016/j.cell.2025.02.017)     | `donor_id`                  | country-qualified donor ID    | supplementary experimental batch | Source `obs.assay` reports 10x 5′ v2. Lonza controls (`DCP_ID` prefixed `LONZA`) are reported as European and retained in harmonized and pseudobulk data for calibration, but excluded from DE models by subject ID. Repeated donor-country records are deduplicated in the supplementary metadata join and their batch is marked unclear. BMI uses the supplementary donor metadata.                                                                                                                                                                                                                                                        |
| AIFI           | [10.1038/s41586-025-09686-5](https://doi.org/10.1038/s41586-025-09686-5)     | specimen GUID prefix        | specimen GUID                 | `pool_id` (`aifi_` prefix)       | The configured Immune Health Atlas reference has one sample per donor. Each `pool_id` combines hashed samples and is distributed across chips and wells; a sample maps to one pool, while it can span wells. We use this lowest shared grouping for technical batch. A broader `batch_id` groups the processing occasion (up to 24 samples across multiple pools), not an exact timestamp; when the same aliquots are processed for scRNA-seq and flow cytometry on the same day, they share a Batch ID. The source fixture confirms one pool per specimen and multiple specimens per pool. Although the source includes `sample.visitName`, it is not retained as a separate harmonized timepoint field. The source H5AD has no assay field. ([AIFI pipeline](https://apps.allenimmunology.org/aifi/pipelines/multiplexed-scrna-seq/)) |
| CIMA           | [10.1126/science.adt3130](https://doi.org/10.1126/science.adt3130)           | `obs.sample`                | `obs.sample`                  | `obs.RNA_experiment`             | The H5AD embeds the retained sample and experiment fields; `RNA_experiment` groups samples by date-coded RNA processing experiment, while `library` can vary within a sample. Each of the 408 retained sample IDs maps to exactly one of 51 RNA experiment groups in local metadata. The separately supplied metadata table is deliberately not a workflow input. Native `cell_type_l4` (70 labels) is retained as the per-study label. Paper methods report DNBelab C Series, STAR v2.7.2b, GRCh38, and exclusion of intronic reads; no demultiplexing method is described, so `none` is inferred. Frozen status is not reported.                                                            |
| OneK1K         | [10.1126/science.abf3041](https://doi.org/10.1126/science.abf3041)           | `donor_id`                  | `donor_id`                    | `pool_number`                    | Donors are the biological sample-level replicates. `pool_number` records the multiplexed 10x processing pool, which contains cells from multiple donors, and is retained as a technical batch. The source has no separate specimen or timepoint field.                                                                                                                                                                                   |
| Terekhova23    | [10.1016/j.immuni.2023.10.013](https://doi.org/10.1016/j.immuni.2023.10.013) | `Donor_id`                  | `Tube_id`                     | barcode-derived batch            | Only each donor's earliest available visit is retained before harmonization; later-visit cells are excluded. Tubes spanning multiple barcode-derived batches receive sample-specific fallback values to preserve one batch assignment per tube. The paper explicitly reports frozen input material.                                                                                                                                                                                                                                    |
| Terekhova26    | [10.1016/j.immuni.2026.09.012](https://doi.org/10.1016/j.immuni.2026.09.012) | `obs.donor_id`              | `obs.donor_id`                | `obs.Pool`                       | The supplied author-reannotated KCL H5AD embeds all runtime metadata. Native `author_celltype` (59 labels) is retained as the per-study label. Genetic demultiplexing (`souporcell_donor`) is paired with HTO antibodies that map genotypes to samples. Paper methods report cryopreserved/thawed PBMCs, Cell Ranger multi v9.0.1, and GRCh38-2020-A; intronic reads are inferred included from the Cell Ranger v6+ default arguments. |
| Wang25         | [10.1038/s41590-024-02059-6](https://doi.org/10.1038/s41590-024-02059-6)     | `Sample ID`                 | `Sample ID`                   | sample-specific fallback        | No separate donor or technical-batch identifier has yet been recovered, so each sample is its own batch. Fresh/frozen status is not provided in the paper or supplementary metadata; the author was contacted but did not respond.                                                                                                                                                                                                                           |
| Nehar-Belaid26 | [10.1038/s41467-026-73729-2](https://doi.org/10.1038/s41467-026-73729-2)     | Supplementary Data 1a `IDs` | Supplementary Data 1a `Names` | Supplementary Data 1a `runs_10x` | Reads the supplied GEO raw 10x tar directly (without extraction); `JB` subject IDs parsed from member names join to Supplementary Data 1a `IDs`, which supplies sample IDs and donor metadata. Published labels transfer only on a unique sample ID plus normalized barcode-sequence match.                                                                                                                                            |
| Fachrul26      | [10.64898/2026.02.15.704933](https://doi.org/10.64898/2026.02.15.704933)     | `donor_id`                  | `sample_id`                   | `library_id`                     | Paper methods state PBMCs were cryopreserved and thawed; this takes precedence over the conflicting CELLxGENE `obs.sample_preservation_method=fresh`.                                                                                                                                                                                                                                                                                  |
| Perez22        | [10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970)           | `donor_id`                  | `sample_uuid`                 | `library_uuid`                   | Paper reports USA frozen PBMCs processed with 10x 3' v2 and Cell Ranger (version unspecified); intronic-read handling is not reported. Uses embedded H5AD metadata; `feature_name` supplies gene symbols. `obs.disease` is normalized and only `normal`/healthy cells are retained; embedded SLE cells are excluded before harmonization.                                                                                                  |

`tests/` contains only automated test code. Reproducible local intermediates
are collected under `cache/` for RDS-to-H5AD conversions. Downsampled smoke-test
fixtures are generated under the ignored `test_data/` directory.

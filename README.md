# PBMC ageing scRNA-seq pipeline

Configuration-driven processing of five PBMC ageing scRNA-seq studies. The
current milestone produces one harmonized H5AD per study with raw counts,
homogeneous metadata, and AIFI L1/L2/L3 cell-type labels. Merging and pseudobulk
analysis are deliberately deferred.

Study-specific decisions live in `config/studies.json`. The output metadata
contract lives in `config/harmonized_obs_schema.json`. See `PLAN.md` for scope,
progress, and unresolved scientific questions.

See `INPUT_FILES.md` for the complete list of external expression files,
supplementary metadata, CellTypist models, expected paths, and checksums.

## Current status

- The complete five-study workflow has run successfully on the 200-cell test
  inputs using the local `.venv`.
- Executed QC notebooks can produce self-contained per-study HTML reports.
- The full datasets have not been run in this workspace.
- The Python harmonization Dockerfile exists, but the image has not yet been
  built or executed.
- RDS conversion is now a parameterized R command. Its R container and the
  Nextflow orchestration layer are still pending.

## Data flow

```text
input_data/ expression objects + metadata
                │
                ├── Seurat RDS ──> cache/converted/*.h5ad
                │
                ▼
      harmonize one configured study
                │
                ├──> output/harmonized/<study>.h5ad
                └──> output/reports/<study>.json
```

All immutable local inputs belong under `input_data/`. Generated conversions,
reports, and harmonized outputs are kept outside it and ignored by Git.

`tests/` contains only automated test code. Reproducible local intermediates are
collected under `cache/`: downsampled fixtures in `cache/test_inputs/` and
RDS-to-H5AD conversions in `cache/converted/`.

## Local Python setup

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[dev,qc]'
```

## Run harmonization

```bash
# Fast validation using the downsampled inputs, without CellTypist
python scripts/harmonize_study.py --study all --test --validate-only

# Complete one study, including AIFI CellTypist labels
python scripts/harmonize_study.py --study onek1k

# Complete all test studies and generate executed notebook/HTML QC reports
python scripts/harmonize_study.py --study all --test --qc

# Recreate downsampled inputs from locally available full inputs
python scripts/create_test_data.py --cells 200
```

## Convert a Seurat RDS

The converter requires R, Seurat, and `sceasy`. It processes exactly one object
and retains only raw counts plus observation metadata:

```bash
Rscript scripts/convert_rds.R \
  --input input_data/scRNA-seqProcessedLabelledObject.rds \
  --output cache/converted/scRNA-seqProcessedLabelledObject.h5ad
```

Nextflow will eventually invoke this only for studies with a `conversion` entry
in `config/studies.json`.

## QC reports

Passing `--qc` executes the tracked `reports/qc_report.ipynb` template after each
harmonized H5AD is written. Each study gets a separate directory:

```text
output/qc/<study>/
├── executed.ipynb
└── report.html
```

The report includes provenance, metadata completeness, cohort
composition, count-depth checks, UMAPs, canonical PBMC marker expression, and
study-label/AIFI-label comparisons. All plots use the complete harmonized study.

To generate a report for an existing harmonized output without rerunning the
study:

```bash
python scripts/generate_qc_report.py \
  --input output/harmonized/onek1k.test.h5ad \
  --run-report output/reports/onek1k.test.json \
  --output-dir output/qc/onek1k.test \
  --study onek1k
```

## Docker

The current `Dockerfile` covers Python harmonization only. Once built, full input
and output directories must be mounted at the paths expected by the
configuration. A tested invocation will be added after the command-line path
overrides, reporting support, and R conversion image are implemented; Docker
execution should not yet be considered complete.

See `input_data/README.md` for the expected input layout and `PLAN.md` for the
remaining reproducibility work.

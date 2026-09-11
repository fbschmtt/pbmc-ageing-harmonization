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
- A DSL2 Nextflow workflow connects conversion, harmonization, and QC, with a
  test profile for the 200-cell inputs.
- Both Docker images build successfully. A fresh containerized Wang25 test run
  completed harmonization and QC, and a synthetic Seurat object passed the RDS
  to H5AD conversion with counts, names, and metadata intact.
- The full datasets and the real Wang25 RDS conversion have not yet been run.

## Data flow

```text
input_data/ expression objects + metadata
                │
                ├── Seurat RDS ──> converted H5AD
                │
                ▼
      harmonize one configured study
                │
                ├──> output/harmonized/<study>.h5ad
                ├──> output/reports/<study>.json
                └──> output/qc/<study>/{executed.ipynb,report.html}
```

All immutable local inputs belong under `input_data/`. Generated conversions,
reports, and harmonized outputs are kept outside it and ignored by Git.

`tests/` contains only automated test code. Reproducible local intermediates
are collected under `cache/` for RDS-to-H5AD conversions. Downsampled smoke-test
fixtures are generated under the ignored `test_data/` directory.

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

The installed commands also accept explicit artifact paths (`--input`,
`--output`, `--report-output`, and the QC command's `--template`). These are
used by Nextflow and are useful when debugging one task outside the workflow.
They write timestamped progress logs to standard error and retain their final
JSON/status output on standard output; Nextflow captures task logs in
`.command.err` and `.command.out` respectively.

## Convert a Seurat RDS

The converter requires R, Seurat, reticulate, and Python anndata. It processes exactly one object
and retains only raw counts plus observation metadata:

```bash
Rscript scripts/convert_rds.R \
  --input input_data/scRNA-seqProcessedLabelledObject.rds \
  --output cache/converted/scRNA-seqProcessedLabelledObject.h5ad
```

Nextflow invokes this for studies with a `conversion` entry in
`config/studies.json`. At present that is Wang25. Test runs use its downsampled
RDS fixture and therefore exercise the same conversion process.

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

## Nextflow and Docker

Nextflow >=25.04 is the canonical workflow runner. It reads studies and inputs from
`config/studies.json`, runs independent studies in parallel, caches successful
processes, and publishes the final H5AD, JSON, and QC artifacts. Limit a run
with a comma-separated study list:

```bash
nextflow run main.nf -profile standard,test --studies wang25 -resume
nextflow run main.nf -profile standard --studies aifi,onek1k -resume
nextflow run main.nf -profile docker,test --studies wang25 -resume
```

The `test` profile lowers resource requests and uses isolated 200-cell fixtures
from `test_data/`. RDS-based studies retain an RDS fixture and therefore still
pass through the same `CONVERT_RDS` process as a production run. Create or refresh fixtures
from the production inputs once with `make test-data`; this never modifies
`input_data/`. Override the fixture size with `make test-data TEST_CELLS=500`
when needed. `--skip_qc` stops after harmonization. Nextflow's `work/`
directory is the cache for converted and other intermediate artifacts; only
final deliverables are copied to `output/`.

Build and use both containers with:

```bash
make docker-build
make test-data
make pipeline-test STUDIES=wang25
make pipeline STUDIES=all
```

`make pipeline`, `make pipeline-test`, and `make pipeline-harmonize` build the
two local images automatically whenever `NF_PROFILE` contains `docker`; use
`make docker-images` to build them without launching a workflow. The image tags
(`pbmc-ageing-python:local` and `pbmc-ageing-r:local`) are intentionally local
Docker tags, not images available from a registry.

For a production run whose repository and task filesystem are separate, keep
the repository path readable by Docker and choose an external work/output
location:

```bash
make pipeline NF_PROFILE=docker STUDIES=all \
  WORK_DIR=/path/to/run/work OUTDIR=/path/to/run/output
```

The Docker profile read-only mounts the repository into each task. This lets the
task use the configured project root for metadata, models, and configuration
while keeping all generated task data in `WORK_DIR`. Do not place `WORK_DIR`
inside the repository when using this profile: the repository is mounted
read-only. Input staging retains Nextflow's default symbolic-link behavior, so
use a work filesystem that supports symbolic links.

The Python image is defined by `docker/python.Dockerfile`; the conversion image is defined by
`docker/r-conversion.Dockerfile`. The latter starts from a digest-pinned Seurat
image and adds the Python/anndata bridge. Metadata and externally
supplied CellTypist models remain outside the images and are declared as
Nextflow task dependencies. Expression files and generated artifacts are staged
through Nextflow. The Python image and a complete Wang25 test workflow have
been verified. The R conversion image is defined separately because it is much
larger and only needed for RDS inputs.

To preserve the Python dependency-install layer during source-only changes, the
Python image temporarily installs a minimal placeholder package, then copies
the real `src/` tree and reinstalls it with `--no-deps`. The cleanup before the
second install is required to prevent stale setuptools build artifacts.

During development, note that Nextflow identifies a container by its configured
image reference. Rebuilding a mutable `:local` tag does not invalidate an
existing `-resume` cache entry; run once without `-resume` after rebuilding, or
use a new versioned image tag.

## Make shortcuts

The `Makefile` is intentionally a thin convenience layer rather than a second
workflow. Run `make help` to see its targets. Common commands include:

```bash
make lint
make workflow-lint
make test
make validate-test STUDIES=wang25
make validate-full STUDIES=terekhova23
make pipeline-test STUDIES=wang25
make pipeline-harmonize STUDIES=wang25
```

For a local Nextflow run, activate `.venv` first so the installed
`pbmc-harmonize` and `pbmc-qc` commands are on `PATH`. The Docker profile does
not require the Python environment on the host.

See `input_data/README.md` for the expected input layout and `PLAN.md` for the
remaining reproducibility work.

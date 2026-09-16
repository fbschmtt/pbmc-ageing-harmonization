# PBMC ageing scRNA-seq pipeline

Configuration-driven processing of five PBMC ageing scRNA-seq studies. It produces
one harmonized H5AD per study, per-study sample × AIFI-L2 pseudobulks, and a
cross-study pseudobulk matrix. An explicitly enabled pathway also
creates a jointly embedded, freshly AIFI-L2-annotated single-cell merge.

Study-specific decisions live in `config/studies.json`. The output metadata
contract lives in `config/harmonized_obs_schema.json`. See `PLAN.md` for scope,
progress, and unresolved scientific questions. See `IMPLEMENTATION.md` for the
high-level architecture.

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
- Pseudobulk aggregation and merge are part of the production Nextflow graph;
  merge QC reports include study/sample contribution, AIFI-L2 overlap, depth,
  gene coverage, and (when applicable) merged UMAP checks.
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
                │
                ▼
  pseudobulk each study by sample/AIFI-L2
                │
                └──> output/pseudobulk/<study>.pseudobulk.h5ad
                ▼
       output/merged/pseudobulk_merged.h5ad
                │
                ├── optional: output/merged/single_cell_merged.h5ad
                └──> output/qc/merged/{executed.ipynb,report.html}
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
pbmc-harmonize --study all --test --validate-only

# Complete one study, including AIFI CellTypist labels
pbmc-harmonize --study onek1k

# Complete all test studies and generate executed notebook/HTML QC reports
pbmc-harmonize --study all --test --qc

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
pbmc-qc \
  --input output/harmonized/onek1k.test.h5ad \
  --run-report output/reports/onek1k.test.json \
  --output-dir output/qc/onek1k.test \
  --study onek1k
```

The merge workflow renders one combined document at `output/qc/merged/`. It
starts with pseudobulk composition and a gene-presence UpSet plot; if the
single-cell branch is enabled it appends its embedding, depth checks, and the
matrix comparing newly merged AIFI L2 calls to each cell's original study-level
call. `output/run_manifest.json` records selected studies, image references,
configuration checksum, and checksums of the merged deliverables.

## Nextflow and Docker

Nextflow >=25.04 is the canonical workflow runner. It reads studies and inputs from
`config/studies.json`, runs independent studies in parallel, caches successful
processes, and publishes the final H5AD, JSON, and QC artifacts. Limit a run
with a comma-separated study list:

```bash
make run-test STUDIES=wang25
make run STUDIES=aifi,onek1k
```

The `test` profile lowers resource requests and uses isolated 200-cell fixtures
from `test_data/`. RDS-based studies retain an RDS fixture and therefore still
pass through the same `CONVERT_RDS` process as a production run. Create or refresh fixtures
from the production inputs once with `make test-data`; this never modifies
`input_data/`. Override the fixture size with `make test-data TEST_CELLS=500`
when needed. `--skip_qc` stops after the H5AD merge artifacts are created (it
skips only notebook/HTML QC). By default the workflow does not create a
whole-dataset single-cell merge, because all studies must fit in RAM at once.
Enable that separate path on a runner sized for the selected inputs:

```bash
make run STUDIES=all MERGE_SINGLE_CELL=true
```

Nextflow's `work/` directory is the cache for
converted and other intermediate artifacts; only
final deliverables are copied to `output/`.

Build and use both containers with:

```bash
make images
make test-data
make run-test STUDIES=wang25
make run STUDIES=all
```

`make run`, `make run-test`, and `make run-no-qc` build the two content-tagged
local images automatically whenever `NF_PROFILE` contains `docker`; use `make
images` to build them without launching a workflow. `IMAGE_TAG` defaults to
the current image build context and may be set explicitly for a release build.

For a production run whose repository and task filesystem are separate, keep
the repository path readable by Docker and choose an external work/output
location:

```bash
make run NF_PROFILE=docker STUDIES=all \
  WORK_DIR=/path/to/run/work OUTDIR=/path/to/run/output
```

Nextflow stages each declared task input into the Docker work directory. Keep
the repository readable to Nextflow, and ensure `WORK_DIR` supports symbolic
links; do not rely on arbitrary host repository paths being available inside a
container.

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

The content-derived image tag participates in Nextflow task identity, so a
changed build context automatically receives separate resumable cache entries.

## Make shortcuts

The `Makefile` is intentionally a thin convenience layer rather than a second
workflow. Run `make help` to see its targets. Common commands include:

```bash
make lint
make workflow-lint
make test
make verify
make validate-test STUDIES=wang25
make validate-full STUDIES=terekhova23
make run-test STUDIES=wang25
make run-no-qc STUDIES=wang25
make run-test STUDIES=all MERGE_SINGLE_CELL=true
```

For a non-Docker Nextflow run, activate `.venv` so the installed
`pbmc-harmonize`, `pbmc-qc`, `pbmc-merge`, `pbmc-merge-qc`, and `pbmc-manifest`
commands are on `PATH`. The
Docker profile does not require the Python environment on the host.

See `input_data/README.md` for the expected input layout and `PLAN.md` for the
remaining reproducibility work.

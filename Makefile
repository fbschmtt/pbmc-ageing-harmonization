.DEFAULT_GOAL := help

STUDIES ?= all
NF_PROFILE ?= docker,harmony_omp16_openblas1
NXF ?= nextflow
WORK_DIR ?= work
OUTDIR ?= output
IMAGE_TAG ?= $(shell git ls-files -co --exclude-standard docker config src scripts reports pyproject.toml README.md | sort | while IFS= read -r file; do test -f "$$file" && sha256sum "$$file"; done | sha256sum | cut -c1-12)
PYTHON_IMAGE ?= pbmc-ageing-python:$(IMAGE_TAG)
R_IMAGE ?= pbmc-ageing-r:$(IMAGE_TAG)
PIPELINE_REVISION ?= $(shell git describe --always --dirty)
TEST_CELLS ?= 200
TEST_SEED ?= 42
VALIDATE_OUTDIR ?= output/validation
MERGE_SINGLE_CELL ?= false
VENV_PYTHON := .venv/bin/python
VENV_DEPS := .venv/.dev-qc-installed

CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),images,)
MERGE_ARGS = $(if $(filter true 1 yes,$(MERGE_SINGLE_CELL)),--merge_single_cell,)

.PHONY: help install lint workflow-lint test test-unit test-integration verify validate-study validate-test validate-full test-data images image-python image-r run run-no-qc run-test pipeline pipeline-harmonize pipeline-test notebook-sync clean-work

help: ## Show available development commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target> [STUDIES=wang25]\n\n"} /^[a-zA-Z_-]+:.*## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

$(VENV_PYTHON):
	python3 -m venv .venv

$(VENV_DEPS): $(VENV_PYTHON) pyproject.toml
	$(VENV_PYTHON) -m pip install -e '.[dev,qc]'
	touch $(VENV_DEPS)

install: $(VENV_DEPS) ## Install development/QC dependencies into the local .venv

lint: $(VENV_DEPS) ## Run Ruff
	$(VENV_PYTHON) -m ruff check src scripts tests

workflow-lint: ## Lint the Nextflow workflow (requires Nextflow >=25.04)
	$(NXF) lint main.nf

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

test: test-unit ## Backward-compatible alias for unit tests

test-integration: run-test ## Run the deterministic Docker integration workflow

verify: lint test-unit workflow-lint test-integration ## Run all checks and a fresh cached Docker test workflow

validate-study:
	@case "$(STUDIES)" in all|*,*) echo "validation requires exactly one STUDIES value" >&2; exit 2;; esac

validate-test: validate-study image-python ## Validate one test fixture in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --test --validate-only --report-output /result/$(STUDIES).test.json

validate-full: validate-study image-python ## Validate one full input in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work -e PYTHONFAULTHANDLER=1 $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --validate-only --report-output /result/$(STUDIES).json

test-data: images ## Create isolated 200-cell H5AD and RDS smoke-test fixtures
	docker run --rm -v "$(CURDIR):/work" -w /work -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --overwrite
	docker run --rm -v "$(CURDIR):/work" -w /work $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED)

images: image-python image-r ## Build cached local workflow images

image-python:
	docker build -f docker/python.Dockerfile -t $(PYTHON_IMAGE) .

image-r:
	docker build -f docker/r-conversion.Dockerfile -t $(R_IMAGE) .

run: $(CONTAINER_PREREQS) ## Resumable production workflow; opt in with MERGE_SINGLE_CELL=true
	$(NXF) run main.nf -profile $(NF_PROFILE) -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) -resume

run-no-qc: $(CONTAINER_PREREQS) ## Resumable workflow without rendered QC reports
	$(NXF) run main.nf -profile $(NF_PROFILE) -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) --skip_qc -resume

run-test: OUTDIR = output/test
run-test: MERGE_SINGLE_CELL = true
run-test: $(CONTAINER_PREREQS) ## Fresh deterministic 200-cell integration workflow, including single-cell merge
	$(NXF) run main.nf -profile $(NF_PROFILE),test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS)

pipeline: run ## Backward-compatible alias for run

pipeline-harmonize: run-no-qc ## Backward-compatible alias for run-no-qc

pipeline-test: run-test ## Backward-compatible alias for run-test

notebook-sync: $(VENV_DEPS) ## Regenerate the merge QC notebook from its paired Python source
	mkdir -p .cache/jupyter
	JUPYTER_DATA_DIR="$(CURDIR)/.cache/jupyter" $(VENV_PYTHON) -m jupytext --to ipynb --output reports/merge_qc_report.ipynb reports/merge_qc_report.py

clean-work: ## Ask Nextflow to remove obsolete cached work directories
	$(NXF) clean -f

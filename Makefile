.DEFAULT_GOAL := help

STUDIES ?= all
NF_PROFILE ?= docker
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
MERGED_INPUT ?= output/merged/single_cell_merged.h5ad
VENV_PYTHON := .venv/bin/python
VENV_DEPS := .venv/.dev-qc-installed

CORE_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),images,)
PYTHON_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python,)
MERGE_ARGS = $(if $(filter true 1 yes,$(MERGE_SINGLE_CELL)),--merge_single_cell,)

.PHONY: help install lint workflow-lint test test-unit test-integration test-cell-type-integration verify validate-study validate-test validate-full test-data download-inputs images image-python image-r run run-no-qc run-test run-cell-types run-cell-types-test clean-work

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
	$(NXF) lint cell_type_analysis.nf

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

test: test-unit ## Backward-compatible alias for unit tests

test-integration: run-test ## Run the deterministic Docker integration workflow

test-cell-type-integration: run-cell-types-test ## Run the downstream cell-type Docker integration workflow

verify: lint test-unit workflow-lint test-integration test-cell-type-integration ## Run all checks and both fresh Docker workflows

validate-study:
	@case "$(STUDIES)" in all|*,*) echo "validation requires exactly one STUDIES value" >&2; exit 2;; esac

validate-test: validate-study image-python ## Validate one test fixture in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --test --output /result/$(STUDIES).test.cells.csv.gz --report-output /result/$(STUDIES).test.prepare.json
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --test --prepared-obs /result/$(STUDIES).test.cells.csv.gz --validate-only --report-output /result/$(STUDIES).test.json

validate-full: validate-study image-python ## Validate one full input in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --input /work/$(shell jq -r '.studies["$(STUDIES)"].input' config/studies.json) --output /result/$(STUDIES).cells.csv.gz --report-output /result/$(STUDIES).prepare.json
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work -e PYTHONFAULTHANDLER=1 $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --prepared-obs /result/$(STUDIES).cells.csv.gz --validate-only --report-output /result/$(STUDIES).json

test-data: images ## Create isolated 200-cell H5AD and RDS smoke-test fixtures
	docker run --rm -v "$(CURDIR):/work" -w /work -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --overwrite
	docker run --rm -v "$(CURDIR):/work" -w /work $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED)

download-inputs: $(VENV_DEPS) ## Best-effort public input download; never runs implicitly
	$(VENV_PYTHON) scripts/download_inputs.py --project-root "$(CURDIR)" --studies "$(STUDIES)"

images: image-python image-r ## Build cached local workflow images

image-python:
	docker build -f docker/python.Dockerfile -t $(PYTHON_IMAGE) .

image-r:
	docker build -f docker/r-conversion.Dockerfile -t $(R_IMAGE) .

run: $(CORE_CONTAINER_PREREQS) ## Resumable production workflow; opt in with MERGE_SINGLE_CELL=true
	$(NXF) run main.nf -profile $(NF_PROFILE),core -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) -resume

run-no-qc: $(CORE_CONTAINER_PREREQS) ## Resumable workflow without rendered QC reports
	$(NXF) run main.nf -profile $(NF_PROFILE),core -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) --skip_qc -resume

run-test: OUTDIR = output/test
run-test: MERGE_SINGLE_CELL = true
run-test: $(CORE_CONTAINER_PREREQS) ## Fresh deterministic 200-cell integration workflow, including single-cell merge
	$(NXF) run main.nf -profile $(NF_PROFILE),core,core_test,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS)

run-cell-types: $(PYTHON_CONTAINER_PREREQS) ## Split a merged H5AD and render one residual-variation report per AIFI L2 type
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

run-cell-types-test: OUTDIR = output/cell-type-test
run-cell-types-test: MERGED_INPUT = output/cell-type-test/merged/single_cell_merged.h5ad
run-cell-types-test: run-test $(PYTHON_CONTAINER_PREREQS) ## Fresh Docker downstream test using the core test merge
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

clean-work: ## Ask Nextflow to remove obsolete cached work directories
	$(NXF) clean -f

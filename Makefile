.DEFAULT_GOAL := help

STUDIES ?= all
NF_PROFILE ?= docker
NXF ?= nextflow
WORK_DIR ?= work
OUTDIR ?= output
PYTHON_IMAGE ?= pbmc-ageing-python:local
R_IMAGE ?= pbmc-ageing-r:local
TEST_CELLS ?= 200
TEST_SEED ?= 42
VALIDATE_OUTDIR ?= output/validation

CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),docker-images,)

.PHONY: help install lint workflow-lint test validate-study validate-test validate-full test-data docker-images docker-build docker-build-python docker-build-r pipeline pipeline-harmonize pipeline-test clean-work

help: ## Show available development commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target> [STUDIES=wang25]\n\n"} /^[a-zA-Z_-]+:.*## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install the editable Python package and development/QC dependencies
	python -m pip install -e '.[dev,qc]'

lint: ## Run Ruff
	.venv/bin/ruff check src scripts tests

workflow-lint: ## Lint the Nextflow workflow (requires Nextflow >=25.04)
	$(NXF) lint main.nf

test: ## Run unit tests
	.venv/bin/pytest -q

validate-study:
	@case "$(STUDIES)" in all|*,*) echo "validation requires exactly one STUDIES value" >&2; exit 2;; esac

validate-test: validate-study docker-build-python ## Validate one test fixture in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --test --validate-only --report-output /result/$(STUDIES).test.json

validate-full: validate-study docker-build-python ## Validate one full input in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work -e PYTHONFAULTHANDLER=1 $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --validate-only --report-output /result/$(STUDIES).json

test-data: docker-images ## Create isolated 200-cell H5AD and RDS smoke-test fixtures
	docker run --rm -v "$(CURDIR):/work" -w /work -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --overwrite
	docker run --rm -v "$(CURDIR):/work" -w /work $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED)

docker-images: docker-build ## Build both workflow images required by the Docker profile

docker-build: docker-build-python docker-build-r ## Build both workflow images

docker-build-python: ## Build the harmonization and QC image
	docker build -f docker/python.Dockerfile -t pbmc-ageing-python:local .

docker-build-r: ## Build the Seurat RDS conversion image
	docker build -f docker/r-conversion.Dockerfile -t pbmc-ageing-r:local .

pipeline: $(CONTAINER_PREREQS) ## Run selected studies; Docker images build automatically with NF_PROFILE=docker
	$(NXF) run main.nf -profile $(NF_PROFILE) -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) -resume

pipeline-harmonize: $(CONTAINER_PREREQS) ## Run through harmonization but skip QC reports
	$(NXF) run main.nf -profile $(NF_PROFILE) -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --skip_qc -resume

pipeline-test: OUTDIR = output/test
pipeline-test: $(CONTAINER_PREREQS) ## Run the selected 200-cell test workflow
	$(NXF) run main.nf -profile $(NF_PROFILE),test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES)

clean-work: ## Ask Nextflow to remove obsolete cached work directories
	$(NXF) clean -f

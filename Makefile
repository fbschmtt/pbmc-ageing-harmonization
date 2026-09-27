.DEFAULT_GOAL := help

STUDIES ?= all
# Fast development tests omit Nehar-Belaid because its much larger fixture
# dominates runtime. Pass STUDIES=all (or a list) to override this selection.
TEST_STUDIES ?= aida25,aifi,fachrul26,onek1k,perez22,terekhova23,wang25
RUN_TEST_STUDIES = $(if $(filter command line environment environment override,$(origin STUDIES)),$(STUDIES),$(TEST_STUDIES))
NF_PROFILE ?= docker
NXF ?= nextflow
WORK_DIR ?= work
OUTDIR ?= output
# IMAGE_TAG remains an optional common override for release builds. By default,
# each image is tagged from only the files copied by its own Docker build context.
IMAGE_TAG ?=
PYTHON_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell git ls-files -co --exclude-standard docker/python.Dockerfile config src pyproject.toml requirements.lock | sort | while IFS= read -r file; do test -f "$$file" && sha256sum "$$file"; done | sha256sum | cut -c1-12))
R_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell git ls-files -co --exclude-standard docker/r-conversion.Dockerfile scripts/convert_rds.R | sort | while IFS= read -r file; do test -f "$$file" && sha256sum "$$file"; done | sha256sum | cut -c1-12))
PYTHON_IMAGE ?= pbmc-ageing-python:$(PYTHON_IMAGE_TAG)
R_IMAGE ?= pbmc-ageing-r:$(R_IMAGE_TAG)
PIPELINE_REVISION ?= $(shell git describe --always --dirty)
TEST_CELLS ?= 200
TEST_SEED ?= 42
VALIDATE_OUTDIR ?= output/validation
MERGE_SINGLE_CELL ?= false
MERGED_INPUT ?= output/merged/single_cell_merged.h5ad
CELL_TYPE ?= cd14-monocyte
CELL_TYPE_INPUT ?= $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
CELL_TYPE_ANALYSIS_INPUT ?= $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
VENV_PYTHON := .venv/bin/python
VENV_DEPS := .venv/.dev-qc-installed

CORE_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),images,)
PYTHON_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python,)
MERGE_ARGS = $(if $(filter true 1 yes,$(MERGE_SINGLE_CELL)),--merge_single_cell,)

.PHONY: help install lint workflow-lint docs-check test test-unit test-integration test-cell-type-integration verify validate-study validate-test validate-full test-data download-inputs images image-python image-r run run-no-qc run-test split-cell-types split-cell-types-test run-cell-types run-cell-types-test run-cell-types-test-existing run-cell-type run-cell-type-test render-cell-type render-cell-type-test check-cell-type-prerequisites check-cell-type-test-prerequisites check-cell-type-render-prerequisites check-cell-type-test-render-prerequisites clean-work

help: ## Show available development commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target> [STUDIES=wang25]\n\n"} /^[a-zA-Z_-]+:.*## / {printf "  %-20s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

$(VENV_PYTHON):
	python3 -m venv .venv

$(VENV_DEPS): $(VENV_PYTHON) pyproject.toml requirements.lock requirements.dev.lock
	$(VENV_PYTHON) -m pip install -r requirements.dev.lock
	$(VENV_PYTHON) -m pip install --no-deps --no-build-isolation -e '.[dev,qc]'
	touch $(VENV_DEPS)

install: $(VENV_DEPS) ## Install development/QC dependencies into the local .venv

lint: $(VENV_DEPS) ## Run Ruff
	$(VENV_PYTHON) -m ruff check src scripts tests

workflow-lint: ## Lint the Nextflow workflow (requires Nextflow >=25.04)
	$(NXF) lint main.nf
	$(NXF) lint cell_type_analysis.nf

docs-check: ## Validate repository Markdown links and documented Make targets
	python3 scripts/check_docs.py

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

test: test-unit ## Backward-compatible alias for unit tests

test-integration: run-test

test-cell-type-integration: run-cell-types-test-existing

verify: lint test-unit workflow-lint docs-check test-integration test-cell-type-integration ## Run all checks and fresh Docker workflows for all studies
verify: TEST_STUDIES = all

validate-study:
	@case "$(STUDIES)" in all|*,*) echo "validation requires exactly one STUDIES value" >&2; exit 2;; esac

validate-test: validate-study image-python ## Validate one test fixture in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --test --output /result/$(STUDIES).test.cells.csv.gz --source-obs-output /result/$(STUDIES).test.source_obs.csv.gz --report-output /result/$(STUDIES).test.prepare.json
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --test --prepared-obs /result/$(STUDIES).test.cells.csv.gz --validate-only --report-output /result/$(STUDIES).test.json

validate-full: validate-study image-python ## Validate one full input in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --input /work/$(shell jq -r '.studies["$(STUDIES)"].input' config/studies.json) --output /result/$(STUDIES).cells.csv.gz --source-obs-output /result/$(STUDIES).source_obs.csv.gz --report-output /result/$(STUDIES).prepare.json
	docker run --rm -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work -e PYTHONFAULTHANDLER=1 $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --prepared-obs /result/$(STUDIES).cells.csv.gz --validate-only --report-output /result/$(STUDIES).json

test-data: images ## Create isolated 200-cell H5AD and RDS smoke-test fixtures
	docker run --rm -v "$(CURDIR):/work" -w /work -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)" --overwrite
	docker run --rm -v "$(CURDIR):/work" -w /work $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)"

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
run-test: $(CORE_CONTAINER_PREREQS) ## Fresh deterministic test workflow; excludes Nehar-Belaid unless STUDIES overrides it
	$(NXF) run main.nf -profile $(NF_PROFILE),core,core_test,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(RUN_TEST_STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS)

run-cell-types: $(PYTHON_CONTAINER_PREREQS) ## Split a merged H5AD and render one residual-variation report per AIFI L2 type
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

split-cell-types: $(PYTHON_CONTAINER_PREREQS) ## Create reusable raw-count cell-type splits from an existing merged H5AD
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --split_only --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

check-cell-type-prerequisites:
	@if test -f "$(CELL_TYPE_INPUT)"; then \
		echo "Ready: $(CELL_TYPE_INPUT)"; \
	elif test -f "$(MERGED_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make split-cell-types CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Missing merge: $(MERGED_INPUT). Run: make run-cell-types CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	fi

run-cell-type: check-cell-type-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-cell-type-render-prerequisites:
	@if test ! -f "$(CELL_TYPE_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make run-cell-type CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	elif test ! -f "$(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; then \
		echo "Missing analysis: $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad. Run: make run-cell-type CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Ready: $(CELL_TYPE_INPUT) and $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; \
	fi

render-cell-type: check-cell-type-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing analysis without recomputing it
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

run-cell-types-test: OUTDIR = output/test
run-cell-types-test: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
run-cell-types-test: run-test run-cell-types-test-existing ## Fresh core and downstream Docker test workflow

run-cell-types-test-existing: OUTDIR = output/test
run-cell-types-test-existing: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
run-cell-types-test-existing: $(PYTHON_CONTAINER_PREREQS)
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

split-cell-types-test: OUTDIR = output/test
split-cell-types-test: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
split-cell-types-test: $(PYTHON_CONTAINER_PREREQS) ## Create reusable test splits from an existing test merge; does not run the core workflow
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --split_only --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

run-cell-type-test: OUTDIR = output/test
run-cell-type-test: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
run-cell-type-test: CELL_TYPE_INPUT = $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
run-cell-type-test: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
check-cell-type-test-prerequisites: OUTDIR = output/test
check-cell-type-test-prerequisites: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
check-cell-type-test-prerequisites: CELL_TYPE_INPUT = $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
check-cell-type-test-prerequisites:
	@if test -f "$(CELL_TYPE_INPUT)"; then \
		echo "Ready: $(CELL_TYPE_INPUT)"; \
	elif test -f "$(MERGED_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make split-cell-types-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Missing test merge: $(MERGED_INPUT). Run: make run-cell-types-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	fi

run-cell-type-test: check-cell-type-test-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing test split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-cell-type-test-render-prerequisites: OUTDIR = output/test
check-cell-type-test-render-prerequisites: CELL_TYPE_INPUT = $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
check-cell-type-test-render-prerequisites: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
check-cell-type-test-render-prerequisites:
	@if test ! -f "$(CELL_TYPE_INPUT)"; then \
		echo "Missing test split: $(CELL_TYPE_INPUT). Run: make run-cell-type-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	elif test ! -f "$(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; then \
		echo "Missing test analysis: $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad. Run: make run-cell-type-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Ready: $(CELL_TYPE_INPUT) and $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; \
	fi

render-cell-type-test: OUTDIR = output/test
render-cell-type-test: CELL_TYPE_INPUT = $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
render-cell-type-test: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
render-cell-type-test: check-cell-type-test-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing test analysis without recomputing it
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

clean-work: ## Ask Nextflow to remove obsolete cached work directories
	$(NXF) clean -f

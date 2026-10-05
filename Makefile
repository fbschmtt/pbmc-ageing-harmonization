.DEFAULT_GOAL := help

STUDIES ?= all
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
TEST_DATA_OVERWRITE ?= false
VALIDATE_OUTDIR ?= output/validation
MERGE_SINGLE_CELL ?= false
MERGED_INPUT ?= $(OUTDIR)/merged/single_cell_merged.h5ad
INTEGRATION_BENCHMARK_INPUT ?= $(MERGED_INPUT)
CELL_TYPE ?= cd14-monocyte
CELL_TYPE_INPUT ?= $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
CELL_TYPE_ANALYSIS_INPUT ?= $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
PSEUDOBULK_INPUT ?= $(OUTDIR)/merged/pseudobulk_merged.h5ad
DE_RESULTS_DIR ?= $(OUTDIR)/differential_expression
SYNTHETIC_DE_DIR ?= $(OUTDIR)/synthetic_de
VENV_PYTHON := .venv/bin/python
VENV_DEPS := .venv/.dev-qc-installed

DOCKER_WORKSPACE = docker run --rm --user "$(shell id -u):$(shell id -g)" -v "$(CURDIR):/work" -w /work
DOCKER_WORKSPACE_RO = docker run --rm --user "$(shell id -u):$(shell id -g)" -v "$(CURDIR):/work:ro" -v "$(abspath $(VALIDATE_OUTDIR)):/result" -w /work
WANG_RDS_SOURCE := input_data/wang25/scRNA-seqProcessedLabelledObject.rds
WANG_TEST_RDS := test_data/scRNA-seqProcessedLabelledObject.rds
SELECTS_WANG = $(if $(filter all,$(1)),yes,$(if $(findstring wang25,$(1)),yes))
PRODUCTION_R_IMAGE_PREREQ = $(if $(and $(call SELECTS_WANG,$(STUDIES)),$(wildcard $(WANG_RDS_SOURCE))),image-r)
TEST_R_IMAGE_PREREQ = $(if $(and $(call SELECTS_WANG,$(STUDIES)),$(wildcard $(WANG_TEST_RDS))),image-r)
CORE_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python $(PRODUCTION_R_IMAGE_PREREQ),)
TEST_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python $(TEST_R_IMAGE_PREREQ),)
TEST_DATA_CONTAINER_PREREQS := image-python $(PRODUCTION_R_IMAGE_PREREQ)
PYTHON_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python,)
MERGE_ARGS = $(if $(filter true 1 yes,$(MERGE_SINGLE_CELL)),--merge_single_cell,)

.PHONY: help install lint workflow-lint docs-check test-unit verify run-all-test run-trajectory-test \
	validate-study validate-test validate-full test-data download-inputs images image-python \
	download-inputs-strict check-download-inputs \
	image-r run run-no-qc run-test run-all split-cell-types split-cell-types-test \
	run-cell-type-analysis run-cell-type-analysis-test \
	run-cell-type run-cell-type-test render-cell-type render-cell-type-test \
	render-cell-type-reports render-cell-type-reports-test \
	run-integration-benchmark run-integration-benchmark-test check-integration-benchmark-prerequisites check-integration-benchmark-test \
	run-de run-de-test run-de-synthetic-test \
	render-trajectory-report render-trajectory-report-test \
	check-synthetic-de-results check-synthetic-de-reports check-de-prerequisites \
	check-cell-type-prerequisites check-cell-type-analysis-test-prerequisites check-cell-type-test-prerequisites \
	check-cell-type-render-prerequisites check-cell-type-test-render-prerequisites \
	check-cell-type-reports-prerequisites clean-core-output clean-analysis-output clean-work

help: ## Show available development commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target> [STUDIES=<comma-separated-list>]\n\n"} /^[a-zA-Z_-]+:.*## / {printf "  %-31s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

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
	$(NXF) lint differential_expression.nf
	$(NXF) lint integration_benchmark.nf

docs-check: ## Validate repository Markdown links and documented Make targets
	python3 scripts/check_docs.py

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

verify: lint test-unit workflow-lint docs-check run-all-test ## Run all checks and fresh Docker workflows for all studies

# Keep benchmark, DE, and reports ordered even under `make -j`: each stage
# consumes the fresh core output or synthetic-DE manifest from the prior stage.
run-all-test: run-test ## Run the ordered core, benchmark, positive-DE, and cell-type Docker test workflows
	$(MAKE) run-integration-benchmark-test
	$(MAKE) check-integration-benchmark-test
	$(MAKE) run-de-synthetic-test
	$(MAKE) check-synthetic-de-results
	$(MAKE) run-cell-type-analysis-test
	$(MAKE) check-synthetic-de-reports

run-trajectory-test: check-cell-type-analysis-test-prerequisites ## Exercise synthetic trajectory fits and both report levels from an existing test merge
	$(MAKE) run-de-synthetic-test
	$(MAKE) check-synthetic-de-results
	$(MAKE) run-cell-type-analysis-test
	$(MAKE) check-synthetic-de-reports

validate-study:
	@case "$(STUDIES)" in all|*,*) echo "validation requires exactly one STUDIES value" >&2; exit 2;; esac

validate-test: validate-study image-python ## Validate one test fixture in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	$(DOCKER_WORKSPACE_RO) $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --test --output /result/$(STUDIES).test.cells.csv.gz --source-obs-output /result/$(STUDIES).test.source_obs.csv.gz --report-output /result/$(STUDIES).test.prepare.json
	$(DOCKER_WORKSPACE_RO) $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --test --prepared-obs /result/$(STUDIES).test.cells.csv.gz --validate-only --report-output /result/$(STUDIES).test.json

validate-full: validate-study image-python ## Validate one full input in the Python container
	mkdir -p "$(VALIDATE_OUTDIR)"
	$(DOCKER_WORKSPACE_RO) $(PYTHON_IMAGE) pbmc-prepare --project-root /work --config config/pipeline.json --study $(STUDIES) --input /work/$(shell jq -r '.studies["$(STUDIES)"].input' config/studies.json) --output /result/$(STUDIES).cells.csv.gz --source-obs-output /result/$(STUDIES).source_obs.csv.gz --report-output /result/$(STUDIES).prepare.json
	$(DOCKER_WORKSPACE_RO) -e PYTHONFAULTHANDLER=1 $(PYTHON_IMAGE) pbmc-harmonize --project-root /work --config config/pipeline.json --study $(STUDIES) --prepared-obs /result/$(STUDIES).cells.csv.gz --validate-only --report-output /result/$(STUDIES).json

test-data: $(TEST_DATA_CONTAINER_PREREQS) ## Create missing 200-cell H5AD and RDS smoke-test fixtures
	$(DOCKER_WORKSPACE) -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)" $(if $(filter true 1 yes,$(TEST_DATA_OVERWRITE)),--overwrite,)
ifneq ($(PRODUCTION_R_IMAGE_PREREQ),)
	$(DOCKER_WORKSPACE) $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)" $(if $(filter true 1 yes,$(TEST_DATA_OVERWRITE)),--overwrite,)
endif

download-inputs: image-python ## Best-effort public input download; never runs implicitly
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --studies "$(STUDIES)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,)

download-inputs-strict: image-python ## Explicit download preflight; fail unless every selected artifact is present and valid
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --studies "$(STUDIES)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,) --require-all

check-download-inputs: image-python ## Run downloader contract checks in the Python container
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_download_inputs.py

images: image-python image-r ## Build cached local workflow images

image-python:
	docker build -f docker/python.Dockerfile -t $(PYTHON_IMAGE) .

image-r:
	docker build -f docker/r-conversion.Dockerfile -t $(R_IMAGE) .

run: $(CORE_CONTAINER_PREREQS) ## Resumable production workflow; opt in with MERGE_SINGLE_CELL=true
	$(NXF) run main.nf -profile $(NF_PROFILE),core -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) -resume

run-all: ## Run the core workflow, pseudobulk DE, and all cell-type analyses/reports in order
	$(MAKE) run MERGE_SINGLE_CELL=true
	$(MAKE) run-de
	$(MAKE) run-cell-type-analysis

run-no-qc: $(CORE_CONTAINER_PREREQS) ## Resumable workflow without rendered QC reports
	$(NXF) run main.nf -profile $(NF_PROFILE),core -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) --skip_qc -resume

run-test: OUTDIR = output/test
run-test: MERGE_SINGLE_CELL = true
run-test: $(TEST_CONTAINER_PREREQS) ## Fresh deterministic test workflow for all configured studies; skips unavailable fixtures
	$(NXF) run main.nf -profile $(NF_PROFILE),core,core_test,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS)

run-cell-type-analysis: $(PYTHON_CONTAINER_PREREQS) ## Analyse all cell types and render their reports, including DE results when present
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

check-de-prerequisites:
	@if test ! -f "$(PSEUDOBULK_INPUT)"; then \
		echo "Missing pseudobulk merge: $(PSEUDOBULK_INPUT). Run the core workflow first." >&2; exit 2; \
	else \
		echo "Ready: $(PSEUDOBULK_INPUT)"; \
	fi

run-de: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Run per-study and combined covariate PyDESeq2 models from a pseudobulk H5AD
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression -work-dir $(WORK_DIR) --outdir $(OUTDIR) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

run-de-test: OUTDIR = output/test
run-de-test: PSEUDOBULK_INPUT = output/test/merged/pseudobulk_merged.h5ad
run-de-test: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Run DE from an existing test pseudobulk merge
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-integration-benchmark-prerequisites:
	@if test ! -f "$(INTEGRATION_BENCHMARK_INPUT)"; then \
		echo "Missing merged single-cell input: $(INTEGRATION_BENCHMARK_INPUT). Run the core workflow first." >&2; exit 2; \
	else \
		echo "Ready: $(INTEGRATION_BENCHMARK_INPUT)"; \
	fi

run-integration-benchmark: check-integration-benchmark-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Run optional global Harmony and unintegrated-PCA label/UMAP diagnostics
	$(NXF) run integration_benchmark.nf -profile $(NF_PROFILE),integration_benchmark -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(INTEGRATION_BENCHMARK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

run-integration-benchmark-test: OUTDIR = output/test
run-integration-benchmark-test: INTEGRATION_BENCHMARK_INPUT = output/test/merged/single_cell_merged.h5ad
run-integration-benchmark-test: check-integration-benchmark-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fresh global integration benchmark from an existing test merge
	$(NXF) run integration_benchmark.nf -profile $(NF_PROFILE),integration_benchmark,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(INTEGRATION_BENCHMARK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-integration-benchmark-test: $(PYTHON_CONTAINER_PREREQS) ## Check thin integration-benchmark test artifact and method contract
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_integration_benchmark.py --outdir output/test

run-de-synthetic-test: OUTDIR = output/test
run-de-synthetic-test: SYNTHETIC_DE_DIR = $(OUTDIR)/synthetic_de
run-de-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Generate a compact positive DE fixture and fit it with production filters enabled
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad"
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression,synthetic_de_test -work-dir $(WORK_DIR) --outdir "$(OUTDIR)" --pseudobulk_input "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-synthetic-de-results: $(PYTHON_CONTAINER_PREREQS) ## Check completed synthetic DE fits and planted age markers
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_synthetic_de.py results --outdir output/test

check-synthetic-de-reports: $(PYTHON_CONTAINER_PREREQS) ## Check synthetic DE warnings, markers, and plots in standard test reports
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_synthetic_de.py reports --outdir output/test

check-trajectory-report-prerequisites:
	@if test ! -f "$(DE_RESULTS_DIR)/differential_expression.json"; then \
		echo "Missing DE manifest: $(DE_RESULTS_DIR)/differential_expression.json. Run make run-de first." >&2; exit 2; \
	else \
		echo "Ready to render the trajectory report from $(DE_RESULTS_DIR)"; \
	fi

render-trajectory-report: check-trajectory-report-prerequisites image-python ## Rerender the cross-cell-type trajectory report from existing DE results
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) pbmc-trajectory-report --differential-expression-dir "$(DE_RESULTS_DIR)" --output-dir "$(DE_RESULTS_DIR)/trajectory_analysis" --config config/pipeline.json --template reports/trajectory_analysis_report.py --project-root /work

render-trajectory-report-test: OUTDIR = output/test
render-trajectory-report-test: DE_RESULTS_DIR = $(OUTDIR)/differential_expression
render-trajectory-report-test: check-trajectory-report-prerequisites image-python ## Rerender the test trajectory report from existing synthetic DE results
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) pbmc-trajectory-report --differential-expression-dir "$(DE_RESULTS_DIR)" --output-dir "$(DE_RESULTS_DIR)/trajectory_analysis" --config config/pipeline.json --template reports/trajectory_analysis_report.py --project-root /work

split-cell-types: $(PYTHON_CONTAINER_PREREQS) ## Create reusable raw-count cell-type splits from an existing merged H5AD
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --split_only --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

check-cell-type-prerequisites:
	@if test -f "$(CELL_TYPE_INPUT)"; then \
		echo "Ready: $(CELL_TYPE_INPUT)"; \
	elif test -f "$(MERGED_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make split-cell-types CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Missing merge: $(MERGED_INPUT). Run: make run MERGE_SINGLE_CELL=true" >&2; exit 2; \
	fi

run-cell-type: check-cell-type-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-cell-type-render-prerequisites:
	@if test ! -f "$(CELL_TYPE_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make run-cell-type CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	elif test ! -f "$(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; then \
		echo "Missing analysis: $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad. Run: make run-cell-type CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Ready: $(CELL_TYPE_INPUT) and $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; \
	fi

render-cell-type: check-cell-type-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing analysis without recomputing it
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-cell-type-reports-prerequisites:
	@if test -z "$(wildcard $(OUTDIR)/cell_type_splits/*.h5ad)"; then \
		echo "Missing cell-type splits under $(OUTDIR)/cell_type_splits. Run make split-cell-types first." >&2; exit 2; \
	elif test -z "$(wildcard $(OUTDIR)/cell_type_analysis/*/analysis.h5ad)"; then \
		echo "Missing cell-type analysis artifacts under $(OUTDIR)/cell_type_analysis. Run make run-cell-type-analysis first." >&2; exit 2; \
	else \
		echo "Ready to render all existing cell-type analyses under $(OUTDIR)"; \
	fi

render-cell-type-reports: check-cell-type-reports-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Rerender all cell-type reports from existing splits and analysis artifacts
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(OUTDIR)/cell_type_splits/*.h5ad" --cell_type_analysis_input "$(OUTDIR)/cell_type_analysis/*" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

check-cell-type-analysis-test-prerequisites: OUTDIR = output/test
check-cell-type-analysis-test-prerequisites: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
check-cell-type-analysis-test-prerequisites:
	@if test ! -f "$(MERGED_INPUT)"; then \
		echo "Missing test merge: $(MERGED_INPUT). Run make run-test first." >&2; exit 2; \
	else \
		echo "Ready: $(MERGED_INPUT)"; \
	fi

run-cell-type-analysis-test: OUTDIR = output/test
run-cell-type-analysis-test: MERGED_INPUT = output/test/merged/single_cell_merged.h5ad
run-cell-type-analysis-test: check-cell-type-analysis-test-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render all cell types from an existing test merge
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

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
		echo "Missing test merge: $(MERGED_INPUT). Run: make run-test" >&2; exit 2; \
	fi

run-cell-type-test: check-cell-type-test-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing test split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

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
render-cell-type-test: DE_RESULTS_DIR = output/test/differential_expression
render-cell-type-test: CELL_TYPE_INPUT = $(OUTDIR)/cell_type_splits/$(CELL_TYPE).h5ad
render-cell-type-test: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/cell_type_analysis/$(CELL_TYPE)
render-cell-type-test: check-cell-type-test-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing test analysis without recomputing it
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

render-cell-type-reports-test: OUTDIR = output/test
render-cell-type-reports-test: DE_RESULTS_DIR = output/test/differential_expression
render-cell-type-reports-test: check-cell-type-reports-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Rerender all test cell-type reports from existing artifacts
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(OUTDIR)/cell_type_splits/*.h5ad" --cell_type_analysis_input "$(OUTDIR)/cell_type_analysis/*" --differential_expression_manifest "$(DE_RESULTS_DIR)/differential_expression.json" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

clean-work: ## Delete completed Nextflow cache entries (destructive)
	$(NXF) clean -f

clean-core-output: ## Delete all published production outputs, preserving output/test and the Nextflow cache
	rm -rf output/prepared output/harmonized output/pseudobulk output/merged output/reports output/qc output/differential_expression output/cell_type_analysis output/cell_type_splits output/run_manifest.json

clean-analysis-output: ## Delete published production DE and cell-type-analysis outputs, preserving core and output/test
	rm -rf output/differential_expression output/cell_type_analysis output/cell_type_splits

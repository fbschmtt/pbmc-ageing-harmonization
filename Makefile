.DEFAULT_GOAL := help

STUDIES ?= all
# Representative default for test fixtures and non-production workflows. Keep
# production `make run` on all studies; override any test selection with
# `STUDIES=<comma-separated-list>` when needed.
TEST_STUDIES ?= aida25,cima,terekhova26
NF_PROFILE ?= docker
NXF ?= nextflow
WORK_DIR ?= work
OUTDIR ?= output
# IMAGE_TAG remains an optional common override for release builds. By default,
# each image is tagged from only the files copied by its own Docker build context.
IMAGE_TAG ?=
PYTHON_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell git ls-files -co --exclude-standard docker/python.Dockerfile config src pyproject.toml requirements.lock | sort | while IFS= read -r file; do test -f "$$file" && sha256sum "$$file"; done | sha256sum | cut -c1-12))
R_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell git ls-files -co --exclude-standard docker/r-conversion.Dockerfile scripts/convert_rds.R | sort | while IFS= read -r file; do test -f "$$file" && sha256sum "$$file"; done | sha256sum | cut -c1-12))
# Freeze the content hashes at parse time. Recursive expansion can otherwise
# recompute a different tag between an image prerequisite and its recipe.
PYTHON_IMAGE_TAG := $(PYTHON_IMAGE_TAG)
R_IMAGE_TAG := $(R_IMAGE_TAG)
PYTHON_IMAGE ?= pbmc-ageing-python:$(PYTHON_IMAGE_TAG)
R_IMAGE ?= pbmc-ageing-r:$(R_IMAGE_TAG)
PIPELINE_REVISION ?= $(shell git describe --always --dirty)
TEST_CELLS ?= 200
TEST_SEED ?= 42
TEST_DATA_OVERWRITE ?= false
VALIDATE_OUTDIR ?= output/validation
MERGE_SINGLE_CELL ?= true
MERGED_INPUT ?= $(OUTDIR)/harmonization/merged/single_cell_merged.h5ad
INTEGRATION_BENCHMARK_INPUT ?= $(MERGED_INPUT)
INTEGRATION_BENCHMARK_LOG_LEVEL ?= INFO
CELL_TYPE ?= cd14-monocyte
GENE ?=
TOP_AGE_GENE_COUNT ?= 10
CELL_TYPE_INPUT ?= $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/split.h5ad
CELL_TYPE_ANALYSIS_INPUT ?= $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/analysis
PSEUDOBULK_INPUT ?= $(OUTDIR)/harmonization/merged/pseudobulk_merged.h5ad
GENE_PROFILE_DIR ?= $(OUTDIR)/analysis/per_gene
GENE_PROFILE_TEST_DIR ?= $(OUTDIR)/analysis/combined/gene_profile_synthetic
GENE_PROFILE_TEST_INPUT ?= $(GENE_PROFILE_TEST_DIR)/pseudobulk_merged.h5ad
PER_TYPE_DIR ?= $(OUTDIR)/analysis/per_type
SYNTHETIC_DE_DIR ?= $(OUTDIR)/analysis/combined/synthetic_de
AGE_CLOCK_DIR ?= $(PER_TYPE_DIR)
AGE_CLOCK_EXTERNAL_STUDIES ?= nehar_belaid26,perez22
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
TEST_CONTAINER_PREREQS = $(if $(findstring docker,$(NF_PROFILE)),image-python,)
TEST_DATA_CONTAINER_PREREQS = image-python
PYTHON_CONTAINER_PREREQS := $(if $(findstring docker,$(NF_PROFILE)),image-python,)
MERGE_ARGS = $(if $(filter true 1 yes,$(MERGE_SINGLE_CELL)),--merge_single_cell,)

.PHONY: help install lint workflow-lint docs-check docs-update test-unit verify run-all-test \
validate-study validate-test validate-full test-data download-inputs download-go-inputs images image-python \
	download-inputs-strict check-download-inputs \
	image-r run run-test run-all split-cell-types split-cell-types-test \
	run-cell-type-analysis run-cell-type-analysis-test \
	run-cell-type run-cell-type-test render-cell-type render-cell-type-test \
	render-cell-type-reports render-cell-type-reports-test \
	render-de-reports render-de-reports-test \
	profile-gene profile-top-age-genes check-top-age-profile-prerequisites \
	profile-top-age-genes-test \
	run-gene-profile-synthetic-test \
	run-age-clock run-age-clocks run-age-clock-synthetic-test check-age-clock-synthetic-test check-age-clock-prerequisites \
	run-integration-benchmark run-integration-benchmark-test check-integration-benchmark-prerequisites check-integration-benchmark-test \
	run-de run-de-test run-de-synthetic-test \
	render-report-index \
	check-synthetic-de-results check-synthetic-de-reports check-de-prerequisites \
	check-cell-type-prerequisites check-cell-type-analysis-test-prerequisites check-cell-type-test-prerequisites \
	check-cell-type-render-prerequisites check-cell-type-test-render-prerequisites \
	check-cell-type-reports-prerequisites

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

docs-update: ## Update generated documentation from config/pipeline.json
	python3 scripts/check_docs.py --write

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

verify: STUDIES = $(TEST_STUDIES)
verify: lint test-unit workflow-lint docs-check run-all-test ## Run all checks and fresh Docker workflows for representative studies

# Keep benchmark, DE, age-clock, and reports ordered even under `make -j`: each stage
# consumes the fresh core output or synthetic-DE manifest from the prior stage.
run-all-test: STUDIES = $(TEST_STUDIES)
run-all-test: run-test ## Run core, benchmark, DE, age clocks, profiles, and cell-type reports
	$(MAKE) run-integration-benchmark-test STUDIES="$(STUDIES)"
	$(MAKE) check-integration-benchmark-test STUDIES="$(STUDIES)"
	$(MAKE) run-de-synthetic-test STUDIES="$(STUDIES)"
	$(MAKE) check-synthetic-de-results STUDIES="$(STUDIES)"
	$(MAKE) run-age-clock-synthetic-test STUDIES="$(STUDIES)"
	$(MAKE) check-age-clock-synthetic-test STUDIES="$(STUDIES)"
	$(MAKE) profile-top-age-genes-test STUDIES="$(STUDIES)"
	$(MAKE) run-gene-profile-synthetic-test STUDIES="$(STUDIES)"
	$(MAKE) run-cell-type-analysis-test STUDIES="$(STUDIES)"
	$(MAKE) check-synthetic-de-reports STUDIES="$(STUDIES)"

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

BUILD_TEST_R_IMAGE = $(if $(TEST_R_IMAGE_PREREQ),$(MAKE) image-r,@:)
BUILD_TEST_DATA_R_IMAGE = $(if $(PRODUCTION_R_IMAGE_PREREQ),$(MAKE) image-r,@:)
CREATE_TEST_RDS = $(if $(and $(call SELECTS_WANG,$(STUDIES)),$(wildcard $(WANG_RDS_SOURCE))),$(DOCKER_WORKSPACE) $(R_IMAGE) Rscript scripts/create_test_rds.R --project-root /work --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)" $(if $(filter true 1 yes,$(TEST_DATA_OVERWRITE)),--overwrite,),@:)

test-data: STUDIES = $(TEST_STUDIES)
test-data: $(TEST_DATA_CONTAINER_PREREQS) ## Create missing 200-cell H5AD and RDS smoke-test fixtures for representative studies
	$(DOCKER_WORKSPACE) -e PYTHONPATH=/work/src $(PYTHON_IMAGE) python scripts/create_test_data.py --cells $(TEST_CELLS) --seed $(TEST_SEED) --studies "$(STUDIES)" $(if $(filter true 1 yes,$(TEST_DATA_OVERWRITE)),--overwrite,)
	$(BUILD_TEST_DATA_R_IMAGE)
	$(CREATE_TEST_RDS)

download-inputs: image-python ## Best-effort public input download; never runs implicitly
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --studies "$(STUDIES)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,)

download-go-inputs: image-python ## Download the GO enrichment inputs only
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --artifacts "go_biological_process_gene_sets,go_ontology" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,)

download-inputs-strict: image-python ## Explicit download preflight; fail unless every selected artifact is present and valid
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --studies "$(STUDIES)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,) --require-all

check-download-inputs: image-python ## Run downloader contract checks in the Python container
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_download_inputs.py

images: image-python image-r ## Build cached local workflow images

image-python:
	docker build -f docker/python.Dockerfile -t $(PYTHON_IMAGE) .

image-r:
	docker build -f docker/r-conversion.Dockerfile -t $(R_IMAGE) .

run: $(CORE_CONTAINER_PREREQS) ## Resumable production workflow; merge single-cell data by default (opt out with MERGE_SINGLE_CELL=false)
	$(NXF) run main.nf -profile $(NF_PROFILE),core -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS) -resume
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-all: ## Run the core workflow, pseudobulk DE, and all cell-type analyses/reports in order
	$(MAKE) run MERGE_SINGLE_CELL=true
	$(MAKE) run-de
	$(MAKE) profile-top-age-genes
	$(MAKE) run-cell-type-analysis

run-test: OUTDIR = output/test
run-test: MERGE_SINGLE_CELL = true
run-test: STUDIES = $(TEST_STUDIES)
run-test: $(TEST_CONTAINER_PREREQS) ## Fresh deterministic test workflow for representative studies; skips unavailable fixtures
	$(BUILD_TEST_R_IMAGE)
	$(NXF) run main.nf -profile $(NF_PROFILE),core,core_test,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --studies $(STUDIES) --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION) $(MERGE_ARGS)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-cell-type-analysis: $(PYTHON_CONTAINER_PREREQS) ## Analyse cell types and render single cell analysis reports
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-de-prerequisites:
	@if test ! -f "$(PSEUDOBULK_INPUT)"; then \
		echo "Missing pseudobulk merge: $(PSEUDOBULK_INPUT). Run the core workflow first." >&2; exit 2; \
	else \
		echo "Ready: $(PSEUDOBULK_INPUT)"; \
	fi

run-de: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit pseudobulk DE models and render per-type DE reports
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression -work-dir $(WORK_DIR) --outdir $(OUTDIR) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-de-test: OUTDIR = output/test
run-de-test: PSEUDOBULK_INPUT = output/test/harmonization/merged/pseudobulk_merged.h5ad
run-de-test: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit DE and render per-type reports from an existing test pseudobulk merge
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

render-de-reports: image-python ## Render per-type DE reports from existing model results
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression -work-dir $(WORK_DIR) --outdir $(OUTDIR) --per_type_de_dir "$(PER_TYPE_DIR)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

render-de-reports-test: OUTDIR = output/test
render-de-reports-test: PER_TYPE_DIR = $(OUTDIR)/analysis/per_type
render-de-reports-test: render-de-reports ## Render per-type DE reports from existing synthetic-DE results

render-report-index: ## Rebuild the central HTML index for the reports currently published under OUTDIR
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

profile-gene: image-python check-de-prerequisites ## Refit merged DESeq2 for one cell type and profile one gene; set CELL_TYPE and GENE
	@test -n "$(CELL_TYPE)" || { echo "Set CELL_TYPE to a cell-type label or slug" >&2; exit 2; }
	@test -n "$(GENE)" || { echo "Set GENE to a gene identifier" >&2; exit 2; }
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/profile_gene.py --input "$(PSEUDOBULK_INPUT)" --config config/pipeline.json --cell-type "$(CELL_TYPE)" --gene "$(GENE)" --output-dir "$(GENE_PROFILE_DIR)" --differential-expression-dir "$(PER_TYPE_DIR)"

check-top-age-profile-prerequisites: check-de-prerequisites
	@if ! find "$(PER_TYPE_DIR)" -mindepth 3 -maxdepth 3 -name merged.csv -print -quit 2>/dev/null | grep -q .; then \
		echo "Missing per-type merged age-DE result files under $(PER_TYPE_DIR). Run make run-de first." >&2; exit 2; \
	else \
		echo "Ready: merged age-DE results under $(PER_TYPE_DIR)"; \
	fi

profile-top-age-genes: image-python check-top-age-profile-prerequisites ## Profile top age-DE genes from the merged model in each cell type
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/profile_top_age_genes.py --input "$(PSEUDOBULK_INPUT)" --config config/pipeline.json --differential-expression-dir "$(PER_TYPE_DIR)" --output-dir "$(GENE_PROFILE_DIR)" --top-n "$(TOP_AGE_GENE_COUNT)"

profile-top-age-genes-test: OUTDIR = output/test
profile-top-age-genes-test: PSEUDOBULK_INPUT = $(OUTDIR)/analysis/combined/synthetic_de/pseudobulk_merged.h5ad
profile-top-age-genes-test: PER_TYPE_DIR = $(OUTDIR)/analysis/per_type
profile-top-age-genes-test: GENE_PROFILE_DIR = $(OUTDIR)/analysis/per_gene
profile-top-age-genes-test: TOP_AGE_GENE_COUNT = 2
profile-top-age-genes-test: image-python check-top-age-profile-prerequisites ## Profile top synthetic merged age-DE genes per cell type
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/profile_top_age_genes.py --input "$(PSEUDOBULK_INPUT)" --config config/pipeline.json --differential-expression-dir "$(PER_TYPE_DIR)" --output-dir "$(GENE_PROFILE_DIR)" --top-n "$(TOP_AGE_GENE_COUNT)"

run-gene-profile-synthetic-test: OUTDIR = output/test
run-gene-profile-synthetic-test: GENE_PROFILE_TEST_DIR = $(OUTDIR)/analysis/combined/gene_profile_synthetic
run-gene-profile-synthetic-test: GENE_PROFILE_TEST_INPUT = $(GENE_PROFILE_TEST_DIR)/pseudobulk_merged.h5ad
run-gene-profile-synthetic-test: GENE_PROFILE_DIR = $(OUTDIR)/analysis/per_gene
run-gene-profile-synthetic-test: image-python ## Fit and render a focused gene profile from a synthetic pseudobulk fixture
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(GENE_PROFILE_TEST_INPUT)"
	$(MAKE) profile-gene OUTDIR="$(OUTDIR)" PSEUDOBULK_INPUT="$(GENE_PROFILE_TEST_INPUT)" GENE_PROFILE_DIR="$(GENE_PROFILE_DIR)" CELL_TYPE="CD14 monocyte" GENE=SYNTH_AGE_MARKER_1
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_gene_profile.py --profile-dir "$(GENE_PROFILE_DIR)/cd14-monocyte/synth-age-marker-1" --gene SYNTH_AGE_MARKER_1

check-integration-benchmark-prerequisites:
	@if test ! -f "$(INTEGRATION_BENCHMARK_INPUT)"; then \
		echo "Missing merged single-cell input: $(INTEGRATION_BENCHMARK_INPUT). Run the core workflow first." >&2; exit 2; \
	else \
		echo "Ready: $(INTEGRATION_BENCHMARK_INPUT)"; \
	fi

run-integration-benchmark: check-integration-benchmark-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Run optional global embedding and type-calling benchmark
	$(NXF) run integration_benchmark.nf -profile $(NF_PROFILE),integration_benchmark -ansi-log false -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(INTEGRATION_BENCHMARK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) --log_level $(INTEGRATION_BENCHMARK_LOG_LEVEL) -resume
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-integration-benchmark-test: OUTDIR = output/test
run-integration-benchmark-test: INTEGRATION_BENCHMARK_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
run-integration-benchmark-test: check-integration-benchmark-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fresh global integration benchmark from an existing test merge
	$(NXF) run integration_benchmark.nf -profile $(NF_PROFILE),integration_benchmark,test -ansi-log false -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(INTEGRATION_BENCHMARK_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) --log_level $(INTEGRATION_BENCHMARK_LOG_LEVEL)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-integration-benchmark-test: $(PYTHON_CONTAINER_PREREQS) ## Check thin integration-benchmark test artifact and method contract
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_integration_benchmark.py --outdir output/test

run-de-synthetic-test: OUTDIR = output/test
run-de-synthetic-test: SYNTHETIC_DE_DIR = $(OUTDIR)/analysis/combined/synthetic_de
run-de-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Generate a positive DE fixture, fit it, and render per-type reports
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad"
	$(NXF) run differential_expression.nf -profile $(NF_PROFILE),differential_expression,synthetic_de_test -work-dir $(WORK_DIR) --outdir "$(OUTDIR)" --pseudobulk_input "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-age-clock-synthetic-test: OUTDIR = output/test
run-age-clock-synthetic-test: SYNTHETIC_DE_DIR = $(OUTDIR)/analysis/combined/synthetic_de
run-age-clock-synthetic-test: AGE_CLOCK_DIR = $(OUTDIR)/analysis/per_type
run-age-clock-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Smoke-test the default CELL_TYPE age-clock fit on synthetic DE pseudobulks
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad"
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock --input "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad" --config config/pipeline.json --output-dir "$(AGE_CLOCK_DIR)" --external-studies synthetic_study_c --cell-type "$(CELL_TYPE)"
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock_report --report-root "$(AGE_CLOCK_DIR)" --cell-type "$(CELL_TYPE)" --template reports/age_clock_report.py --comparison-template reports/age_clock_comparison_report.py --project-root /work
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-age-clock-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Check that expected synthetic age-clock output files exist
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_age_clock.py --outdir output/test --cell-type "$(CELL_TYPE)"

check-age-clock-prerequisites:
	@if test ! -f "$(PSEUDOBULK_INPUT)"; then \
		echo "Missing pseudobulk input: $(PSEUDOBULK_INPUT). Run the core workflow first." >&2; exit 2; \
	elif test -z "$(AGE_CLOCK_EXTERNAL_STUDIES)"; then \
		echo "Set AGE_CLOCK_EXTERNAL_STUDIES to one or more comma-separated study IDs." >&2; exit 2; \
	else \
		echo "Ready to train on all pseudobulk studies except $(AGE_CLOCK_EXTERNAL_STUDIES)"; \
	fi

run-age-clock: check-age-clock-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit and externally evaluate one cell-type age clock; set CELL_TYPE and AGE_CLOCK_EXTERNAL_STUDIES
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock --input "$(PSEUDOBULK_INPUT)" --config config/pipeline.json --output-dir "$(AGE_CLOCK_DIR)" --external-studies "$(AGE_CLOCK_EXTERNAL_STUDIES)" --cell-type "$(CELL_TYPE)"
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock_report --report-root "$(AGE_CLOCK_DIR)" --cell-type "$(CELL_TYPE)" --template reports/age_clock_report.py --comparison-template reports/age_clock_comparison_report.py --project-root /work
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

run-age-clocks: check-age-clock-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit and externally evaluate age clocks for every pseudobulk cell type
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock --input "$(PSEUDOBULK_INPUT)" --config config/pipeline.json --output-dir "$(AGE_CLOCK_DIR)" --external-studies "$(AGE_CLOCK_EXTERNAL_STUDIES)" --all-cell-types
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.age_clock_report --report-root "$(AGE_CLOCK_DIR)" --template reports/age_clock_report.py --comparison-template reports/age_clock_comparison_report.py --project-root /work
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-synthetic-de-results: $(PYTHON_CONTAINER_PREREQS) ## Check completed synthetic DE fits and planted age markers
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_synthetic_de.py results --outdir output/test

check-synthetic-de-reports: $(PYTHON_CONTAINER_PREREQS) ## Check synthetic DE and residual reports
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_synthetic_de.py reports --outdir output/test

split-cell-types: $(PYTHON_CONTAINER_PREREQS) ## Create reusable raw-count cell-type splits from an existing merged H5AD
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --split_only --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION) -resume

check-cell-type-prerequisites:
	@if test -f "$(CELL_TYPE_INPUT)"; then \
		echo "Ready: $(CELL_TYPE_INPUT)"; \
	elif test -f "$(MERGED_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make split-cell-types CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Missing merge: $(MERGED_INPUT). Run: make run (single-cell merge is enabled by default)" >&2; exit 2; \
	fi

run-cell-type: check-cell-type-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

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
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-cell-type-reports-prerequisites:
	@if test -z "$(wildcard $(OUTDIR)/analysis/per_type/*/split.h5ad)"; then \
		echo "Missing cell-type splits under $(OUTDIR)/analysis/per_type. Run make split-cell-types first." >&2; exit 2; \
	elif test -z "$(wildcard $(OUTDIR)/analysis/per_type/*/analysis/analysis.h5ad)"; then \
		echo "Missing cell-type analysis artifacts under $(OUTDIR)/analysis/per_type. Run make run-cell-type-analysis first." >&2; exit 2; \
	else \
		echo "Ready to render all existing cell-type analyses under $(OUTDIR)"; \
	fi

render-cell-type-reports: check-cell-type-reports-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Rerender all cell-type reports from existing splits and analysis artifacts
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(OUTDIR)/analysis/per_type/*/split.h5ad" --cell_type_analysis_input "$(OUTDIR)/analysis/per_type/*/analysis" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-cell-type-analysis-test-prerequisites: OUTDIR = output/test
check-cell-type-analysis-test-prerequisites: MERGED_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
check-cell-type-analysis-test-prerequisites:
	@if test ! -f "$(MERGED_INPUT)"; then \
		echo "Missing test merge: $(MERGED_INPUT). Run make run-test first." >&2; exit 2; \
	else \
		echo "Ready: $(MERGED_INPUT)"; \
	fi

run-cell-type-analysis-test: OUTDIR = output/test
run-cell-type-analysis-test: MERGED_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
run-cell-type-analysis-test: check-cell-type-analysis-test-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render all cell types from an existing test merge
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

split-cell-types-test: OUTDIR = output/test
split-cell-types-test: MERGED_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
split-cell-types-test: $(PYTHON_CONTAINER_PREREQS) ## Create reusable test splits from an existing test merge; does not run the core workflow
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --merged_input "$(MERGED_INPUT)" --split_only --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)

run-cell-type-test: OUTDIR = output/test
run-cell-type-test: MERGED_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
run-cell-type-test: CELL_TYPE_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/split.h5ad
run-cell-type-test: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/analysis
check-cell-type-test-prerequisites: OUTDIR = output/test
check-cell-type-test-prerequisites: MERGED_INPUT = output/test/harmonization/merged/single_cell_merged.h5ad
check-cell-type-test-prerequisites: CELL_TYPE_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/split.h5ad
check-cell-type-test-prerequisites:
	@if test -f "$(CELL_TYPE_INPUT)"; then \
		echo "Ready: $(CELL_TYPE_INPUT)"; \
	elif test -f "$(MERGED_INPUT)"; then \
		echo "Missing split: $(CELL_TYPE_INPUT). Run: make split-cell-types-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Missing test merge: $(MERGED_INPUT). Run: make run-test" >&2; exit 2; \
	fi

run-cell-type-test: check-cell-type-test-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render one existing test split; default CELL_TYPE=cd14-monocyte
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

check-cell-type-test-render-prerequisites: OUTDIR = output/test
check-cell-type-test-render-prerequisites: CELL_TYPE_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/split.h5ad
check-cell-type-test-render-prerequisites: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/analysis
check-cell-type-test-render-prerequisites:
	@if test ! -f "$(CELL_TYPE_INPUT)"; then \
		echo "Missing test split: $(CELL_TYPE_INPUT). Run: make run-cell-type-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	elif test ! -f "$(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; then \
		echo "Missing test analysis: $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad. Run: make run-cell-type-test CELL_TYPE=$(CELL_TYPE)" >&2; exit 2; \
	else \
		echo "Ready: $(CELL_TYPE_INPUT) and $(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad"; \
	fi

render-cell-type-test: OUTDIR = output/test
render-cell-type-test: CELL_TYPE_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/split.h5ad
render-cell-type-test: CELL_TYPE_ANALYSIS_INPUT = $(OUTDIR)/analysis/per_type/$(CELL_TYPE)/analysis
render-cell-type-test: check-cell-type-test-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing test analysis without recomputing it
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

render-cell-type-reports-test: OUTDIR = output/test
render-cell-type-reports-test: check-cell-type-reports-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Rerender all test cell-type reports from existing artifacts
	$(NXF) run cell_type_analysis.nf -profile $(NF_PROFILE),cell_type_analysis,test -work-dir $(WORK_DIR) --outdir $(OUTDIR) --cell_type_input "$(OUTDIR)/analysis/per_type/*/split.h5ad" --cell_type_analysis_input "$(OUTDIR)/analysis/per_type/*/analysis" --python_image $(PYTHON_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

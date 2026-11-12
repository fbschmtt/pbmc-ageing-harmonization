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
RUN_CLOCKS ?= false
# Empty uses the mode default: benchmark in full, skip in analysis.
RUN_BENCHMARK ?=
# IMAGE_TAG remains an optional common override for release builds. By default,
# each image is tagged from only the files copied by its own Docker build context.
IMAGE_TAG ?=
PYTHON_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell python3 scripts/image_context_hash.py python))
R_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell python3 scripts/image_context_hash.py r))
SITE_IMAGE_TAG ?= $(if $(IMAGE_TAG),$(IMAGE_TAG),$(shell python3 scripts/image_context_hash.py site))
SITE_BROWSER_IMAGE_TAG ?= $(SITE_IMAGE_TAG)
# Freeze the content hashes at parse time. Recursive expansion can otherwise
# recompute a different tag between an image prerequisite and its recipe.
PYTHON_IMAGE_TAG := $(PYTHON_IMAGE_TAG)
R_IMAGE_TAG := $(R_IMAGE_TAG)
SITE_IMAGE_TAG := $(SITE_IMAGE_TAG)
SITE_BROWSER_IMAGE_TAG := $(SITE_BROWSER_IMAGE_TAG)
PYTHON_IMAGE ?= pbmc-ageing-python:$(PYTHON_IMAGE_TAG)
R_IMAGE ?= pbmc-ageing-r:$(R_IMAGE_TAG)
SITE_IMAGE ?= pbmc-ageing-site:$(SITE_IMAGE_TAG)
SITE_BROWSER_IMAGE ?= pbmc-ageing-site-browser:$(SITE_BROWSER_IMAGE_TAG)
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
ALLOW_PARTIAL_RESULTS ?= false
GENE_PROFILE_DIR ?= $(OUTDIR)/analysis/per_gene
GENE_PROFILE_TEST_DIR ?= $(OUTDIR)/analysis/combined/gene_profile_synthetic
GENE_PROFILE_TEST_INPUT ?= $(GENE_PROFILE_TEST_DIR)/pseudobulk_merged.h5ad
PER_TYPE_DIR ?= $(OUTDIR)/analysis/per_type
SYNTHETIC_DE_DIR ?= $(OUTDIR)/analysis/combined/synthetic_de
AGE_CLOCK_DIR ?= $(PER_TYPE_DIR)
AGE_CLOCK_EXTERNAL_STUDIES ?= nehar_belaid26,perez22
AGE_CLOCK_N_JOBS ?= 20
AGE_CLOCK_METHODS ?=
AGE_CLOCK_METHOD_ARGS = $(if $(strip $(AGE_CLOCK_METHODS)),--methods "$(AGE_CLOCK_METHODS)",)
VERIFY_AGE_CLOCK_METHODS ?= elastic_net
VERIFY_AGE_CLOCK_ALPHAS ?= 0.1
VERIFY_AGE_CLOCK_L1_RATIOS ?= 0.5
VERIFY_CELL_TYPES ?= cd14-monocyte,naive-cd4-t-cell,cd56dim-nk-cell
CELL_TYPE_FILTER_ARGS = $(if $(strip $(CELL_TYPES)),--cell_types "$(CELL_TYPES)",)
DE_CELL_TYPE_ARGS = $(if $(strip $(CELL_TYPE)),--de_cell_type "$(CELL_TYPE)",)
AGE_CLOCK_SEARCH_ARGS = $(if $(strip $(AGE_CLOCK_ALPHAS)),--alphas "$(AGE_CLOCK_ALPHAS)",) $(if $(strip $(AGE_CLOCK_L1_RATIOS)),--l1-ratios "$(AGE_CLOCK_L1_RATIOS)",)
AGE_CLOCK_CALIBRATION_SAMPLES ?= 10
AGE_CLOCK_CALIBRATION_SEED ?= 42
AGE_CLOCK_CALIBRATION_ARGS = --calibration-samples "$(AGE_CLOCK_CALIBRATION_SAMPLES)" --calibration-seed "$(AGE_CLOCK_CALIBRATION_SEED)"
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

.PHONY: run-analysis run-analysis-test check-cell-type-analysis-prerequisites

.PHONY: help install lint workflow-lint docs-check docs-update test-unit verify run-all-test \
image-site image-site-browser \
validate-study validate-test validate-full test-data download-inputs download-gene-set-inputs images image-python \
	download-inputs-strict check-download-inputs \
	image-r run run-test run-all split-cell-types split-cell-types-test \
	run-cell-type-analysis run-cell-type-analysis-test \
	run-cell-type run-cell-type-test render-cell-type render-cell-type-test \
	render-cell-type-reports render-cell-type-reports-test \
	render-de-reports render-de-reports-test \
	profile-gene profile-top-age-genes check-top-age-profile-prerequisites \
	profile-top-age-genes-test \
	run-gene-profile-synthetic-test \
	run-age-clock run-age-clocks aggregate-age-clocks run-age-clock-synthetic-test check-age-clock-synthetic-test check-age-clock-prerequisites \
	run-integration-benchmark run-integration-benchmark-test check-integration-benchmark-prerequisites check-integration-benchmark-test \
	run-de run-de-test run-de-synthetic-test \
	render-report-index \
	check-synthetic-de-results check-synthetic-de-reports check-de-prerequisites \
	check-cell-type-prerequisites check-cell-type-analysis-test-prerequisites check-cell-type-test-prerequisites \
	check-cell-type-render-prerequisites check-cell-type-test-render-prerequisites \
	check-cell-type-reports-prerequisites

help: ## Show available development commands
	@awk 'BEGIN {FS = ":.*## "; printf "Usage: make <target> [STUDIES=<comma-separated-list>]\n\n"} /^[a-zA-Z_-]+( [a-zA-Z_-]+)*:.*## / {printf "  %-31s %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@printf '\nFrontend server: make serve-frontend SITE_PORT=8000 SITE_BIND=127.0.0.1\n  Set SITE_BIND=0.0.0.0 to listen on all host interfaces.\n'

include make/checks.mk
include make/workflows.mk
include make/site.mk

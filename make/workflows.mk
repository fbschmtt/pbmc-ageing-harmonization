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

download-gene-set-inputs: image-python ## Download GO enrichment inputs and static atlas example pathways
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --artifacts "go_biological_process_gene_sets,go_ontology" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,)
	python3 scripts/download_site_pathways.py --output-dir "$(SITE_PATHWAY_DIR)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,)

download-inputs-strict: image-python ## Explicit download preflight; fail unless every selected artifact is present and valid
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/download_inputs.py --project-root /work --studies "$(STUDIES)" $(if $(filter 1 true yes,$(DOWNLOAD_FORCE)),--force,) --require-all

check-download-inputs: image-python ## Run downloader contract checks in the Python container
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_download_inputs.py

images: image-python image-r image-site ## Build cached local workflow images

image-python:
	docker build -f docker/python.Dockerfile -t $(PYTHON_IMAGE) .

image-site: ## Build the locked frontend runtime image (Node, npm, and frontend packages)
	docker build -f docker/site.Dockerfile --target runtime -t $(SITE_IMAGE) .

image-site-browser: ## Build the frontend image with the matching Chromium and system libraries
	docker build -f docker/site.Dockerfile --target browser -t $(SITE_BROWSER_IMAGE) .

image-r:
	docker build -f docker/r-conversion.Dockerfile -t $(R_IMAGE) .

# Test aliases share their production recipe while selecting fresh runs and test paths.
TEST_WORKFLOW_TARGETS := run-test run-de-test run-cell-type-analysis-test split-cell-types-test \
    run-cell-type-test render-cell-type-test render-cell-type-reports-test render-de-reports-test \
    run-integration-benchmark-test run-all-test run-analysis-test
$(TEST_WORKFLOW_TARGETS): OUTDIR = output/test
$(TEST_WORKFLOW_TARGETS): NF_TEST_SUFFIX = ,test
$(TEST_WORKFLOW_TARGETS): NF_RESUME =
NF_TEST_SUFFIX ?=
NF_RESUME ?= -resume
RUN_NEXTFLOW = python3 -m scripts.run_workflow --outdir "$(OUTDIR)" --name "$@" -- $(NXF)
NF_COMMON_ARGS = -profile $(NF_PROFILE),$(WORKFLOW_STAGE)$(NF_TEST_SUFFIX) -work-dir "$(WORK_DIR)" --outdir "$(OUTDIR)" --python_image $(PYTHON_IMAGE) --r_image $(R_IMAGE) --pipeline_revision $(PIPELINE_REVISION)
CLOCK_WORKFLOW_ARGS = --clock_external_studies "$(AGE_CLOCK_EXTERNAL_STUDIES)" $(if $(strip $(AGE_CLOCK_METHODS)),--clock_methods "$(AGE_CLOCK_METHODS)",) --clock_alphas "$(or $(AGE_CLOCK_ALPHAS),0.001,0.01,0.1,1,10)" --clock_l1_ratios "$(or $(AGE_CLOCK_L1_RATIOS),0.1,0.5,0.9,1)" --clock_n_jobs "$(AGE_CLOCK_N_JOBS)" --clock_calibration_samples "$(AGE_CLOCK_CALIBRATION_SAMPLES)" --clock_calibration_seed "$(AGE_CLOCK_CALIBRATION_SEED)"

run run-test: WORKFLOW_STAGE = core
run: $(CORE_CONTAINER_PREREQS) ## Resumable production harmonization and merge workflow
run-test: STUDIES = $(TEST_STUDIES)
run-test: MERGE_SINGLE_CELL = true
run-test: $(TEST_CONTAINER_PREREQS) ## Fresh core test; skips unavailable fixtures
run run-test:
	$(if $(NF_TEST_SUFFIX),$(BUILD_TEST_R_IMAGE),@:)
	$(RUN_NEXTFLOW) run main.nf $(NF_COMMON_ARGS) --studies "$(STUDIES)" $(MERGE_ARGS) $(NF_RESUME)

run-all run-all-test run-analysis-test: WORKFLOW_STAGE = core,cell_type_analysis,differential_expression,integration_benchmark,analysis
run-all: $(CORE_CONTAINER_PREREQS) ## Compose core, DE, top-gene profiles and cell reports; optional RUN_CLOCKS/RUN_BENCHMARK
run-all:
	$(RUN_NEXTFLOW) run pipeline.nf $(NF_COMMON_ARGS) --merge_single_cell --studies "$(STUDIES)" \
	  --top_age_gene_count "$(TOP_AGE_GENE_COUNT)" $(CELL_TYPE_FILTER_ARGS) --run_clocks "$(RUN_CLOCKS)" --run_benchmark "$(RUN_BENCHMARK)" \
	  --allow_partial_results "$(ALLOW_PARTIAL_RESULTS)" --gene_profile_output "$(GENE_PROFILE_DIR)" --age_clock_output "$(AGE_CLOCK_DIR)" $(CLOCK_WORKFLOW_ARGS) $(NF_RESUME)

run-all-test run-analysis-test: STUDIES = $(TEST_STUDIES)
run-all-test run-analysis-test: AGE_CLOCK_METHODS = $(VERIFY_AGE_CLOCK_METHODS)
run-all-test run-analysis-test: AGE_CLOCK_ALPHAS = $(VERIFY_AGE_CLOCK_ALPHAS)
run-all-test run-analysis-test: AGE_CLOCK_L1_RATIOS = $(VERIFY_AGE_CLOCK_L1_RATIOS)
run-all-test run-analysis-test: AGE_CLOCK_EXTERNAL_STUDIES = synthetic_study_c
run-all-test run-analysis-test: $(TEST_CONTAINER_PREREQS)
run-all-test: ## Fresh composed core plus benchmark, synthetic DE, clocks, profiles and selected cell reports
run-analysis-test: check-cell-type-analysis-prerequisites check-de-prerequisites ## Exercise the composed downstream graph using existing test merges
run-analysis-test: ANALYSIS_ONLY_ARGS = --analysis_only --merged_input "$(MERGED_INPUT)" --pseudobulk_input "$(PSEUDOBULK_INPUT)"
run-all-test run-analysis-test:
	$(if $(filter run-all-test,$@),$(BUILD_TEST_R_IMAGE),@:)
	$(RUN_NEXTFLOW) run pipeline.nf $(NF_COMMON_ARGS) --merge_single_cell --studies "$(STUDIES)" \
	  $(ANALYSIS_ONLY_ARGS) --synthetic_analysis --verify_gene_profile --run_benchmark --run_clocks --clock_cell_type "$(CELL_TYPE)" \
	  --cell_types "$(VERIFY_CELL_TYPES)" --top_age_gene_count 2 $(CLOCK_WORKFLOW_ARGS)
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_workflow_suite.py --outdir "$(OUTDIR)" --cell-type "$(CELL_TYPE)" --cell-types "$(VERIFY_CELL_TYPES)" --methods "$(AGE_CLOCK_METHODS)"

check-de-prerequisites:
	@test -f "$(PSEUDOBULK_INPUT)" || { echo "Missing pseudobulk merge: $(PSEUDOBULK_INPUT). Run make run[-test] first." >&2; exit 2; }

run-de run-de-test render-de-reports render-de-reports-test: WORKFLOW_STAGE = differential_expression
run-de run-de-test: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit DE and render per-type and combined reports
run-de run-de-test:
	$(RUN_NEXTFLOW) run differential_expression.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --allow_partial_results "$(ALLOW_PARTIAL_RESULTS)" $(NF_RESUME)

render-de-reports render-de-reports-test: CELL_TYPE =
render-de-reports render-de-reports-test: $(PYTHON_CONTAINER_PREREQS) ## Rerender existing per-type DE results; optional CELL_TYPE
render-de-reports render-de-reports-test:
	$(RUN_NEXTFLOW) run differential_expression.nf $(NF_COMMON_ARGS) --per_type_de_dir "$(PER_TYPE_DIR)" $(DE_CELL_TYPE_ARGS) $(NF_RESUME)

render-report-index: ## Rebuild the workspace report index from currently published artifacts
	python3 scripts/render_report_index.py --outdir "$(OUTDIR)"

profile-gene profile-top-age-genes profile-top-age-genes-test: WORKFLOW_STAGE = analysis
profile-gene: check-de-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Profile one gene through the shared model, calculation and rendering graph
profile-gene:
	@test -n "$(CELL_TYPE)" -a -n "$(GENE)" || { echo "Set CELL_TYPE and GENE" >&2; exit 2; }
	$(RUN_NEXTFLOW) run gene_profiles.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --per_type_de_dir "$(PER_TYPE_DIR)" \
	  --profile_cell_type "$(CELL_TYPE)" --profile_gene "$(GENE)" --gene_profile_output "$(GENE_PROFILE_DIR)" $(NF_RESUME)

check-top-age-profile-prerequisites: check-de-prerequisites
	@test -n "$(wildcard $(PER_TYPE_DIR)/*/differential_expression/merged.csv)" || { echo "Missing DE models under $(PER_TYPE_DIR). Run make run-de[-test]." >&2; exit 2; }

profile-top-age-genes-test: OUTDIR = output/test
profile-top-age-genes-test: PSEUDOBULK_INPUT = $(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad
profile-top-age-genes-test: TOP_AGE_GENE_COUNT = 2
profile-top-age-genes-test: NF_RESUME =
profile-top-age-genes profile-top-age-genes-test: check-top-age-profile-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Profile top merged age-DE genes using one model preparation per type
profile-top-age-genes profile-top-age-genes-test:
	$(RUN_NEXTFLOW) run gene_profiles.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --per_type_de_dir "$(PER_TYPE_DIR)" \
	  --gene_profile_output "$(GENE_PROFILE_DIR)" --top_age_gene_count "$(TOP_AGE_GENE_COUNT)" $(NF_RESUME)

run-gene-profile-synthetic-test: OUTDIR = output/test
run-gene-profile-synthetic-test: image-python ## Fit and render one planted synthetic gene through Nextflow
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(GENE_PROFILE_TEST_INPUT)"
	$(MAKE) profile-gene OUTDIR="$(OUTDIR)" PSEUDOBULK_INPUT="$(GENE_PROFILE_TEST_INPUT)" CELL_TYPE="CD14 monocyte" GENE=SYNTH_AGE_MARKER_1 NF_RESUME= PYTHON_CONTAINER_PREREQS=
	python3 scripts/check_workflow_status.py --outdir "$(OUTDIR)"
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_gene_profile.py --profile-dir "$(GENE_PROFILE_DIR)/cd14-monocyte/synth-age-marker-1" --gene SYNTH_AGE_MARKER_1

check-integration-benchmark-prerequisites:
	@test -f "$(INTEGRATION_BENCHMARK_INPUT)" || { echo "Missing single-cell merge: $(INTEGRATION_BENCHMARK_INPUT)" >&2; exit 2; }

run-integration-benchmark run-integration-benchmark-test: WORKFLOW_STAGE = integration_benchmark
run-integration-benchmark run-integration-benchmark-test: check-integration-benchmark-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Run the global integration and annotation benchmark
run-integration-benchmark run-integration-benchmark-test:
	$(RUN_NEXTFLOW) run integration_benchmark.nf $(NF_COMMON_ARGS) -ansi-log false --merged_input "$(INTEGRATION_BENCHMARK_INPUT)" --log_level "$(INTEGRATION_BENCHMARK_LOG_LEVEL)" $(NF_RESUME)

check-integration-benchmark-test: OUTDIR = output/test
check-integration-benchmark-test: $(PYTHON_CONTAINER_PREREQS) ## Validate integration benchmark artifacts
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_integration_benchmark.py --outdir "$(OUTDIR)"

# Standalone fixture checks remain convenient; the composed suite generates its fixture once in Nextflow.
run-de-synthetic-test: OUTDIR = output/test
run-de-synthetic-test: WORKFLOW_STAGE = differential_expression,analysis
run-de-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Fit positive synthetic DE with production inclusion filters
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad"
	$(RUN_NEXTFLOW) run differential_expression.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad" --de_test_mode false --de_synthetic_test_data true

run-age-clock run-age-clocks run-age-clock-synthetic-test: WORKFLOW_STAGE = analysis
run-age-clock: CLOCK_SELECTION_ARGS = --clock_cell_type "$(CELL_TYPE)"
run-age-clock-synthetic-test: CLOCK_SELECTION_ARGS = --clock_cell_type "$(CELL_TYPE)"
run-age-clock-synthetic-test: OUTDIR = output/test
run-age-clock-synthetic-test: PSEUDOBULK_INPUT = $(SYNTHETIC_DE_DIR)/pseudobulk_merged.h5ad
run-age-clock-synthetic-test: AGE_CLOCK_EXTERNAL_STUDIES = synthetic_study_c
run-age-clock-synthetic-test: NF_RESUME =
run-age-clock-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Fit and report the selected synthetic clock through Nextflow
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python -m pbmc_pipeline.synthetic_de --config config/pipeline.json --output "$(PSEUDOBULK_INPUT)"
	$(RUN_NEXTFLOW) run age_clocks.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --age_clock_output "$(AGE_CLOCK_DIR)" $(CLOCK_SELECTION_ARGS) $(CLOCK_WORKFLOW_ARGS)

check-age-clock-prerequisites: check-de-prerequisites
	@test -n "$(AGE_CLOCK_EXTERNAL_STUDIES)" || { echo "Set AGE_CLOCK_EXTERNAL_STUDIES" >&2; exit 2; }

run-age-clock run-age-clocks: check-age-clock-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Fit one CELL_TYPE clock or every type, with aggregation and reports
run-age-clock run-age-clocks:
	$(RUN_NEXTFLOW) run age_clocks.nf $(NF_COMMON_ARGS) --pseudobulk_input "$(PSEUDOBULK_INPUT)" --age_clock_output "$(AGE_CLOCK_DIR)" $(CLOCK_SELECTION_ARGS) $(CLOCK_WORKFLOW_ARGS) $(NF_RESUME)

check-age-clock-synthetic-test: OUTDIR = output/test
check-age-clock-synthetic-test: $(PYTHON_CONTAINER_PREREQS) ## Check synthetic clock artifacts
	python3 scripts/check_workflow_status.py --outdir "$(OUTDIR)"
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_age_clock.py --outdir "$(OUTDIR)" --cell-type "$(CELL_TYPE)" $(AGE_CLOCK_METHOD_ARGS)

aggregate-age-clocks: WORKFLOW_STAGE = analysis
aggregate-age-clocks: $(PYTHON_CONTAINER_PREREQS) ## Aggregate existing clock predictions through Nextflow without refitting
	$(RUN_NEXTFLOW) run age_clocks.nf $(NF_COMMON_ARGS) --aggregate_clocks_only --age_clock_output "$(AGE_CLOCK_DIR)" $(NF_RESUME)

check-synthetic-de-results check-synthetic-de-reports: OUTDIR = output/test
check-synthetic-de-results check-synthetic-de-reports: $(PYTHON_CONTAINER_PREREQS) ## Check synthetic DE models or downstream reports
	$(DOCKER_WORKSPACE) $(PYTHON_IMAGE) python scripts/check_synthetic_de.py $(if $(filter check-synthetic-de-results,$@),results,reports) --outdir "$(OUTDIR)"

run-cell-type-analysis run-cell-type-analysis-test split-cell-types split-cell-types-test run-cell-type run-cell-type-test render-cell-type render-cell-type-test render-cell-type-reports render-cell-type-reports-test: WORKFLOW_STAGE = cell_type_analysis
run-cell-type-analysis run-cell-type-analysis-test: check-cell-type-analysis-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render from a merge; optional CELL_TYPES selection
run-cell-type-analysis run-cell-type-analysis-test:
	@set -e; started_after=$$(python3 -c 'import time; print(time.time())'); \
	$(RUN_NEXTFLOW) run cell_type_analysis.nf $(NF_COMMON_ARGS) --merged_input "$(MERGED_INPUT)" $(CELL_TYPE_FILTER_ARGS) $(NF_RESUME); \
	if [ -n "$(NF_TEST_SUFFIX)" ]; then python3 scripts/check_cell_type_reports.py --outdir "$(OUTDIR)" --cell-types "$(CELL_TYPES)" --started-after "$$started_after"; fi

check-cell-type-analysis-test-prerequisites: OUTDIR = output/test
check-cell-type-analysis-test-prerequisites: check-cell-type-analysis-prerequisites

check-cell-type-analysis-prerequisites:
	@test -f "$(MERGED_INPUT)" || { echo "Missing single-cell merge: $(MERGED_INPUT). Run make run[-test]." >&2; exit 2; }

split-cell-types split-cell-types-test: $(PYTHON_CONTAINER_PREREQS) ## Write reusable raw-count splits from an existing merge
split-cell-types split-cell-types-test:
	$(RUN_NEXTFLOW) run cell_type_analysis.nf $(NF_COMMON_ARGS) --merged_input "$(MERGED_INPUT)" --split_only $(NF_RESUME)

check-cell-type-prerequisites check-cell-type-test-prerequisites:
	@test -f "$(CELL_TYPE_INPUT)" || { echo "Missing split: $(CELL_TYPE_INPUT). Run make split-cell-types[-test]." >&2; exit 2; }

check-cell-type-test-prerequisites check-cell-type-test-render-prerequisites: OUTDIR = output/test
run-cell-type run-cell-type-test: NF_RESUME =
run-cell-type run-cell-type-test: check-cell-type-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Analyse and render an existing split; default CELL_TYPE=cd14-monocyte
run-cell-type run-cell-type-test:
	$(RUN_NEXTFLOW) run cell_type_analysis.nf $(NF_COMMON_ARGS) --cell_type_input "$(CELL_TYPE_INPUT)" $(NF_RESUME)

check-cell-type-render-prerequisites check-cell-type-test-render-prerequisites: check-cell-type-prerequisites
	@test -f "$(CELL_TYPE_ANALYSIS_INPUT)/analysis.h5ad" || { echo "Missing analysis: $(CELL_TYPE_ANALYSIS_INPUT). Run make run-cell-type[-test]." >&2; exit 2; }

render-cell-type render-cell-type-test: NF_RESUME =
render-cell-type render-cell-type-test: check-cell-type-render-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render one existing cell-type analysis without recomputing it
render-cell-type render-cell-type-test:
	$(RUN_NEXTFLOW) run cell_type_analysis.nf $(NF_COMMON_ARGS) --cell_type_input "$(CELL_TYPE_INPUT)" --cell_type_analysis_input "$(CELL_TYPE_ANALYSIS_INPUT)"

check-cell-type-reports-prerequisites:
	@test -n "$(wildcard $(OUTDIR)/analysis/per_type/*/split.h5ad)" || { echo "Missing splits under $(OUTDIR)" >&2; exit 2; }
	@test -n "$(wildcard $(OUTDIR)/analysis/per_type/*/analysis/analysis.h5ad)" || { echo "Missing analyses under $(OUTDIR)" >&2; exit 2; }

render-cell-type-reports render-cell-type-reports-test: NF_RESUME =
render-cell-type-reports render-cell-type-reports-test: check-cell-type-reports-prerequisites $(PYTHON_CONTAINER_PREREQS) ## Render existing cell-type analyses; optional CELL_TYPES
render-cell-type-reports render-cell-type-reports-test:
	$(RUN_NEXTFLOW) run cell_type_analysis.nf $(NF_COMMON_ARGS) --cell_type_input "$(OUTDIR)/analysis/per_type/*/split.h5ad" --cell_type_analysis_input "$(OUTDIR)/analysis/per_type/*/analysis" $(CELL_TYPE_FILTER_ARGS)

$(VENV_PYTHON):
	python3 -m venv .venv

$(VENV_DEPS): $(VENV_PYTHON) pyproject.toml requirements.lock requirements.dev.lock
	$(VENV_PYTHON) -m pip install -r requirements.dev.lock
	$(VENV_PYTHON) -m pip install --no-deps --no-build-isolation -e '.[dev,qc]'
	touch $(VENV_DEPS)

install: $(VENV_DEPS) ## Install development/QC dependencies into the local .venv

lint: $(VENV_DEPS) ## Run Ruff
	$(VENV_PYTHON) -m ruff check src scripts tests
	$(VENV_PYTHON) -m ruff check --select F821 reports

workflow-lint: ## Lint the Nextflow workflow (requires Nextflow >=25.04)
	$(NXF) lint tests/workflows/configuration_preflight.nf
	$(NXF) lint main.nf
	$(NXF) lint cell_type_analysis.nf
	$(NXF) lint differential_expression.nf
	$(NXF) lint integration_benchmark.nf
	$(NXF) lint gene_profiles.nf
	$(NXF) lint age_clocks.nf
	$(NXF) lint tests/workflows/de_publication.nf
	$(NXF) lint tests/workflows/cell_type_selection.nf
	$(NXF) lint tests/workflows/analysis_failures.nf
	$(NXF) lint tests/workflows/gene_resume.nf

docs-check: ## Validate repository Markdown links and documented Make targets
	python3 scripts/check_docs.py

docs-update: ## Update generated documentation from config/pipeline.json
	python3 scripts/check_docs.py --write

# Select a test file, node, or marker for focused iteration; default runs the suite.
PYTEST_ARGS ?=
test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q $(PYTEST_ARGS)

.PHONY: test-models
test-models: $(VENV_DEPS) ## Check compatibility with locally acquired AIFI models
	$(VENV_PYTHON) -m pytest -q -m external_models

verify: STUDIES = $(TEST_STUDIES)
verify: lint test-unit test-models workflow-lint docs-check workflow-test-publication workflow-test-selection workflow-test-modes workflow-test-preflight workflow-test-failures workflow-test-main-failures workflow-test-resume site-test site-test-browser run-all-test ## Run static, frontend, and fresh Docker workflow checks
	$(MAKE) site-test-workflow

.PHONY: workflow-test-publication
workflow-test-publication: ## Exercise DE publication and staging with lightweight Nextflow stubs
	@mkdir -p output/tests
	@set -e; publication_dir=$$(mktemp -d "$(CURDIR)/output/tests/de-publication.XXXXXX"); \
	$(NXF) run tests/workflows/de_publication.nf -c nextflow.config -profile standard -stub-run \
	  -work-dir "$$publication_dir/work" --outdir "$$publication_dir/published"; \
	python3 scripts/check_de_publication.py --outdir "$$publication_dir/published"

.PHONY: workflow-test-selection
workflow-test-selection: ## Check cell-type selection and split/analysis pairing without Docker
	$(NXF) run tests/workflows/cell_type_selection.nf -c nextflow.config -profile standard -ansi-log false
	python3 scripts/check_cell_type_selection.py --nextflow "$(NXF)"

.PHONY: workflow-test-failures
workflow-test-failures: ## Prove that independent failed analysis tasks preserve sibling outputs and partial summaries
	python3 -m scripts.check_workflow_failures --nextflow "$(NXF)"

.PHONY: workflow-test-resume
workflow-test-resume: ## Check scientific caching and report-only invalidation with the real profile graph
	python3 -m scripts.check_workflow_resume --nextflow "$(NXF)"

.PHONY: workflow-test-modes
workflow-test-modes: ## Check main workflow modes, branch defaults and input preflight without scientific data
	python3 scripts/check_workflow_modes.py --nextflow "$(NXF)"

.PHONY: workflow-test-main-failures check-workflow-inventory
workflow-test-main-failures: ## Prove single-task failure isolation and DE handoffs through main.nf with stubs
	python3 -m scripts.check_main_workflow_failures --nextflow "$(NXF)"

WORKFLOW_INVENTORY_CHECK = python3 -m scripts.check_workflow_inventory --outdir "$(OUTDIR)" --mode "$(or $(PIPELINE_MODE),full)" --cell-types "$(VERIFY_CELL_TYPES)" --clock-cell-type "$(CELL_TYPE)"
check-workflow-inventory: OUTDIR = output/test
check-workflow-inventory: ## Validate current composed test tasks and every selected gene profile without rerunning science
	$(WORKFLOW_INVENTORY_CHECK)

.PHONY: workflow-failures
workflow-failures: ## Show task errors and failed DE models from the latest saved workflow; optional OUTDIR
	python3 -m scripts.show_workflow_failures --outdir "$(OUTDIR)"

.PHONY: workflow-test-preflight
workflow-test-preflight: $(VENV_DEPS) ## Check real configuration validation and fatal fan-out gating without scientific fits
	PATH="$(CURDIR)/.venv/bin:$$PATH" python3 -m scripts.check_configuration_preflight --nextflow "$(NXF)"

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
	$(NXF) lint main.nf
	$(NXF) lint cell_type_analysis.nf
	$(NXF) lint differential_expression.nf
	$(NXF) lint integration_benchmark.nf
	$(NXF) lint pipeline.nf
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

test-unit: $(VENV_DEPS) ## Run deterministic unit and output-contract tests
	$(VENV_PYTHON) -m pytest -q

.PHONY: test-models
test-models: $(VENV_DEPS) ## Check compatibility with locally acquired AIFI models
	$(VENV_PYTHON) -m pytest -q -m external_models

verify: STUDIES = $(TEST_STUDIES)
verify: lint test-unit test-models workflow-lint docs-check workflow-test-publication workflow-test-selection workflow-test-failures workflow-test-resume site-test site-test-browser run-all-test ## Run static, frontend, and fresh Docker workflow checks
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
	@mkdir -p output/tests
	@set -e; for phase in fits reports all-cells; do \
	  failure_dir=$$(mktemp -d "$(CURDIR)/output/tests/analysis-failures.XXXXXX"); \
	  case "$$phase" in \
	    fits) failure_args="--stub_fail_gene --stub_fail_clock --stub_fail_cell_type --stub_fail_de";; \
	    reports) failure_args="--stub_fail_gene_report --stub_fail_clock_report --stub_fail_cell_report --stub_fail_de_report";; \
	    all-cells) failure_args="--stub_fail_all_cell_types";; \
	  esac; \
	  python3 -m scripts.run_workflow --outdir "$$failure_dir/published" --name analysis-failures -- \
	    $(NXF) run tests/workflows/analysis_failures.nf -c nextflow.config -profile standard -stub-run -ansi-log false \
	    -work-dir "$$failure_dir/work" --outdir "$$failure_dir/published" --project_dir "$(CURDIR)" $$failure_args; \
	  python3 scripts/check_analysis_failures.py --outdir "$$failure_dir/published" --phase "$$phase"; \
	done

.PHONY: workflow-test-resume
workflow-test-resume: ## Check scientific caching and report-only invalidation with the real profile graph
	python3 -m scripts.check_workflow_resume --nextflow "$(NXF)"

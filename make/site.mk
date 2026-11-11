SITE_INPUT ?= $(OUTDIR)/analysis/per_type
SITE_PATHWAY_DIR ?= site/public/pathways
SITE_PORT ?= 8000
SITE_BIND ?= 127.0.0.1
SITE_ANALYSIS_DIR = $(patsubst %/,%,$(dir $(abspath $(SITE_INPUT))))
SITE_OUTPUT_DIR = $(patsubst %/,%,$(dir $(SITE_ANALYSIS_DIR)))
SITE_ANALYSIS_NAME = $(notdir $(SITE_ANALYSIS_DIR))
SITE_INPUT_NAME = $(notdir $(patsubst %/,%,$(SITE_INPUT)))
SITE_EXPORTER_REVISION = $(shell git rev-parse HEAD)
SITE_CONTAINER_ROOT = /opt/pbmc-ageing/site
SITE_DOCKER_USER = --user "$(shell id -u):$(shell id -g)"
SITE_DOCKER_RUN = docker run --rm $(SITE_DOCKER_USER) -e HOME=/tmp

.PHONY: site-install site-browser-install site-export site-build site-test site-test-unit site-test-browser site-test-browser-navigation site-test-browser-release site-test-browser-focus site-test-workflow serve-frontend
site-install: ## Install locked frontend dependencies locally
	cd site && npm ci --no-bin-links --cache .npm-cache --no-audit --no-fund

site-browser-install: ## Install the locked browser locally for frontend browser tests
	cd site && PLAYWRIGHT_BROWSERS_PATH="$(CURDIR)/site/.playwright-browsers" node node_modules/@playwright/test/cli.js install chromium

site-export: image-site ## Export completed per-cell-type DE results in the frontend image
	$(SITE_DOCKER_RUN) -v "$(SITE_OUTPUT_DIR):/analysis:ro" -v "$(CURDIR)/site/public:$(SITE_CONTAINER_ROOT)/public" -v "$(CURDIR)/site/data-inputs:$(SITE_CONTAINER_ROOT)/data-inputs:ro" -e SITE_EXPORTER_REVISION="$(SITE_EXPORTER_REVISION)" $(SITE_IMAGE) node scripts/release.mjs "/analysis/$(SITE_ANALYSIS_NAME)/$(SITE_INPUT_NAME)" "$(SITE_CONTAINER_ROOT)/public/data"

site-build: image-site ## Validate the existing export and build the static frontend in the image
	@echo '[site-build] Validating the existing release'
	$(SITE_DOCKER_RUN) -v "$(CURDIR)/site/public:$(SITE_CONTAINER_ROOT)/public:ro" $(SITE_IMAGE) node scripts/verify-release.mjs public/data
	@echo '[site-build] Building the frontend (type checking, then Vite)'
	mkdir -p site/dist
	$(SITE_DOCKER_RUN) -v "$(CURDIR)/site/public:$(SITE_CONTAINER_ROOT)/public:ro" -v "$(CURDIR)/site/dist:$(SITE_CONTAINER_ROOT)/dist" $(SITE_IMAGE) npm run build
	@echo '[site-build] Frontend build complete: site/dist'

site-test: image-site ## Test the frontend and export contract in the frontend image
	$(SITE_DOCKER_RUN) $(SITE_IMAGE) node scripts/export.contract.mjs
	$(SITE_DOCKER_RUN) $(SITE_IMAGE) node scripts/release.contract.mjs
	$(SITE_DOCKER_RUN) $(SITE_IMAGE) npm test

site-test-unit: image-site ## Run a focused Vitest file: SITE_TEST_FILE=src/path/file.test.ts
	@test -n "$(SITE_TEST_FILE)" || { echo 'Set SITE_TEST_FILE to a Vitest file' >&2; exit 2; }
	$(SITE_DOCKER_RUN) $(SITE_IMAGE) npm test -- "$(SITE_TEST_FILE)"

site-test-browser: image-site-browser ## Build a synthetic release and run browser journeys in the browser image
	docker run --rm $(SITE_BROWSER_IMAGE) npm run test:browser

site-test-browser-navigation: image-site-browser ## Run browser navigation and visualization journeys
	docker run --rm -e SITE_BROWSER_GROUP=navigation $(SITE_BROWSER_IMAGE) npm run test:browser

site-test-browser-release: image-site-browser ## Run browser release, loading, and failure-mode checks
	docker run --rm -e SITE_BROWSER_GROUP=release $(SITE_BROWSER_IMAGE) npm run test:browser

site-test-browser-focus: image-site-browser ## Run browser journeys matching SITE_BROWSER_FILTER
	@test -n "$(SITE_BROWSER_FILTER)" || { echo 'Set SITE_BROWSER_FILTER to part of a journey name' >&2; exit 2; }
	docker run --rm -e SITE_BROWSER_FILTER="$(SITE_BROWSER_FILTER)" $(SITE_BROWSER_IMAGE) npm run test:browser

site-test-workflow: OUTDIR = output/test
site-test-workflow: image-site-browser ## Exercise existing test workflow outputs in the browser image
	docker run --rm -v "$(abspath $(OUTDIR)):/analysis:ro" -v "$(CURDIR)/site/public/pathways:$(SITE_CONTAINER_ROOT)/public/pathways:ro" -v "$(CURDIR)/site/data-inputs:$(SITE_CONTAINER_ROOT)/data-inputs:ro" -e SITE_TEST_INPUT=/analysis/analysis/per_type $(SITE_BROWSER_IMAGE) node tests/browser/workflow.mjs

serve-frontend: image-site ## Serve the frontend with live reload; set SITE_PORT or SITE_BIND=0.0.0.0 for network access
	docker run --rm --init --interactive --tty -p "$(SITE_BIND):$(SITE_PORT):$(SITE_PORT)" -v "$(CURDIR)/site:$(SITE_CONTAINER_ROOT)" -v "$(SITE_CONTAINER_ROOT)/node_modules" -w $(SITE_CONTAINER_ROOT) $(SITE_IMAGE) node node_modules/vite/bin/vite.js --host 0.0.0.0 --port $(SITE_PORT) --strictPort

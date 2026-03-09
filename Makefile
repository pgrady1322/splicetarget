# ═══════════════════════════════════════════════════════════════════════
# splicetarget — Makefile
# ═══════════════════════════════════════════════════════════════════════

.PHONY: help install dev test test-cov lint typecheck format clean docker demo

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install: ## Install package in editable mode
	pip install -e .

dev: ## Install with dev dependencies
	pip install -e ".[dev]"

test: ## Run test suite
	pytest tests/ -v --tb=short

test-cov: ## Run tests with coverage report
	pytest tests/ -v --tb=short --cov=splicetarget --cov-report=html

lint: ## Lint with ruff
	ruff check splicetarget/ tests/

typecheck: ## Run mypy static type checking
	mypy splicetarget/

format: ## Auto-format code
	ruff format splicetarget/ tests/
	ruff check --fix splicetarget/ tests/

clean: ## Remove build artifacts and caches
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache .ruff_cache htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +

docker: ## Build Docker image
	docker build -t splicetarget:latest .

# ── Quick-start demo ──────────────────────────────────────────────
demo: ## Run pipeline demo on example data
	splicetarget run \
		--reads examples/patient_isoseq.bam \
		--reference examples/GRCh38_subset.fa \
		--annotation examples/gencode_v44_subset.gtf \
		--gene DMD \
		--outdir results/demo_run

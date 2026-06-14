.PHONY: all install test lint format build clean release

# Modern uv-based workflow. Trunk branch is `main`; versions come from git tags
# via setuptools-scm (no manual version bumps).

all: lint test

install:
	uv sync

test:
	uv run pytest -m "not perf" --cov --cov-report term-missing

perftest:
	uv run pytest -m "perf"

lint:
	uv run ruff check ngsindex tests

format:
	uv run ruff format ngsindex tests
	uv run ruff check --fix ngsindex tests

build:
	uv build

clean:
	rm -rf build dist *.egg-info .pytest_cache .ruff_cache htmlcov
	rm -f .coverage coverage.xml
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +

# Release: tag the commit, then build + publish.
#   make release version=0.2.0
release: clean
	git tag $(version)
	uv build
	uv run pytest -m "not perf"
	@echo "Built dist/. Publish with: uv publish  (and: git push origin $(version))"

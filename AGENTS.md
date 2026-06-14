# AGENTS.md — ngsindex

Utilities for reading NGS index formats (BAI, TBI, CSI, SBI). Pure-Python library.

## Layout
- `ngsindex/__init__.py` — public API: index classes, `parse_index`, `resolve_index_file`.
- `ngsindex/utils.py` — `BinReader` (binary struct reader), `DefaultDict`.
- `ngsindex/cli.py` — `summarize-index` console script.
- `tests/` — pytest suite. Binary fixtures live in `tests/test_index/` and
  `tests/test_utils/` (consumed via the `pytest-datadir` `datadir` fixture).

## Build / test / release (uv)
- Install dev env:   `uv sync`
- Run tests + cov:   `uv run pytest --cov --cov-report term-missing`
- Lint:              `uv run ruff check ngsindex tests`
- Format / autofix:  `uv run ruff format ngsindex tests && uv run ruff check --fix ngsindex tests`
- Build sdist+wheel: `uv build`

Convenience targets exist in the `Makefile` (`make test`, `make lint`, `make build`).

## Packaging
- `pyproject.toml` is the single source of build config (setuptools backend).
- Version is dynamic via **setuptools-scm** from git tags (bare numeric tags,
  e.g. `0.1.7`). Do not hand-edit a version anywhere; tag the release commit.
- `requires-python = ">=3.10"`; classifiers cover 3.10–3.13.

## Conventions
- Trunk branch is `main`.
- Type hints use PEP 604 syntax (`X | None`) and PEP 585 generics (`list[...]`).
- `pkg_resources` was removed in favor of `importlib.metadata`.

## Known dev-env caveat
The pinned `xphyle` dependency still imports `pkg_resources` at runtime, which
setuptools 81+ no longer ships, and its transitive `pokrok` sdist fails to build
on Python 3.12+ (uses the removed `configparser.SafeConfigParser`). The dev env
therefore pins Python 3.11 (`.python-version`) and `setuptools<81`
(`[dependency-groups].dev`). Remove both pins once xphyle/pokrok are modernized.

We welcome any contributions via pull requests.

All code must be written in idiomatic python 3.10+. We use [PEP 484](https://peps.python.org/pep-0484/) type hints with modern syntax (PEP 604 `X | None`, PEP 585 `list[...]`). Linting and formatting are handled by [ruff](https://docs.astral.sh/ruff/).

The project is managed with [uv](https://docs.astral.sh/uv/):

```bash
uv sync                                          # set up the dev environment
uv run pytest --cov --cov-report term-missing    # run the tests
uv run ruff check ngsindex tests                 # lint
uv run ruff format ngsindex tests                # format
```

Style-wise, we try to adhere to the Google python style guidelines. We use Google-style docstrings, which are formatted by the [Napoleon Sphinx Plugin](https://pypi.python.org/pypi/sphinxcontrib-napoleon).

We enforce the [Contributor Covenant](http://contributor-covenant.org/) code of conduct.

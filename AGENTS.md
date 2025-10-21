# Repository Guidelines

## Project Structure & Module Organization
Python interfaces live under `lib/healpy`, organized by functionality (`pixelfunc.py`, `visufunc.py`, `rotator.py`, etc.). Cython and C++ performance kernels are in `src/` (e.g., `_pixelfunc.pyx`, `_healpy_sph_transform_lib.cc`). End-to-end tests sit in `test/`, with shared fixtures in `lib/healpy/conftest.py`. Documentation sources are in `doc/`, and vendored HEALPix dependencies reside in `cextern/`. Keep ancillary scripts in `bin/` and avoid introducing new top-level directories without discussion.

## Build, Test, and Development Commands
Install editable dependencies with `python -m pip install -e .[all,test]` and ensure `cython` is available before building extensions. Produce a distribution bundle via `python -m build`. Run the full suite locally with `pytest --doctest-plus --doctest-cython --pyargs healpy && pytest test/`. To iterate on a single module, target the file: `pytest test/test_pixelfunc.py -k ring2nest`. Regenerate documentation using `make -C doc html`, which writes to `doc/.build/html`.

When working with `uv`, create the project environment with `uv venv && source .venv/bin/activate`. Install the core package via `uv pip install .` (or `uv pip install --no-editable .` if C extensions must be imported during tests). Pull in test extras with `uv pip install '.[test]'`. Run the standard suite using `uv run --no-editable --extra test pytest --doctest-plus --doctest-cython --pyargs healpy` followed by `uv run --no-editable --extra test pytest test/`.

## Coding Style & Naming Conventions
Follow PEP 8 with 4-space indentation and descriptive, lower_snake_case function names. Public APIs exported from `lib/healpy/__init__.py` should use concise, astronomy-friendly names; private helpers stay prefixed with `_`. Maintain NumPy-style docstrings, including clear parameter units. Cython modules should mirror existing naming (`_module.cpp`, `_module.pyx`) and keep the Python shim in `lib/healpy/`. When touching compiled code, match the brace and spacing conventions already present in the surrounding file.

## Testing Guidelines
Extend or create tests alongside code under `test/`, mirroring the module name (`test_pixelfunc.py` exercises `pixelfunc`). Prefer `pytest.mark.parametrize` for coverage over loops. New features should include doctest-ready examples so the `--doctest-plus` stage exercises them. When a change alters numerical tolerance, adjust the relevant assertion message and justify the threshold in a code comment. Keep optional dependencies (matplotlib, scipy) guarded with `pytest.importorskip`.

## Commit & Pull Request Guidelines
Commits follow a short, imperative summary (`Fix query_disc strip for NEST`). Group related file changes together and update `CHANGELOG.rst` for user-visible behavior. Before opening a pull request, rerun the full pytest matrix and document output or screenshots for visual routines. Create and manage PRs via the GitHub CLI (`gh pr create`, `gh pr status`) so reviewers get a consistent template. Reference related issues with GitHub keywords, describe the scientific motivation, and mention any doc or data updates. Maintain PRs in sync with `main` to avoid stale generated C files.

## Release Process
Follow the checklist in `RELEASE.md` for tagging, wheel builds, and PyPI uploads. Confirm CI artifacts match the documented version bump and update communication channels as described there.

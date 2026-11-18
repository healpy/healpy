# Steps to release `healpy`

## Github

* Review recent pull requests and update `CHANGELOG.rst`
* Create a git tag
* Draft a new release on Github using the same version name of the tag, in the description just put `See [CHANGELOG.rst](./CHANGELOG.rst)` so we don't duplicate

## PyPI - binary wheels, source

Once you publish a release in GitHub, the GitHub Actions workflow will automatically build the source package and binary wheels for all supported platforms and upload them to the Python Package Index.

## Conda packages

Conda forge should automatically detect the PyPI package and try to build the conda package,
review and merge the Pull Request at <https://github.com/conda-forge/healpy-feedstock/pulls>

## Track release with an issue

Template:

Release:
* [ ] Github: https://github.com/healpy/healpy/releases/
* [ ] PyPI: https://pypi.org/project/healpy
* [ ] conda-forge: https://anaconda.org/conda-forge/healpy
* [ ] wheels on PyPI for linux and Mac OS https://pypi.org/project/healpy/#files
* [ ] Debian packages https://packages.debian.org/sid/python3-healpy https://packages.debian.org/sid/healpy-data

@zonca @lpsinger @hivon @mreineck 

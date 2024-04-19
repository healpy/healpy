# Steps to release `healpy`

## Synchronize the C++ library

The HEALPix C++ library is available on Sourceforge under SVN, we maintain a read only mirror at <https://github.com/healpy/healpixmirror> so that we can embed into `healpy` with `git submodule`.

We should **only update the C++ sources after HEALPix has been released**, otherwise there could be incompatibilities for users that compile HEALPix separately.

Once new version of HEALPix C++ has been released, we can update `healpixmirror` with:

    git svn rebase
    git push

then in `healpy`:

    cd cextern/healpix
    git pull master
    cd ..
    git add cextern/healpix
    git commit -m "Updated HEALPix C++ to 3.5.0"

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
* [ ] wheels on PyPI for linux and Mac OS 

@zonca @lpsinger @hivon @mreineck @lpsinger 

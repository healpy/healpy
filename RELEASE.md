# Steps to release `healpy`

## Synchronize the C++ library

The HEALPix C++ library is available on Sourceforge under SVN, we maintain a read only mirror at <https://github.com/healpy/healpixmirror> so that we can embed into `healpy` with `git submodule`.

We should **only update the C++ sources after HEALPix has been released**, otherwise there could be incompatibilities for users that compile HEALPix separately.

Once new version of HEALPix C++ has been released, we can update `healpixmirror` with:

    git svn rebase
    git svn push

then in `healpy`:

    cd healpixsubmodule
    git pull master
    cd ..
    git add healpixsubmodule
    git commit -m "Updated HEALPix C++ to 3.5.0"
    
## Github

* Review recent pull requests and update `CHANGELOG.rst`
* Edit `healpy/version.py` and create a git tag
* Draft a new release on Github using the same version name of the tag, add the same details added in the Changelog to the release description

## PyPI - source

* `python setup.py build sdist`
* `twine upload dist/*`
* Attach the PyPI source package to the Github release (because that also includes the submodules and the compiled cython files, otherwise people might download the archive automatically created by Github that does not contain those)

## Conda packages

Conda forge should automatically detect the PyPI package and try to build the conda package,
review and merge the Pull Request at <https://github.com/conda-forge/healpy-feedstock/pulls>

## PyPI - Wheels

### Linux

Edit the version number in the line below and run on machine with Docker:

    mkdir -p wheelhouse && docker run --rm -v $(pwd)/wheelhouse:/wheelhouse quay.io/pypa/manylinux1_x86_64 bash -c 'for PIP in /opt/python/*/bin/pip; do $PIP install numpy==1.13.3\;python_version\<\"3.7\" numpy==1.14.3\;python_version\>=\"3.7\" && $PIP wheel --no-deps healpy==1.12.10; done; for WHEEL in *.whl; do auditwheel repair $WHEEL; done'
    
Once in a while, update the `manulinux1` docker container with:

    docker pull quay.io/pypa/manylinux1_x86_64

### macOS + MacPorts

    sudo port -N install py{27,35,36,37}-{matplotlib,numpy,six,astropy,scipy,pytest-runner,six,setuptools,pip,wheel,virtualenv} gcc8 clang-6.0
    sudo port select --set clang mp-clang-6.0
    export CC=gcc-mp-8
    export CXX=g++-mp-8
    export CFLAGS=-Wa,-q
    export CXXFLAGS=-Wa,-q
    for VERS in {2.7,3.5,3.6,3.7}; do rm -rf env && virtualenv-$VERS --system-site-packages env && env/bin/pip install --upgrade pip setuptools wheel && env/bin/pip install "numpy==1.13.3;python_version<'3.7'" "numpy==1.14.3;python_version>='3.7'" && env/bin/pip wheel --verbose --no-deps healpy==1.12.10; done
    python3.7 -m venv --system-site-packages delocate
    delocate/bin/pip install delocate
    for WHEEL in *.whl; do delocate/bin/delocate-wheel -w wheelhouse $WHEEL; done
    
## Track release with an issue

Template:

Release:
* [ ] Github: https://github.com/healpy/healpy/releases/tag/1.12.10
* [ ] PyPI: https://pypi.org/project/healpy/1.12.10/
* [ ] conda-forge: https://anaconda.org/conda-forge/healpy
* [ ] wheels on PyPI for linux and Mac OS 

```
healpy-1.12.10-cp27-cp27m-manylinux1_x86_64.whl   healpy-1.12.10-cp35-cp35m-manylinux1_x86_64.whl
healpy-1.12.10-cp27-cp27mu-manylinux1_x86_64.whl  healpy-1.12.10-cp36-cp36m-manylinux1_x86_64.whl
healpy-1.12.10-cp34-cp34m-manylinux1_x86_64.whl   healpy-1.12.10-cp37-cp37m-manylinux1_x86_64.whl
```

@zonca @DanielLenz @lpsinger @hivon @mreineck @lpsinger 

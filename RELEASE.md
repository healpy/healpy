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

## Conda packages

Conda forge should automatically detect the PyPI package and try to build the conda package,
review and merge the Pull Request at <https://github.com/conda-forge/healpy-feedstock/pulls>

## PyPI - Wheels

### Linux

Edit the version number in the line below and run on machine with Docker:

    mkdir -p wheelhouse && docker run --rm -v $(pwd)/wheelhouse:/wheelhouse quay.io/pypa/manylinux1_x86_64 bash -c 'for PIP in /opt/python/*/bin/pip; do $PIP install numpy==1.13.3\;python_version\<\"3.7\" numpy==1.14.3\;python_version\>=\"3.7\" && $PIP wheel --no-deps healpy==1.12.10; done; for WHEEL in *.whl; do auditwheel repair $WHEEL; done'

### macOS + MacPorts

    sudo port -N install py{27,35,36,37}-{matplotlib,numpy,six,astropy,scipy,pytest-runner,six,setuptools,pip,wheel,virtualenv} gcc8 clang-6.0
    sudo port select --set clang mp-clang-6.0
    export CC=gcc-mp-8
    export CXX=g++-mp-8
    export CFLAGS=-Wa,-q
    export CXXFLAGS=-Wa,-q
    for VERS in {2.7,3.5,3.6,3.7}; do rm -rf env && virtualenv-$VERS --system-site-packages env && env/bin/pip install --upgrade pip setuptools wheel && env/bin/pip install "numpy==1.13.3;python_version<'3.7'" "numpy==1.14.3;python_version>='3.7'" && env/bin/pip wheel --verbose --no-deps healpy==1.12.9; done
    python3.7 -m venv --system-site-packages delocate
    delocate/bin/pip install delocate
    for WHEEL in *.whl; do delocate/bin/delocate-wheel -w wheelhouse $WHEEL; done

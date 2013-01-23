Installation procedure for Healpy
=================================

* Note: Healpix is included, so you don't need to get it separately

Requirements
------------

Healpy needs HEALPix. You can either:
* let Healpy build its own HEALPix library from the source code included in
this package (the default behavior)
* use an existing installation :
Define the environment variable HEALPIX_EXT_PREFIX where to find the
healpix libraries and include files (eg /usr/local, so that
/usr/local/include/healpix_base.h and /usr/local/lib/libhealpix_cxx.a
exist)

Healpix needs cfitsio. You can either:
* use an existing installation :
Either define the environment variable CFITSIO_EXT_PREFIX where to find the
cfitsio library and include file (eg /usr/local, so that
/usr/local/include/fitsio.h and /usr/local/lib/libcfitsio.a exists),
or CFITSIO_EXT_INC with the include folder, e.g. /usr/local/include and 
CFITSIO_EXT_LIB with full path to the libcfitsio.* library with full filename
e.g. /usr/local/lib/libcfitsio.a or .so
* compile a specific cfitsio lib:
Define EXTERNAL_CFITSIO=no, place the  cfitsioXXXX.tar.gz in
hpbeta/libcfitsio before installing. The cfitsio version XXXX must
match the version in hpbeta/planck.make (or you need to modify it there).

the psht spherical transform library is integrated into healpy with the
pshyt cython wrapper. a pregenerated c code is included in the repository, but
if you have cython installed it will be run in the build phase.

Installation
------------

    cd healpy
    python setup.py build

OR, if you do not want OpenMP support (sometimes, it causes SegFault)

    python setup.py build --without-openmp

(alternatively, you can define the environment variable HEALPY_WITHOUT_OPENMP)

If you do not want the "-march=native" flag (if your g++ is too old)

    python setup.py build --without-native

(alternatively, you can define the environment variable HEALPY_WITHOUT_NATIVE)

If everything goes fine, you can test it:

    cd build/lib*
    ipython -pylab

>>> import healpy as H
>>> H.mollview(arange(12))
>>> pylab.show()

or run the test suite with nose:

nosetests -v

If the plot looks good, you can install:
$ sudo python setup.py install  # install in default location, need root rights
or
$ python setup.py install --install-lib=~/Softs/Python # will install healpy in directory ~/Softs/Python, which then must be in your PYTHONPATH
or
$ python setup.py install --user # will install it in your User python directory (python >= 2.6)

Compile on OSX
--------------

Suggested compilation on OSX Lion is installing pyfits, cython and cfitsio using mac ports and run:
>>> python setup.py --without-openmp

Known issues
------------

* Incompatibility with cfitisio from HEASOFT: due to a conflict of header file names it is currently not possible to use the cfitsio library provided with the HEASOFT package for compilation of Healpix C++. HEASOFT's include directory contains a file called "rotmatrix.h" which clashes with Healpix's own rotmatrix.h.

Development install
-------------------

Developers building from a snapshot of the github repository need:
  * cython > 0.14 
  * run `git submodule init` and `git submodule update` to get the healpix sources

the best way to install healpy if you plan to develop is to build the C++ extensions in place with:

python setup.py build_ext --inplace

the add the healpy/healpy folder to your PYTHONPATH

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the "build" directory,
then run:

python setup.py clean --all

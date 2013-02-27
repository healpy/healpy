Installation procedure for Healpy
=================================

Requirements
------------

healpy comes with source code for ``HEALPix`` and ``cfitsio``. If you already
have these libraries installed, then healpy will detect them using
``pkg-config`` and build and link against them. If you do not have one of them
installed, then healpy will build them itself from the bundled sources.

Installation
------------
::

    cd healpy
    python setup.py build
    sudo python setup.py install

If everything goes fine, you can test it::

    cd build/lib*
    ipython -pylab

>>> import healpy as H
>>> H.mollview(arange(12))
>>> pylab.show()

or run the test suite with nose::

    nosetests -v

If the plot looks good, you can install::

    sudo python setup.py install  # install in default location, need root rights
or::

    python setup.py install --install-lib=~/Softs/Python # will install healpy in directory ~/Softs/Python, which then must be in your PYTHONPATH
or::

    python setup.py install --user # will install it in your User python directory (python >= 2.6)

Known issues
------------

* Incompatibility with ``cfitisio`` from ``HEASOFT``: due to a conflict of header file names it is currently not possible to use the cfitsio library provided with the HEASOFT package for compilation of Healpix C++. HEASOFT's include directory contains a file called "rotmatrix.h" which clashes with Healpix's own rotmatrix.h.

* Compilation problems in the C++ package: some gcc versions (we have reports for 4.4.5 and 4.4.6) crash with an internal compiler error during compilation of libsharp. Unfortunately we have not found a workaround for this compiler problem. To our knowledge, it has been fixed in gcc 4.4.7 and in the 4.5.x and newer versions.

Development install
-------------------

Developers building from a snapshot of the github repository need:
  * cython > 0.14 
  * run ``git submodule init`` and ``git submodule update`` to get the healpix sources

the best way to install healpy if you plan to develop is to build the C++ extensions in place with::

    python setup.py build_ext --inplace

the add the healpy/healpy folder to your PYTHONPATH

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the "build" directory,
then run::

    python setup.py clean --all

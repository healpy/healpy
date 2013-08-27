Installation procedure for Healpy
=================================

Requirements
------------

Healpy depends on the Healpix C++ and cfitsio C libraries. Source code is
include with Healpy and you do not have to install them separately.

Building against external Healpix and cfitsio
---------------------------------------------

Healpy uses pkg-config to detect the presence of the Healpix and cfitsio
libraries. pkg-config is available on most systems. If you do not have
pkg-config installed, then Healpy will download and use (but not install) a
Python clone called pykg-config.

If you want to provide your own external builds of Healpix and cfitsio, then
download the following packages:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix <http://sourceforge.net/projects/healpix/>`_

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

If you are going to install the packages in a nonstandard location (say,
--prefix=/path/to/local), then you should set the environment variable
PKG_CONFIG_PATH=/path/to/local/lib/pkgconfig when building. No other
environment variable settings are necessary, and you do not need to set
PKG_CONFIG_PATH to use Healpy after you have built it.

Then, unpack each of the above packages and build them with the usual
'configure; make; make install' recipe.

Installation
------------

healpy is available on pipy, you can install it with:

::

    pip install healpy
    
otherwise, you can download a source tarball from:

https://pypi.python.org/pypi/healpy

*DO NOT DOWNLOAD* from github, github does not include the dependencies.

and build it with:

::

    cd healpy
    python setup.py build
    sudo python setup.py install

If everything goes fine, you can test it::

    cd build/lib*
    python

>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> import healpy as H
>>> H.mollview(np.arange(12))
>>> plt.show()

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
  * `autoconf` (in Ubuntu: sudo apt-get install autoconf automake libtool pkg-config)
  * `cython` > 0.14 
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

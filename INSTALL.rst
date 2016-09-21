Installation procedure for Healpy
=================================

Requirements
------------

Healpy depends on the HEALPix C++ and cfitsio C libraries. Source code for both
is included with Healpy and is built automatically, so you do not need to
install them yourself.
Only Linux and MAC OS X are supported, not Windows.

Binary installation with conda
-----------------------

The `OpenAstronomy <https://github.com/OpenAstronomy>`_ collaboration provides a `conda
channel <https://anaconda.org/openastronomy/repo>`_ with a pre-compiled version of ``healpy``
for linux 64bit and MAC OS X platforms, you can install it in Anaconda with:

    conda install -c openastronomy healpy

Source installation with Pip
---------------------------

It is possible to build the latest ``healpy`` with `pip <http://www.pip-installer.org>`_ ::

    pip install --user healpy

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade healpy

Installation on Mac OS with MacPorts
-------------------------------------------------

If you are using a Mac and have the `MacPorts <https://www.macports.org>`_
package manager, it's even easer to install Healpy with::

    sudo port install py27-healpy

Binary `apt-get` style packages are also available in the development versions of 
`Debian (sid) <https://packages.debian.org/sid/python-healpy>`_ and
`Ubuntu (utopic) <http://packages.ubuntu.com/utopic/python-healpy>`_.

Almost-as-quick installation from official source release
---------------------------------------------------------

Healpy is also available in the
`Python Package Index (PyPI) <https://pypi.python.org/pypi/healpy>`_. You can
download it with::

    curl -O https://pypi.python.org/packages/source/h/healpy/healpy-1.7.4.tar.gz

and build it with::

    tar -xzf healpy-1.7.4.tar.gz
    pushd healpy-1.7.4
    python setup.py install --user
    popd

If everything goes fine, you can test it::

    python
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> import healpy as hp 
>>> hp.mollview(np.arange(12))
>>> plt.show()

or run the test suite with nose::

    cd healpy-1.7.4 && python setup.py test

Building against external Healpix and cfitsio
---------------------------------------------

Healpy uses pkg-config to detect the presence of the Healpix and cfitsio
libraries. pkg-config is available on most systems. If you do not have
pkg-config installed, then Healpy will download and use (but not install) a
Python clone called pykg-config.

If you want to provide your own external builds of Healpix and cfitsio, then
download the following packages:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix
  <http://sourceforge.net/projects/healpix/files/Healpix_3.11/autotools_packages/>`_
  autotools-style C++ package

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

If you are going to install the packages in a nonstandard location (say,
``--prefix=/path/to/local``), then you should set the environment variable
``PKG_CONFIG_PATH=/path/to/local/lib/pkgconfig`` when building. No other
environment variable settings are necessary, and you do not need to set
``PKG_CONFIG_PATH`` to use Healpy after you have built it.

Then, unpack each of the above packages and build them with the usual
``configure; make; make install`` recipe.

Development install
-------------------

Developers building from a snapshot of the github repository need:

* ``autoconf`` and ``libtool`` (in Debian or Ubuntu:
  ``sudo apt-get install autoconf automake libtool pkg-config``)

* `cython` > 0.16

* run ``git submodule init`` and ``git submodule update`` to get the bundled
  HEALPix sources

the best way to install healpy if you plan to develop is to build the C++
extensions in place with::

    python setup.py build_ext --inplace

then add the ``healpy/healpy`` folder to your ``PYTHONPATH``.

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the ``build`` directory,
then run::

    python setup.py clean --all

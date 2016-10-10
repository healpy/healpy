====================================
Healpy, a python wrapper for healpix
====================================

.. image:: https://travis-ci.org/healpy/healpy.png?branch=master
   :target: https://travis-ci.org/healpy/healpy

.. image:: https://readthedocs.org/projects/healpy/badge/?version=latest
   :target: https://readthedocs.org/projects/healpy/?badge=latest
   :alt: Documentation Status

Description
-----------

Healpy provides a python package to manipulate healpix maps. It is based on the
standard numeric and visualisation tools for Python, Numpy and matplotlib.

To find more information about HEALPix, please visit its home page at
http://healpix.sourceforge.net/.

The documentation can be found at https://healpy.readthedocs.io, tutorial at
https://healpy.readthedocs.io/en/latest/tutorial.html.

Characteristics
---------------

* pixelisation manipulation (ang2pix, pix2ang, *etc.*)

* spherical harmonic transforms (map2alm, alm2map, synfast, anafast, *etc.*
  both for temperature and polarisation)

* plotting capabilities (mollweide and gnomonic projection)

* reading and writing of HEALPix FITS maps and alm

Requirements
------------

* `Python <http://www.python.org>`_ 2.7, 3.2, 3.3, 3.4 or 3.5

* `Numpy <http://numpy.scipy.org/>`_ (tested with version >=1.5.0)

* `Matplotlib <http://matplotlib.sourceforge.net/>`_ 

* Python development package is required for some distribution (e.g.,
  python-dev package for Ubuntu)

* `Astropy <http://www.astropy.org>`_

Quick installation with Pip
---------------------------

The quickest way to install Healpy is with `pip <http://www.pip-installer.org>`_
(>= 1.4.2), which automatically fetches the latest version of Healpy and any
missing dependencies::

    pip install --user healpy

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade healpy

See `INSTALL.rst <https://github.com/healpy/healpy/blob/master/INSTALL.rst>`_
for further details and other installation options.

Optional
--------

Healpy depends on the HEALPix C++ and cfitsio C libraries. Source code is
include with Healpy and you do not have to install them separately.

However, if you have them installed already, Healpy should detect and reuse
them instead of building them from source. To use your own installations of
HEALPix and cfitsio, you will also need:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix
  <http://sourceforge.net/projects/healpix/files/Healpix_3.11/autotools_packages/>`_
  autotools-style C++ package

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

See `INSTALL.rst <https://github.com/healpy/healpy/blob/master/INSTALL.rst>`_
for further instructions.

Known issues
------------

* Building with OpenMP support: the underlying HEALPix C++ library can be built
  to use `OpenMP <http://openmp.org/wp/>`_ to speed up some operations on
  systems with multiple cores. Most, but not all, modern C/C++ compilers support
  OpenMP, `the notable exception being clang <http://openmp.llvm.org>`_.

  If your Healpy build fails with an error message about being unable to link
  against `-lgomp`, then this typically means that Healpy detected an
  already-installed HEALPix C++ library that was built with OpenMP support, but
  you are trying to build Healpy with a compiler that does not support OpenMP.
  Try cleaning the build with `python setup.py clean --all`, and set the
  environment variables `CC` and `CXX` to point to an OpenMP-capable compiler,
  such as gcc/g++.

* Healpy does not currently support Windows.
  See https://github.com/healpy/healpy/issues/25.

* Incompatibility with ``cfitisio`` from ``HEASOFT``: due to a conflict of
  header file names it is currently not possible to use the cfitsio library
  provided with the HEASOFT package for compilation of HEALPix C++. HEASOFT's
  include directory contains a file called "rotmatrix.h" which clashes with
  HEALPix's own rotmatrix.h.

* Compilation problems in the C++ package: some gcc versions (we have reports
  for 4.4.5 and 4.4.6) crash with an internal compiler error during compilation
  of libsharp. Unfortunately we have not found a workaround for this compiler
  problem. To our knowledge, it has been fixed in gcc 4.4.7 and in the 4.5.x
  and newer versions.

* Healpy pixel functions, e.g. ``ang2pix``, do not support 32-bit platforms.
  See https://github.com/healpy/healpy/issues/194.

Support
-------

For specific *HOWTO* questions please create a question on StackOverflow_ and
tag it with the `healpy` tag, so that answers will be easily searchable on
google.

If you think you found a bug or you have install issues, open an issue on GitHub:
https://github.com/healpy/healpy/issues

For more general discussion, you can write to the healpy mailing list:
https://groups.google.com/d/forum/healpy

.. _StackOverflow: http://stackoverflow.com/questions/ask

Contribute
----------

Project development takes place on github, http://github.com/healpy/healpy,
please open an issue over there for reporting bugs or suggest improvements.
Collaboration is very welcome, just fork the project on github and send pull
requests back to the main repository.

Developers
----------
Core developers:

* Cyrille Rosset
* Andrea Zonca
* Martin Reinecke
* Leo Singer 

List of contributors: https://github.com/healpy/healpy/graphs/contributors

Acknowledgements
----------------

Note that, as stated `here <http://healpix.sourceforge.net/downloads.php>`_
publications based on work using the HEALPix software package should include
both of the following:

1. an acknowledgment statement: "Some of the results in this paper have been
   derived using the HEALPix (Górski et al., 2005) package". The complete
   reference is:

      Górski, K.M., E. Hivon, A.J. Banday, B.D. Wandelt, F.K. Hansen,
      M. Reinecke, and M. Bartelmann, HEALPix: A Framework for
      High-resolution Discretization and Fast Analysis of Data
      Distributed on the Sphere, Ap.J., 622, 759-771, 2005.

2. at the first use of the HEALPix acronym, a footnote placed in the main body
   of the paper referring to the HEALPix web site, currently
   http://healpix.sf.net

As healpy is based on HEALPix Software (the C++ library), the same
condition applies to it.

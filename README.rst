====================================
Healpy, a python wrapper for healpix
====================================

.. image:: https://badge.fury.io/py/healpy.svg
    :target: https://badge.fury.io/py/healpy

.. image:: https://anaconda.org/conda-forge/healpy/badges/version.svg
    :target: https://anaconda.org/conda-forge/healpy

.. image:: https://github.com/healpy/healpy/actions/workflows/cibuildwheel.yml/badge.svg
   :target: https://github.com/healpy/healpy/actions/workflows/cibuildwheel.yml

.. image:: https://readthedocs.org/projects/healpy/badge/?version=latest
   :target: https://readthedocs.org/projects/healpy/?badge=latest
   :alt: Documentation Status

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gist/zonca/9c114608e0903a3b8ea0bfe41c96f255/master

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2605425.svg
   :target: https://doi.org/10.5281/zenodo.2605425

.. image:: http://joss.theoj.org/papers/10.21105/joss.01298/status.svg
   :target: https://doi.org/10.21105/joss.01298

Description
-----------

`healpy` is a Python package to handle pixelated data on the sphere. It is based on the
`Hierarchical Equal Area isoLatitude Pixelization (HEALPix) <https://healpix.jpl.nasa.gov/>`_
scheme and bundles the `HEALPix` C++ library.

`HEALPix` was developed to efficiently process Cosmic Microwave Background data from Cosmology
experiments like BOOMERANG and WMAP but it is now used in other branches of Astrophysics to
store data from all-sky surveys. The target audience used to be primarily the Cosmology
scientific community but currently anyone interested in handling pixelated data on the sphere
is very welcome to propose new features.

Capabilities
------------

`healpy` provides utilities to:

* convert between sky coordinates and pixel indices in `HEALPix` nested and ring schemes
* find pixels within a disk, a polygon or a strip in the sky
* apply coordinate transformations between Galactic, Ecliptic and Equatorial reference frames
* apply custom rotations either to vectors or full maps
* read and write `HEALPix` maps to disk in FITS format
* upgrade and downgrade the resolution of existing `HEALPix` maps
* visualize maps in Mollweide, Gnomonic and Cartographic projections
* transform maps to Spherical Harmonics space and back using multi-threaded C++ routines
* compute Auto and Cross Power Spectra from maps and create map realizations from spectra

The documentation can be found at https://healpy.readthedocs.io, tutorial at
https://healpy.readthedocs.io/en/latest/tutorial.html, or execute it on `mybinder <https://mybinder.org/v2/gist/zonca/9c114608e0903a3b8ea0bfe41c96f255/master>`_

Requirements
------------

* `Python <http://www.python.org>`_ 3.10, 3.11, or 3.12

* `Numpy <http://numpy.scipy.org/>`_ (tested with version >=1.19)

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

Healpy depends on the `HEALPix` C++ and cfitsio C libraries. Source code is
include with Healpy and you do not have to install them separately.

However, if you have them installed already, Healpy should detect and reuse
them instead of building them from source. To use your own installations of
`HEALPix` and cfitsio, you will also need:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix
  <http://sourceforge.net/projects/healpix/files/Healpix_3.11/autotools_packages/>`_
  autotools-style C++ package

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

See `INSTALL.rst <https://github.com/healpy/healpy/blob/master/INSTALL.rst>`_
for further instructions.

Known issues
------------

* Building with OpenMP support: the underlying `HEALPix` C++ library can be built
  to use `OpenMP <http://openmp.org/wp/>`_ to speed up some operations on
  systems with multiple cores. Most, but not all, modern C/C++ compilers support
  OpenMP, `the notable exception being clang <http://openmp.llvm.org>`_.

  If your Healpy build fails with an error message about being unable to link
  against `-lgomp`, then this typically means that Healpy detected an
  already-installed `HEALPix` C++ library that was built with OpenMP support, but
  you are trying to build Healpy with a compiler that does not support OpenMP.
  Try cleaning the build with `python setup.py clean --all`, and set the
  environment variables `CC` and `CXX` to point to an OpenMP-capable compiler,
  such as gcc/g++.

* Healpy does not currently support Windows.
  See https://github.com/healpy/healpy/issues/25.

* Incompatibility with ``cfitisio`` from ``HEASOFT``: due to a conflict of
  header file names it is currently not possible to use the cfitsio library
  provided with the HEASOFT package for compilation of `HEALPix` C++. HEASOFT's
  include directory contains a file called "rotmatrix.h" which clashes with
  `HEALPix`'s own rotmatrix.h.

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
* Daniel Lenz

List of contributors: https://github.com/healpy/healpy/graphs/contributors

Acknowledgements
----------------

1. Cite the HEALPix and `healpy` papers, see the `CITATION file <https://github.com/healpy/healpy/blob/master/CITATION>`_ in the repository.

2. Add an acknowledgment statement: "Some of the results in this paper have been
   derived using the `healpy` and `HEALPix` packages".

3. at the first use of the `HEALPix` acronym, a footnote placed in the main body
   of the paper referring to the `HEALPix` web site, currently
   http://healpix.sf.net

Welcome to the healpy documentation
===================================

`healpy` is a Python package to handle pixelated data on the sphere. It is based on the
`Hierarchical Equal Area isoLatitude Pixelization (HEALPix) <https://healpix.jpl.nasa.gov/>`_
scheme and bundles the `HEALPix` C++ library.

`HEALPix` was developed to efficiently process Cosmic Microwave Background data from Cosmology
experiments like BOOMERANG and WMAP but it is now used in other branches of Astrophysics to
store data from all-sky surveys. The target audience used to be primarily the Cosmology
scientific community but currently anyone interested in handling pixelated data on the sphere
is very welcome to propose new features.

`healpy` provides utilities to:

* convert between sky coordinates and pixel indices in HEALPix nested and ring schemes
* find pixels within a disk, a polygon or a strip in the sky
* apply coordinate transformations between Galactic, Ecliptic and Equatorial reference frames
* apply custom rotations either to vectors or full maps
* read and write HEALPix maps to disk in FITS format
* upgrade and downgrade the resolution of existing HEALPix maps
* visualize maps in Mollweide, Gnomonic and Cartographic projections
* transform maps to Spherical Harmonics space and back using multi-threaded C++ routines
* compute Auto and Cross Power Spectra from maps and create map realizations from spectra

Changelog
---------

Review the changes in each release in the `CHANGELOG on Github <https://github.com/healpy/healpy/blob/main/CHANGELOG.rst>`_.

Citing
------

1. Cite the HEALPix and `healpy` papers, see the `CITATION file <https://github.com/healpy/healpy/blob/master/CITATION>`_ in the repository.

2. Add an acknowledgment statement: "Some of the results in this paper have been
   derived using the `healpy` and `HEALPix` packages".

3. at the first use of the `HEALPix` acronym, a footnote placed in the main body
   of the paper referring to the `HEALPix` web site, currently
   http://healpix.sf.net

Tutorial
--------

.. toctree::
   :maxdepth: 1

   tutorial
   other_tutorials

Installation
------------

.. toctree::
   :maxdepth: 1

   install

Reference
---------

.. toctree::
   :maxdepth: 2

   healpy_pix
   healpy_spht
   healpy_visu
   healpy_newvisu
   healpy_fits
   healpy_query

   healpy_rotator
   healpy_projector
   healpy_zoomtool
   healpy_line_integral_convolution

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

License
-------

.. toctree::
   :maxdepth: 1

   license

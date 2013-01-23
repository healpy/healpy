====================================
Healpy, a python wrapper for healpix
====================================

Description
-----------

Healpy provides a python package to manipulate healpix maps. It is
based on the standard numeric and visualisation tools for Python,
Numpy and matplotlib.

To find find more information about Healpix, please visit its home
page at http://healpix.jpl.nasa.gov/.

The documentation can be found at http://healpy.readthedocs.org, 
tutorial at http://healpy.readthedocs.org/en/latest/tutorial.html.

Characteristics
---------------

* pixelisation manipulation (ang2pix, pix2ang, *etc.*)

* spherical harmonic transforms (map2alm, alm2map, synfast, anafast,
  *etc.* both for temperature and polarisation)

* plotting capabilities (mollweide and gnomonic projection)

* reading and writing of Healpix FITS maps and alm

Requirements
------------

* `Python <http://www.python.org>`_, tested with 2.4, 2.5, 2.6 and
  2.7; however see `bug <http://code.google.com/p/healpy/issues/detail?id=19>`_ 
  for Python 2.4.

* `Numpy <http://numpy.scipy.org/>`_ (tested with version >=1.0.1)

* `Matplotlib <http://matplotlib.sourceforge.net/>`_ (tested with
  version >= 0.91.2 up to 1.0.1, please use latest version)

* Python development package is required for some distribution (e.g.,
  python-dev package for Ubuntu)

* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

The Healpix C++ library (from HEALPix 2.11 or HEALPix 2.20 for healpy >=
0.10) is included in the healpy package, so you don't need to get it
separately. If you do want to use your own build of HEALPix instead of
having Healpy build its bundled source, see INSTALL for further instructions.

Download
--------

The latest released version is available as a source
package at:
https://github.com/healpy/healpy/tags

Support
-------

For specific *HOWTO* questions please create a question on StackOverflow_ and tag it with the `healpy` tag, so that answers will be easily searchable on google.

If you think you found a bug or you have install issues, open an issue on github:
https://github.com/healpy/healpy/issues?state=open

For more general discussion, you can write to the healpy mailing list: https://groups.google.com/d/forum/healpy

.. _StackOverflow: http://stackoverflow.com/questions/ask

Contribute
----------

Project development takes place on github, http://github.com/healpy,
please open an issue over there for reporting bugs or suggest improvements.
Collaboration is very welcome, just fork the project on github and 
send pull requests back to the main repository.

Installation
------------

see INSTALL

Developers
----------
Core developers:

* Cyrille Rosset
* Andrea Zonca
* Martin Reinecke

Contributors:

* Yu Feng
* Duncan Hanson
* Sergey Koposov
* Maude Martin Lejeune
* Leo Singer 
* Maurizio Tomasi

Acknowledgements
----------------

Note that, as stated `here
<http://healpix.jpl.nasa.gov/healpixSoftwareGetHealpix.shtml>`_
publications based on work using the HEALPix software package should
include both of the following:

1. an acknowledgment statement: "Some of the results in this paper
   have been derived using the HEALPix (Górski et al., 2005)
   package". The complete reference is:

      Górski, K.M., E. Hivon, A.J. Banday, B.D. Wandelt, F.K. Hansen,
      M. Reinecke, and M. Bartelmann, HEALPix: A Framework for
      High-resolution Discretization and Fast Analysis of Data
      Distributed on the Sphere, Ap.J., 622, 759-771, 2005.

2. at the first use of the HEALPix acronym, a footnote placed in the
   main body of the paper referring to the HEALPix web site,
   currently http://healpix.jpl.nasa.gov

As healpy is based on HEALPix Software (the C++ library), the same
condition applies to it.

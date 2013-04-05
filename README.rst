====================================
Healpy, a python wrapper for healpix
====================================

.. image:: https://travis-ci.org/healpy/healpy.png?branch=master
   :target: https://travis-ci.org/healpy/healpy

Description
-----------

Healpy provides a python package to manipulate healpix maps. It is
based on the standard numeric and visualisation tools for Python,
Numpy and matplotlib.

To find find more information about Healpix, please visit its home
page at http://healpix.sourceforge.net/.

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

* `Numpy <http://numpy.scipy.org/>`_ (tested with version >=1.5.0)

* `Matplotlib <http://matplotlib.sourceforge.net/>`_ 

* Python development package is required for some distribution (e.g.,
  python-dev package for Ubuntu)

* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

Optional
--------

Healpy depends on the Healpix C++ and cfitsio C libraries. Source code is
include with Healpy and you do not have to install them separately.

However, if you have them installed already, Healpy should detect and reuse
them instead of building them from source. To use your own installations of
HEALPix and cfitsio, you will also need:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix <http://sourceforge.net/projects/healpix/>`_

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

See INSTALL for further instructions.

Download
--------

The latest released version is available as a source
package at:
https://pypi.python.org/pypi/healpy

`healpy` can also be automatically installed on most systems using:

$ pip install healpy --user

and upgraded with:

$ pip install --upgrade healpy --user

Known Issues
------------

healpy pixel functions, e.g. ang2pix, do not support 32bit platforms, we are working
on fixing this issue.

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

* Pierre Chanial
* Yu Feng
* Duncan Hanson
* Sergey Koposov
* Maude Martin Lejeune
* Paul Price
* Leo Singer 
* Maurizio Tomasi

Acknowledgements
----------------

Note that, as stated `here
<http://healpix.sourceforge.net/downloads.php>`_
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
   currently http://healpix.sf.net

As healpy is based on HEALPix Software (the C++ library), the same
condition applies to it.

Description
-----------

Healpy provides a python package to manipulate healpix maps. It is
based on the standard numeric and visualisation tools for Python,
Numpy and matplotlib.

To find find more information about Healpix, please visit its home
page at http://healpix.jpl.nasa.gov/.

Characteristics
---------------

* pixellisation manipulation (ang2pix, pix2ang, *etc.*)

* spherical harmonic transforms (map2alm, alm2map, synfast, anafast,
  *etc.* both for temperature and polarisation)

* plotting capabilities (mollweide and gnomonic projection)

* reading and writing of healpix fits maps and alm

Requirements
------------

* `Python <http://www.python.org>`_ (tested with 2.4, 2.5 and 2.6,
    however see `bug
    <http://code.google.com/p/healpy/issues/detail?id=19>`_ for Python
    2.4)

* `Numpy <http://numpy.scipy.org/>`_ (tested with version >=1.0.1)

* `Matplotlib <http://matplotlib.sourceforge.net/>`_ (tested with
  version >= 0.91.2 up to 1.0.1, please use latest version)

  - **Note:** matplotlib 0.98.3 has a bug preventing mollview to work
      (infinite recursion in cbook.flatten). Either correct the bug as
      indicated `there
      <http://sourceforge.net/mailarchive/message.php?msg_id=E1Kjmcj-0001UI-Ey%40dn4whf1.ch3.sourceforge.com>`_
      or update to more recent version (>= 0.98.4).


* Python development package is required for some distribution (e.g.,
  python-dev package for Ubuntu)

* `PyFITS <http://www.stsci.edu/resources/software_hardware/pyfits>`_

  - **Note:** For healpix map Nside=8192, there is a bug in Numpy
    preventing pyfits to work. You should upgrade to latest pyfits
    version (>= 2.1) which provide a workaround.

Healpix C++ library (from HEALPix 2.11 or HEALPix 2.20 for healpy >=
0.10) is included in the healpy package, so you don't need to get it
separately.


Download
--------

The latest released version is 0.10.2 and is available as a source
package at
http://code.google.com/p/healpy/downloads/list.

Installation
------------

For version >= 0.10, you must provide the cfitsio library. You can either:

* use an existing installation :

    Define the environment variable CFITSIO_EXT_PREFIX where to find the
    cfitsio library and include file (eg /usr/local, so that
    /usr/local/include/fitsio.h and /usr/local/lib/libcfitsio.a exists),
    or define CFITSIO_EXT_INC (with the cfitsio include
    directory) and CFITSIO_EXT_LIB (with the cfitsio library file name)

* compile a specific cfitsio lib:

    Define EXTERNAL_CFITSIO=no, place the  cfitsioXXXX.tar.gz in
    hpbeta/libcfitsio before installing. The cfitsio version XXXX must
    match the version in hpbeta/planck.make (or you need to modify it there).

**Note:**
  compilation on Windows has not been tested with version >= 0.10. You
  may need to tweak the setup.py file to adjust compilation
  options. If you are a Windows user, please submit a bug to describe
  a setup.py working on Windows.

Then you can follow the following instructions.

Download the last version from download tab, for example
healpy-0.9.4.tar.gz. Then::

   $ tar zxf healpy-0.9.4.tar.gz
   $ cd healpy-0.9.4
   $ python setup.py install

Use the setuptools install option to install healpy in non-standard
place, for example::

   $ python setup.py install --install-lib=~/Softs/Python

will install it in ~/Softs/Python, which you need to put into your
``PYTHONPATH``.

**Note:** 
  With python >=2.6, you can install package in a user specific
  directory with::

    $ python setup.py install --user

  This will make the package directly usable from python without any
  further configuration.

Test the installation with::

   >>> import healpy
   >>> healpy.mollview(arange(12))


Documentation
-------------

Use the online documentation, accessible with help() (or using the
 ``?`` in ipython).

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
      \M. Reinecke, and M. Bartelmann, HEALPix: A Framework for
      High-resolution Discretization and Fast Analysis of Data
      Distributed on the Sphere, Ap.J., 622, 759-771, 2005.

2. at the first use of the HEALPix acronym, a footnote placed in the
   main body of the paper referring to the HEALPix web site -
   currently http://healpix.jpl.nasa.gov

As healpy is based on HEALPix Software (the C++ library), the same
condition applies to it.


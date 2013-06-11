Release 1.6.2, 11 Jun 2013:

* ez_setup, switch from distribute to the new setuptools

Release 1.6.0, 15th March 2013:

* support for NSIDE>8192, this broke compatibility with 32bit systems
* using the new autotools based build system of healpix_cxx
* pkg-config based install for cfitsio and healpix_cxx
* common definition file for cython modules
* test build script
* new matplotlib based mollview in healpy.newvisufunc

Release 1.5.0, 16th January 2013:

* Healpix C++ sources and cython compiled files removed from the repository,
they are however added for the release tarballs
* Added back support for CFITSIO_EXT_INC and CFITSIO_EXT_LIB, but with
same definition of HealPix
* gauss_beam: gaussian beam transfer function

Release 1.4.1, 5th November 2012:

* Removed support for CFITSIO_EXT_INC and CFITSIO_EXT_LIB
* Support for linking with libcfitsio.so or libcfitsio.dyn

Release 1.4, 4th September 2012:

* Support for building using an external HealPix library, by Leo Singer
* fixes on masked array maps

Release 1.3, 21th August 2012:

* all functions covered with unit testing or doctests
* rewrote setup.py using distutils, by Leo Singer
* all functions accept and return masked arrays created with `hp.ma`
* `read_cl` and `write_cl` support polarization
* matplotlib imported only after first plotting function is called

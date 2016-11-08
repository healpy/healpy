Release 1.10.1, 8 Nov 2016

* Removed support for Python 2.6
* Implemented Lambert azimuthal equal-area projection <https://github.com/healpy/healpy/pull/354>
* Bugfix: write multiple alms <https://github.com/healpy/healpy/pull/342>
* Depend on `astropy` instead of `pyfits` <https://github.com/healpy/healpy/pull/337>

Release 1.9.1, 17 Nov 2015, Last version to support Python 2.6

* Remove C++ 11 features <https://github.com/healpy/healpy/pull/297>
* Streamlined setup.py <https://github.com/healpy/healpy/pull/298>
* Plotting fixes for Python 3 <https://github.com/healpy/healpy/pull/303>, <https://github.com/healpy/healpy/pull/304>
* Numpy 1.10 fix <https://github.com/healpy/healpy/pull/305>

Release 1.9.0, 17 Sep 2015

* updated healpix CXX to 786 (trunk) <https://github.com/healpy/healpy/pull/280>
* drop support for Python 2.6 <https://github.com/healpy/healpy/pull/268>
* option to read all fields with `read_map` <https://github.com/healpy/healpy/pull/258>
* `write_map` and `read_map` support for partial sky maps <https://github.com/healpy/healpy/pull/254>
* Allow `read_map` to also take an HDUList or HDU instance <https://github.com/healpy/healpy/issues/249>

Release 1.8.6, 23 Apr 2015

* Renamed `get_neighbours` to `get_interp_weights` <https://github.com/healpy/healpy/issues/240>
* Updated HEALPix C++ to fix bug in `query_disc` <https://github.com/healpy/healpy/issues/229>

Release 1.8.4, 16 Jan 2015

* Fixed another permission issue on install-sh

Release 1.8.3, 16 Jan 2015

* Fix permission issue in the release tarball <https://github.com/healpy/healpy/issues/220>

Release 1.8.2, 13 Jan 2015

* Several fixes in the build process
* Support for `astropy.fits` <https://github.com/healpy/healpy/pull/213>

Release 1.8.1, 22 Jun 2014 

* Added `common.pxd` to source tarball
* Check that nside is less than 2^30 <https://github.com/healpy/healpy/pull/193>

Release 1.8.0, 21 Jun 2014 

* Python 3 support <https://github.com/healpy/healpy/pull/186>
* Fixed bug in `get_interpol_ring`: <https://github.com/healpy/healpy/pull/189>
* Performance improvements in `_query_disc.pyx`: <https://github.com/healpy/healpy/pull/184>

Release 1.7.4, 26 Feb 2014 

* Fix bug for MAC OS X build <https://github.com/healpy/healpy/pull/159>

Release 1.7.3, 28 Jan 2014 

* Minor cleanup for submitting debian package

Release 1.7.2, 27 Jan 2014 

* now package does not require autotools, fixes #155

Release 1.7.1, 23 Jan 2014 

* bugfix for Anaconda/Canopy on MAC OSX #152, #153
* fixed packaging issue #154

Release 1.7.0, 14 Jan 2014 

* rewritten spherical armonics unit tests, now it uses low res maps included in the repository
* fix in HEALPix C++ build flags allows easier install on MAC-OSX and other python environments (e.g. anaconda)
* orthview: orthografic projection
* fixed bug in monopole removal in anafast

Release 1.6.3, 26 Aug 2013:

* updated C++ sources to 3.11
* verbose=True default for most functions

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

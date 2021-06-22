Release 1.15.0 22 June 2021, included in HEALPix 3.8.0

* `write_map` keeps dtype of input map array instead of float32 <https://github.com/healpy/healpy/pull/688>
* `read_map` keeps dtype of FITS file instead of upcasting to float64 <https://github.com/healpy/healpy/pull/688>
* `write_cl` uses dtype of input cl instead of float64 <https://github.com/healpy/healpy/pull/688>
* Changed all warnings to using the `logging` module, deprecated all `verbose` keywords <https://github.com/healpy/healpy/pull/693>
* Experimental `projview` function to plot maps using projections from `matplotlib` <https://github.com/healpy/healpy/pull/695>
* Flip sign for spin-0 `alm2map_spin` and `map2alm_spin` <https://github.com/healpy/healpy/issues/707>
* Support transparency in plotting with the `alpha` parameter <https://github.com/healpy/healpy/pull/696>
* Removed the note that we will change order of cl in `synfast` and `synalm`, we will leave `new=False` default <https://github.com/healpy/healpy/pull/687>
* Added convenice functions `order2npix` and `npix2order` <https://github.com/healpy/healpy/pull/685>
* Support nested maps `hp.smoothing` <https://github.com/healpy/healpy/pull/678>
* Improvements of the build system <https://github.com/healpy/healpy/pull/660> <https://github.com/healpy/healpy/pull/661>
* Automatically build wheels for Linux/MacOS on Github actions <https://github.com/healpy/healpy/pull/656>
* Drop support for Python 2.7-3.5 <https://github.com/healpy/healpy/pull/658>
* Allow OBJECT FITS header not to be a string <https://github.com/healpy/healpy/pull/665>
* Fixed indexing issue in `bl2beam` <https://github.com/healpy/healpy/pull/667>
* Fixed `map2alm_spin` bug for masked input <https://github.com/healpy/healpy/pull/651>
* Minor bugfixes: Accept None for cls in `synalm` <https://github.com/healpy/healpy/pull/711>, Get nside from length of array in `read_map` <https://github.com/healpy/healpy/pull/710>, Fix spin 0 transforms in `alm2map_spin` <https://github.com/healpy/healpy/pull/708>, Raise exception for `rotate_alm` with `complex64` inputs <https://github.com/healpy/healpy/pull/704>, Replace deprecated numpy aliases <https://github.com/healpy/healpy/pull/698>

Release 1.14.0 22 July 2020, included in HEALPix 3.70, Last release with Python 2 support

* Fixed FITS files that were left open <https://github.com/healpy/healpy/pull/631>
* Line Integral Convolution plots to plot polarization <https://github.com/healpy/healpy/pull/617>
* reworked verbose, see `hp.disable_warnings` <https://github.com/healpy/healpy/pull/630>
* increased precision in coordinate transforms <https://github.com/healpy/healpy/pull/633>
* colormaps now are not overwritten by plotting functions <https://github.com/healpy/healpy/pull/627>
* fix propagation on `mmax` in smoothing <https://github.com/healpy/healpy/pull/612>
* updated HEALPix C++ to 3.70 <https://github.com/healpy/healpy/pull/632>
* Updated to cfitsio 3.48 (used only if missing) <https://github.com/healpy/healpy/pull/597>
* Local datapath for pixel weights <https://github.com/healpy/healpy/pull/611>
* Support pixel weights for NSIDE 8192 <https://github.com/healpy/healpy/pull/595>
* Minor bugfixes <https://github.com/healpy/healpy/pull/626>, <https://github.com/healpy/healpy/pull/624>, <https://github.com/healpy/healpy/pull/618>, <https://github.com/healpy/healpy/pull/614>

Release 1.13.0 3 Dec 2019, included in HEALPix 3.60

* updated HEALPix C++ to 3.60 <https://github.com/healpy/healpy/pull/589>
* different handling of default dtype in `read_cl`, `write_cl` and `read_map` <https://github.com/healpy/healpy/pull/586>
* implemented `dist2holes`, distance from pixel center to closest invalid pixel <https://github.com/healpy/healpy/pull/581>
* allow not-power-of-2 NSIDE for RING <https://github.com/healpy/healpy/pull/584>

Release 1.12.10 9 Sep 2019

* fix overflow in nside2npix at NSIDE>8192 <https://github.com/healpy/healpy/pull/573>
* option to set UNSEEN color in plots <https://github.com/healpy/healpy/pull/551>
* option to rotate alms in place <https://github.com/healpy/healpy/pull/555>
* option to keep the FITS dtype in `read_map` <https://github.com/healpy/healpy/pull/554>
* fix compatibility with matplotlib 3 <https://github.com/healpy/healpy/pull/563> and <https://github.com/healpy/healpy/pull/566>

Release 1.12.9 21 Mar 2019, related to the `healpy` JOSS paper

* `lmax` support in `hp.pixwin` <https://github.com/healpy/healpy/pull/544>
* `use_pixel_weights` support in `hp.smoothing` <https://github.com/healpy/healpy/pull/545>
* improved test coverage <https://github.com/healpy/healpy/pull/541>
* tutorial as a Jupyter Notebook <https://github.com/healpy/healpy/blob/master/doc/healpy_tutorial.ipynb>

Release 1.12.8 7 Dec 2018, included in HEALPix 3.5.0

* Update HEALPix C++ to latest 3.5.0 commits

Release 1.12.7 6 Dec 2018

* Rebuild of broken release 1.12.6, it was built with Cython 0.26 instead of a newer version needed for Python 3.7 support

Release 1.12.6 5 Dec 2018

* Broken release due to a packaging issue
* Important bugfix that affected only 1.12.5, synfast had a fixed seed <https://github.com/healpy/healpy/pull/510>
* Updated HEALPix C++ to 3.5.0, dynamic AVX support <https://github.com/healpy/healpy/pull/514>

Release 1.12.5 13 Nov 2018

* Explicitely set Numpy version requirement to >= 1.13 <https://github.com/healpy/healpy/pull/506>
* Implemented `hp.Rotator.rotate_map_alms` and `hp.Rotator.rotate_map_pixel` to rotate maps in spherical harmonics and pixel domain <https://github.com/healpy/healpy/pull/489>

Release 1.12.4, 25 Ago 2018

* Support for Python 3.7 on PyPi
* Update minimum `healpix-cxx` version required <https://github.com/healpy/healpy/pull/478>

Release 1.12.3, 30 Giu 2018

* No changes, just fixed Unicode Error on README.rst

Release 1.12.2, 29 Giu 2018

* No changes, just fixed upload issue to PyPI

Release 1.12.1, 29 Giu 2018

* Fixed bug in polarization rotation in `hp.Rotator.rotate_map` <https://github.com/healpy/healpy/pull/459>
* Fixed packaging issue: Add six to `setup_requires` <https://github.com/healpy/healpy/pull/457>

Release 1.12.0, 12 Giu 2018

* New `hp.Rotator.rotate_map` function to change reference frame of a full map <https://github.com/healpy/healpy/pull/450>
* Implementation of pixel weights for map2alm that makes transform exact <https://github.com/healpy/healpy/pull/442>
* Change default output FITS column names to agree with other HEALPix packages <https://github.com/healpy/healpy/pull/446>
* Reformatted the Python code with black, this made a huge changeset  <https://github.com/healpy/healpy/pull/454>

Release 1.11.0, 8 Aug 2017

* Remove NSIDE restriction to be a power of 2 for RING <https://github.com/healpy/healpy/pull/377>
* Implement Coordsys2euler zyz <https://github.com/healpy/healpy/pull/399>
* Return multiple maps as a single 2D array instead of a tuple of 1D arrays <https://github.com/healpy/healpy/pull/400>
* Support for galactic cut in anafast and map2alm <https://github.com/healpy/healpy/pull/406>
* Change in write_map default behavior: <https://github.com/healpy/healpy/pull/379> and <https://github.com/healpy/healpy/pull/386>

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

* rewritten spherical harmonics unit tests, now it uses low res maps included in the repository
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

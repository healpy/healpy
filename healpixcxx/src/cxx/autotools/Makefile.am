ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libhealpix_cxx.la

src_c_utils = \
  c_utils/c_utils.c \
  c_utils/c_utils.h \
  c_utils/walltime_c.c \
  c_utils/walltime_c.h
src_libfftpack = \
  libfftpack/bluestein.h \
  libfftpack/bluestein.c \
  libfftpack/fftpack.h \
  libfftpack/fftpack.c \
  libfftpack/ls_fft.h \
  libfftpack/ls_fft.c
src_libsharp = \
  libsharp/sharp.c \
  libsharp/sharp.h \
  libsharp/sharp_cxx.h \
  libsharp/sharp_almhelpers.c \
  libsharp/sharp_almhelpers.h \
  libsharp/sharp_complex_hacks.h \
  libsharp/sharp_core.c \
  libsharp/sharp_core.h \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_geomhelpers.h \
  libsharp/sharp_internal.h \
  libsharp/sharp_lowlevel.h \
  libsharp/sharp_vecsupport.h \
  libsharp/sharp_vecutil.h \
  libsharp/sharp_ylmgen_c.c \
  libsharp/sharp_ylmgen_c.h
src_cxxsupport = \
  cxxsupport/announce.cc \
  cxxsupport/fitshandle.cc \
  cxxsupport/geom_utils.cc \
  cxxsupport/string_utils.cc \
  cxxsupport/ls_image.cc \
  cxxsupport/paramfile.cc \
  cxxsupport/pointing.cc \
  cxxsupport/rotmatrix.cc \
  cxxsupport/trafos.cc \
  cxxsupport/walltimer.cc \
  cxxsupport/wigner.cc \
  cxxsupport/error_handling.cc \
  cxxsupport/alloc_utils.h \
  cxxsupport/announce.h \
  cxxsupport/arr.h \
  cxxsupport/bstream.h \
  cxxsupport/datatypes.h \
  cxxsupport/error_handling.h \
  cxxsupport/fftpack_support.h \
  cxxsupport/fitshandle.h \
  cxxsupport/geom_utils.h \
  cxxsupport/levels_facilities.h \
  cxxsupport/ls_image.h \
  cxxsupport/lsconstants.h \
  cxxsupport/math_utils.h \
  cxxsupport/openmp_support.h \
  cxxsupport/paramfile.h \
  cxxsupport/planck_rng.h \
  cxxsupport/pointing.h \
  cxxsupport/rangeset.h \
  cxxsupport/rotmatrix.h \
  cxxsupport/safe_cast.h \
  cxxsupport/safe_ptr.h \
  cxxsupport/share_utils.h \
  cxxsupport/sse_utils_cxx.h \
  cxxsupport/string_utils.h \
  cxxsupport/trafos.h \
  cxxsupport/vec3.h \
  cxxsupport/walltimer.h \
  cxxsupport/wigner.h \
  cxxsupport/xcomplex.h
src_healpix_cxx= \
  Healpix_cxx/alm.cc \
  Healpix_cxx/alm.h \
  Healpix_cxx/alm2map_cxx_module.cc \
  Healpix_cxx/alm_fitsio.cc \
  Healpix_cxx/alm_fitsio.h \
  Healpix_cxx/alm_healpix_tools.cc \
  Healpix_cxx/alm_healpix_tools.h \
  Healpix_cxx/alm_powspec_tools.cc \
  Healpix_cxx/alm_powspec_tools.h \
  Healpix_cxx/anafast_cxx_module.cc \
  Healpix_cxx/calc_powspec_module.cc \
  Healpix_cxx/healpix_base.cc \
  Healpix_cxx/healpix_base.h \
  Healpix_cxx/healpix_data_io.cc \
  Healpix_cxx/healpix_data_io.h \
  Healpix_cxx/healpix_map.cc \
  Healpix_cxx/healpix_map.h \
  Healpix_cxx/healpix_map_fitsio.cc \
  Healpix_cxx/healpix_map_fitsio.h \
  Healpix_cxx/healpix_tables.cc \
  Healpix_cxx/healpix_tables.h \
  Healpix_cxx/hotspots_cxx_module.cc \
  Healpix_cxx/map2tga_module.cc \
  Healpix_cxx/median_filter_cxx_module.cc \
  Healpix_cxx/mult_alm_module.cc \
  Healpix_cxx/powspec.cc \
  Healpix_cxx/powspec.h \
  Healpix_cxx/powspec_fitsio.cc \
  Healpix_cxx/powspec_fitsio.h \
  Healpix_cxx/smoothing_cxx_module.cc \
  Healpix_cxx/syn_alm_cxx_module.cc \
  Healpix_cxx/udgrade_cxx_module.cc \
  Healpix_cxx/udgrade_harmonic_cxx_module.cc

EXTRA_DIST = \
  libfftpack/fftpack_inc.c \
  libsharp/sharp_core_inc.c \
  libsharp/sharp_core_inc2.c \
  libsharp/sharp_core_inchelper.c \
  cxxsupport/font_data.inc
include_HEADERS = \
  libsharp/sharp.h \
  libsharp/sharp_lowlevel.h \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_geomhelpers.h \
  cxxsupport/arr.h \
  cxxsupport/xcomplex.h \
  cxxsupport/datatypes.h \
  cxxsupport/pointing.h \
  cxxsupport/rangeset.h \
  cxxsupport/fitshandle.h \
  cxxsupport/alloc_utils.h \
  cxxsupport/math_utils.h \
  cxxsupport/error_handling.h \
  cxxsupport/lsconstants.h \
  cxxsupport/safe_cast.h \
  cxxsupport/vec3.h \
  cxxsupport/wigner.h \
  cxxsupport/rotmatrix.h \
  cxxsupport/sse_utils_cxx.h \
  cxxsupport/string_utils.h \
  Healpix_cxx/alm.h \
  Healpix_cxx/alm_fitsio.h \
  Healpix_cxx/alm_healpix_tools.h \
  Healpix_cxx/alm_powspec_tools.h \
  Healpix_cxx/healpix_base.h \
  Healpix_cxx/healpix_data_io.h \
  Healpix_cxx/healpix_map.h \
  Healpix_cxx/healpix_map_fitsio.h \
  Healpix_cxx/healpix_tables.h \
  Healpix_cxx/powspec.h \
  Healpix_cxx/powspec_fitsio.h
libhealpix_cxx_la_SOURCES = $(src_c_utils) $(src_libfftpack) $(src_libsharp) $(src_cxxsupport) $(src_healpix_cxx)
libhealpix_cxx_la_LIBADD = $(CFITSIO_LIBS)

bin_PROGRAMS = alm2map_cxx anafast_cxx calc_powspec hotspots_cxx map2tga \
  median_filter_cxx mult_alm rotalm_cxx smoothing_cxx syn_alm_cxx udgrade_cxx \
  udgrade_harmonic_cxx hpxtest
alm2map_cxx_SOURCES = Healpix_cxx/alm2map_cxx.cc
alm2map_cxx_LDADD = libhealpix_cxx.la
alm2map_cxx_LDFLAGS = -static-libtool-libs
anafast_cxx_SOURCES = Healpix_cxx/anafast_cxx.cc
anafast_cxx_LDADD = libhealpix_cxx.la
anafast_cxx_LDFLAGS = -static-libtool-libs
calc_powspec_SOURCES = Healpix_cxx/calc_powspec.cc
calc_powspec_LDADD = libhealpix_cxx.la
calc_powspec_LDFLAGS = -static-libtool-libs
hotspots_cxx_SOURCES = Healpix_cxx/hotspots_cxx.cc
hotspots_cxx_LDADD = libhealpix_cxx.la
hotspots_cxx_LDFLAGS = -static-libtool-libs
map2tga_SOURCES = Healpix_cxx/map2tga.cc
map2tga_LDADD = libhealpix_cxx.la
map2tga_LDFLAGS = -static-libtool-libs
median_filter_cxx_SOURCES = Healpix_cxx/median_filter_cxx.cc
median_filter_cxx_LDADD = libhealpix_cxx.la
median_filter_cxx_LDFLAGS = -static-libtool-libs
mult_alm_SOURCES = Healpix_cxx/mult_alm.cc
mult_alm_LDADD = libhealpix_cxx.la
mult_alm_LDFLAGS = -static-libtool-libs
rotalm_cxx_SOURCES = Healpix_cxx/rotalm_cxx.cc
rotalm_cxx_LDADD = libhealpix_cxx.la
rotalm_cxx_LDFLAGS = -static-libtool-libs
smoothing_cxx_SOURCES = Healpix_cxx/smoothing_cxx.cc
smoothing_cxx_LDADD = libhealpix_cxx.la
smoothing_cxx_LDFLAGS = -static-libtool-libs
syn_alm_cxx_SOURCES = Healpix_cxx/syn_alm_cxx.cc
syn_alm_cxx_LDADD = libhealpix_cxx.la
syn_alm_cxx_LDFLAGS = -static-libtool-libs
udgrade_cxx_SOURCES = Healpix_cxx/udgrade_cxx.cc
udgrade_cxx_LDADD = libhealpix_cxx.la
udgrade_cxx_LDFLAGS = -static-libtool-libs
udgrade_harmonic_cxx_SOURCES = Healpix_cxx/udgrade_harmonic_cxx.cc
udgrade_harmonic_cxx_LDADD = libhealpix_cxx.la
udgrade_harmonic_cxx_LDFLAGS = -static-libtool-libs
hpxtest_SOURCES=Healpix_cxx/hpxtest.cc
hpxtest_LDADD = libhealpix_cxx.la
hpxtest_LDFLAGS = -static-libtool-libs

TESTS=hpxtest

AM_CPPFLAGS = $(CFITSIO_CFLAGS) -I$(top_srcdir)/c_utils -I$(top_srcdir)/libfftpack -I$(top_srcdir)/libsharp -I$(top_srcdir)/cxxsupport

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh

TEMP1	= $(SRCROOT)/Healpix_cxx
VPATH	= $(TEMP1) $(LIBDIR)

LIBRARIES= libhealpix_cxx.a

BINARIES= syn_alm_cxx alm2map_cxx anafast_cxx map2tga udgrade_cxx \
	  hotspots_cxx calc_powspec median_filter hpxtest smoothing_cxx \
	  rotalm_cxx mult_alm

SPHERE_HEADERS= $(TEMP1)/alm.h $(TEMP1)/alm_fitsio.h \
	$(TEMP1)/alm_powspec_tools.h $(TEMP1)/powspec.h \
	$(TEMP1)/powspec_fitsio.h \
	$(TEMP1)/ylmgen.h $(TEMP1)/alm_map_tools.h

HEADERS= $(TEMP1)/healpix_base.h $(TEMP1)/healpix_map.h \
	$(TEMP1)/healpix_map_fitsio.h $(TEMP1)/alm_healpix_tools.h \
	$(TEMP1)/healpix_data_io.h $(TEMP1)/healpix_base2.h \
	$(SPHERE_HEADERS)

include $(PARAMFILE)

SPHERE_OBJ= alm_fitsio.o powspec_fitsio.o alm_powspec_tools.o powspec.o \
	alm_map_tools.o

HEALPIX_OBJ= healpix_base.o healpix_map.o healpix_map_fitsio.o \
	alm_healpix_tools.o healpix_data_io.o healpix_base2.o $(SPHERE_OBJ)

healpix_base.o: healpix_base.h libcxxsupport.a
healpix_base2.o: healpix_base.h healpix_base2.h libcxxsupport.a
healpix_map.o: healpix_base.h healpix_map.h libcxxsupport.a
healpix_map_fitsio.o: healpix_map_fitsio.h healpix_map.h healpix_base.h \
	libcxxsupport.a
alm_fitsio.o: alm_fitsio.h alm.h \
	libcxxsupport.a
healpix_data_io.o: healpix_data_io.h libcxxsupport.a
powspec_fitsio.o: powspec.h libcxxsupport.a
powspec.o: powspec.h libcxxsupport.a
alm_healpix_tools.o: healpix_base.h healpix_map.h alm.h \
	ylmgen.h alm_healpix_tools.h alm_map_tools.h \
	libfftpack.a libcxxsupport.a
alm_map_tools.o: alm.h ylmgen.h alm_map_tools.h libfftpack.a libcxxsupport.a
alm_powspec_tools.o: powspec.h alm.h alm_powspec_tools.h libcxxsupport.a

libhealpix_cxx.a: $(HEALPIX_OBJ) $(HEADERS)
	$(ARCREATE) libhealpix_cxx.a $(HEALPIX_OBJ)

syn_alm_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
syn_alm_cxx: syn_alm_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ syn_alm_cxx.o -lhealpix_cxx \
	  -lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

alm2map_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a libfftpack.a
alm2map_cxx: alm2map_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ alm2map_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

anafast_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a libfftpack.a
anafast_cxx: anafast_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ anafast_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

map2tga.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
map2tga: map2tga.o
	$(CXXL) $(CXXLFLAGS) -o $@ map2tga.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

udgrade_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
udgrade_cxx: udgrade_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ udgrade_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

hotspots_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
hotspots_cxx: hotspots_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ hotspots_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

calc_powspec.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
calc_powspec: calc_powspec.o
	$(CXXL) $(CXXLFLAGS) -o $@ calc_powspec.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

median_filter.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
median_filter: median_filter.o
	$(CXXL) $(CXXLFLAGS) -o $@ median_filter.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

hpxtest.o: libhealpix_cxx.a libcxxsupport.a
hpxtest: hpxtest.o
	$(CXXL) $(CXXLFLAGS) -o $@ hpxtest.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

smoothing_cxx.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a libfftpack.a
smoothing_cxx: smoothing_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ smoothing_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

rotalm_cxx.o: libcxxsupport.a libhealpix_cxx.a
rotalm_cxx: rotalm_cxx.o
	$(CXXL) $(CXXLFLAGS) -o $@ rotalm_cxx.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

mult_alm.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
mult_alm: mult_alm.o
	$(CXXL) $(CXXLFLAGS) -o $@ mult_alm.o -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

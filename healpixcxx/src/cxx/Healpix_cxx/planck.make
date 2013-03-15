PKG:=Healpix_cxx

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libhealpix_cxx.a
CXXBIN:=syn_alm_cxx alm2map_cxx anafast_cxx map2tga udgrade_cxx udgrade_harmonic_cxx hotspots_cxx calc_powspec median_filter_cxx hpxtest smoothing_cxx mult_alm rotalm_cxx
SPHERE_OBJ:=  powspec.o alm.o alm_powspec_tools.o
HEALPIX_OBJ:= healpix_tables.o healpix_base.o healpix_map.o alm_healpix_tools.o
IO_OBJ:= alm_fitsio.o powspec_fitsio.o healpix_data_io.o healpix_map_fitsio.o
MOD_OBJ:= syn_alm_cxx_module.o anafast_cxx_module.o alm2map_cxx_module.o map2tga_module.o udgrade_cxx_module.o udgrade_harmonic_cxx_module.o smoothing_cxx_module.o hotspots_cxx_module.o calc_powspec_module.o median_filter_cxx_module.o mult_alm_module.o
LIBOBJ:= $(SPHERE_OBJ) $(HEALPIX_OBJ) $(IO_OBJ) $(MOD_OBJ)
ALLOBJ:=$(LIBOBJ) $(CXXBIN:%=%.o)
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(LIB_$(PKG)) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)
$(CXXBIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
all_cxxbin+=$(CXXBIN)

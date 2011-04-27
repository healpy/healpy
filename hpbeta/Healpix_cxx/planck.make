PKG:=Healpix_cxx

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libhealpix_cxx.a
LIBOBJ:=alm_powspec_tools.o alm.o powspec.o healpix_tables.o healpix_base.o healpix_map.o alm_healpix_tools.o healpix_data_io.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_cxxsupport) $(HDR_libpsht) $(HDR_libfftpack) $(HDR_c_utils)

$(LIB_$(PKG)): $(LIBOBJ)

$(LIBOBJ): $(ODEP) | $(OD)_mkdir

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))

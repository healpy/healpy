PKG:=alice

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
CXXBIN:=generateTexture alice2 alice_test
CXXBIN:=$(CXXBIN:%=$(BINDIR)/%)

OBJ:=PolarizationHolder.o TextureHolder.o SoSSkyMap.o MollweideSkyMap.o color.o
ALLOBJ:=$(OBJ) generateTexture.o alice2.o alice_test.o
OBJ:=$(OBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)


ODEP:=$(HDR_$(PKG)) $(HDR_Healpix_cxx) $(HDR_cxxsupport) $(HDR_libsharp) $(HDR_libfftpack) $(HDR_c_utils)
BDEP:=$(OBJ) $(LIB_Healpix_cxx) $(LIB_cxxsupport) $(LIB_libsharp) $(LIB_libfftpack) $(LIB_c_utils) $(LIB_libcfitsio)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir

$(BINDIR)/generateTexture: $(OD)/generateTexture.o $(BDEP)
$(BINDIR)/alice2: $(OD)/alice2.o $(BDEP)
$(BINDIR)/alice_test: $(OD)/alice_test.o $(BDEP)

all_cxxbin+=$(CXXBIN)

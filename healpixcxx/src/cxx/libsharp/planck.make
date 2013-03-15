PKG:=libsharp

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libsharp.a
LIBOBJ:=sharp_ylmgen_c.o sharp.o sharp_geomhelpers.o sharp_almhelpers.o sharp_core.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_libfftpack) $(HDR_c_utils)
$(OD)/sharp_core.o: $(SD)/sharp_core_inchelper.c $(SD)/sharp_core_inc.c $(SD)/sharp_core_inc2.c

$(LIB_$(PKG)): $(LIBOBJ)

$(LIBOBJ): $(ODEP) | $(OD)_mkdir

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))

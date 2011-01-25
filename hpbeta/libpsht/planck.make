PKG:=libpsht

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libpsht.a
LIBOBJ:=ylmgen_c.o psht.o psht_geomhelpers.o psht_almhelpers.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_libfftpack) $(HDR_c_utils) $(SD)/psht_inc.c

$(LIB_$(PKG)): $(LIBOBJ)

$(LIBOBJ): $(ODEP) | $(OD)_mkdir

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))

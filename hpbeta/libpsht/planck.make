PKG:=libpsht

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libpsht.a
BIN:=psht_test psht_perftest
LIBOBJ:=ylmgen_c.o psht.o psht_geomhelpers.o psht_almhelpers.o
ALLOBJ:=$(LIBOBJ) psht_test.o psht_perftest.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_libfftpack) $(HDR_c_utils) $(SD)/psht_inc.c
BDEP:=$(LIB_$(PKG)) $(LIB_libfftpack) $(LIB_c_utils)

$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
BIN:=$(BIN:%=$(BINDIR)/%)
$(BIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
all_cbin+=$(BIN)

PKG:=c_utils

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libc_utils.a

OBJ:=c_utils.o walltime_c.o
OBJ:=$(OBJ:%=$(OD)/%)

$(OBJ): $(HDR_$(PKG)) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))

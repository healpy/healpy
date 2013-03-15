PKG:=cxxsupport

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libcxxsupport.a
SUPPORT_OBJ:= error_handling.o string_utils.o paramfile.o walltimer.o ls_image.o
MATH_OBJ:= rotmatrix.o trafos.o wigner.o pointing.o geom_utils.o
FITS_OBJ:= fitshandle.o
OBJ:=$(SUPPORT_OBJ) $(MATH_OBJ) $(FITS_OBJ) announce.o
OBJ:=$(OBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_libfftpack) $(HDR_c_utils)

$(OBJ): $(ODEP) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)

$(OD)/ls_image.o: $(SD)/font_data.inc

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))

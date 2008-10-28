TEMP1	= $(SRCROOT)/cxxsupport
VPATH	= $(TEMP1) $(INCDIR)

LIBRARIES= libcxxsupport.a

HEADERS=$(TEMP1)/cxxutils.h $(TEMP1)/arr.h $(TEMP1)/paramfile.h \
	$(TEMP1)/fitshandle.h $(TEMP1)/message_error.h $(TEMP1)/vec3.h \
	$(TEMP1)/lsconstants.h $(TEMP1)/rotmatrix.h $(TEMP1)/pointing.h \
	$(TEMP1)/planck_rng.h $(TEMP1)/geom_utils.h $(TEMP1)/simparams.h \
	$(TEMP1)/tga_image.h $(TEMP1)/xcomplex.h $(TEMP1)/trafos.h \
	$(TEMP1)/fftpack_support.h $(TEMP1)/datatypes.h \
	$(TEMP1)/openmp_support.h

include $(PARAMFILE)

MESSAGE_ERROR_H := message_error.h
LSCONSTANTS_H := lsconstants.h
DATATYPES_H := datatypes.h $(MESSAGE_ERROR_H)
CXXUTILS_H := cxxutils.h $(MESSAGE_ERROR_H) $(LSCONSTANTS_H)
SIMPARAMS_H := simparams.h $(CXXUTILS_H)
PARAMFILE_H := paramfile.h $(SIMPARAMS_H) $(CXXUTILS_H)
ARR_H := arr.h $(CXXUTILS_H)
FITSHANDLE_H := fitshandle.h fitsio.h $(ARR_H) $(DATATYPES_H)
VEC3_H := vec3.h
ROTMATRIX_H := rotmatrix.h $(CXXUTILS_H) $(VEC3_H)
POINTING_H := pointing.h $(VEC3_H) $(CXXUTILS_H)
PLANCK_RNG_H := planck_rng.h $(CXXUTILS_H)
GEOM_UTILS_H := geom_utils.h $(CXXUTILS_H) $(VEC3_H)
TGA_IMAGE_H := tga_image.h $(ARR_H)
XCOMPLEX_H := xcomplex.h
TRAFOS_H := trafos.h $(ROTMATRIX_H) $(GEOM_UTILS_H) $(POINTING_H)
FFTPACK_SUPPORT_H := fftpack_support.h $(ARR_H) $(XCOMPLEX_H)

SUPPORT_OBJ= cxxutils.o fitshandle.o rotmatrix.o simparams.o tga_image.o \
	trafos.o

cxxutils.o: $(CXXUTILS_H)
fitshandle.o: $(FITSHANDLE_H) $(CXXUTILS_H)
rotmatrix.o: $(ROTMATRIX_H) $(VEC3_H) $(LSCONSTANTS_H)
simparams.o: $(SIMPARAMS_H) $(FITSHANDLE_H)
tga_image.o: $(TGA_IMAGE_H) font_data.inc
trafos.o: $(TRAFOS_H) $(LSCONSTANTS_H)

libcxxsupport.a: $(SUPPORT_OBJ) $(HEADERS)
	$(ARCREATE) libcxxsupport.a $(SUPPORT_OBJ)

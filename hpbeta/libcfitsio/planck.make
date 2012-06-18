PKG:=libcfitsio

ifeq ($(EXTERNAL_CFITSIO),yes)

LIB_$(PKG):=$(CFITSIO_EXT_LIB)
FULL_INCLUDE+= -I$(CFITSIO_EXT_INC)

else

PACKAGE:=$(SRCROOT)/libcfitsio/cfitsio3300.tar.gz

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(BLDROOT)/libcfitsio/cfitsio

HDR_$(PKG):=fitsio.h longnam.h
HDR_$(PKG):=$(HDR_$(PKG):%=$(BLDROOT)/libcfitsio/cfitsio/%)
LIB_$(PKG):=$(LIBDIR)/libcfitsio.a

$(LIB_$(PKG)): $(SD)/* | $(OD)_mkdir $(LIBDIR)_mkdir
	cd $(BLDROOT)/libcfitsio/ && \
	rm -rf cfitsio && \
	gunzip -c $(PACKAGE) | tar xvf - && \
	cd cfitsio && \
	MAKE="$(MAKE)" FC="$(FC)" CC="$(CC)" CFLAGS="$(CCFLAGS_NO_C)" ./configure && \
	$(MAKE) && \
	cp -p libcfitsio.a $(LIBDIR)

$(HDR_$(PKG)): $(LIB_$(PKG))

endif

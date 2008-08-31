include $(PARAMFILE)

all: fitslib

fitslib: $(LIBDIR)/libcfitsio.a

PACKAGE = $(SRCROOT)/libcfitsio/cfitsio3040.tar.gz

$(LIBDIR)/libcfitsio.a: $(PACKAGE)
	rm -rf cfitsio
	gunzip -c $(PACKAGE) | tar xvf -
	cd cfitsio; \
	mv compress_alternate.c compress.c; \
	MAKE="$(MAKE)" FC="$(FC)" CC="$(CC)" CFLAGS="$(CCFLAGS_NO_C)" \
	  ./configure --includedir="$(INCDIR)" --libdir="$(LIBDIR)"; \
	$(MAKE)
	cp cfitsio/libcfitsio.a $(LIBDIR)
	cp cfitsio/fitsio.h cfitsio/longnam.h $(INCDIR)

clean:
	if test \( -d cfitsio \); then rm -rf cfitsio; fi

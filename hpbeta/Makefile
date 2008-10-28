ifndef HEALPIX_TARGET
  HEALPIX_TARGET:=$(error HEALPIX_TARGET undefined)UNDEFINED
endif

ACTIVE_SUBDIRS := Healpix_cxx alice

SUBDIRS	:= $(ACTIVE_SUBDIRS) libcfitsio cxxsupport libfftpack

LEVELS_SRC:=$(shell pwd)
export LEVELS_SRC

PARAMFILE:=$(LEVELS_SRC)/config/config.$(HEALPIX_TARGET)
export PARAMFILE

include $(PARAMFILE)

cxxsupport.all: libcfitsio.all
Healpix_cxx.all: cxxsupport.all libfftpack.all
alice.all: Healpix_cxx.all

.PHONY: doc
doc:
	cd docsrc; $(MAKE) -f planck.make cxx_doc clean

.PHONY: test
test: Healpix_cxx.all
	cd test; ./runtest.sh

ifndef HEALPIX_TARGET
  HEALPIX_TARGET:=$(error HEALPIX_TARGET undefined. Please see README.compilation for help)UNDEFINED
endif

default: compile_all
SRCROOT:=$(shell pwd)
include $(SRCROOT)/config/config.$(HEALPIX_TARGET)
include $(SRCROOT)/config/rules.common

#.SILENT:

.PHONY: test

all_hdr:=
all_lib:=
all_cbin:=
all_cxxbin:=

FULL_INCLUDE:=

include libcfitsio/planck.make
include c_utils/planck.make
include libfftpack/planck.make
include libpsht/planck.make
include cxxsupport/planck.make
include Healpix_cxx/planck.make
include docsrc/planck.make

$(all_lib): %: | $(LIBDIR)_mkdir
	@echo "#  creating library $*"
	$(ARCREATE) $@ $^

$(all_cxxbin): %: | $(BINDIR)_mkdir
	@echo "#  linking C++ binary $*"
	$(CXXL) $(CXXLFLAGS) -o $@ $^ $(CXX_EXTRALIBS)

$(all_cbin): %: | $(BINDIR)_mkdir
	@echo "#  linking C binary $*"
	$(CXXL) $(CXXLFLAGS) -o $@ $^ $(CXX_EXTRALIBS)

compile_all: $(all_lib) $(LIB_libcfitsio) hdrcopy

hdrclean:
	@if [ -d $(INCDIR) ]; then rm -rf $(INCDIR)/* ; fi

hdrcopy: | $(INCDIR)_mkdir
	@if [ "$(all_hdr)" ]; then cp -p $(all_hdr) $(INCDIR); fi

$(notdir $(all_cbin) $(all_cxxbin)) : % : $(BINDIR)/%
$(notdir $(all_lib)) libcfitsio.a : % : $(LIBDIR)/%

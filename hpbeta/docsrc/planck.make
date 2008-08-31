TEMP1	= $(SRCROOT)/docsrc
VPATH	= $(TEMP1)

include $(PARAMFILE)

cxx_doc: prep
	rm -rf $(DOCDIR)/*
	doxygen libfftpack.dox
	mv htmldoc $(DOCDIR)/libfftpack
	doxygen cxxsupport.dox
	mv htmldoc $(DOCDIR)/cxxsupport
	doxygen Healpix_cxx.dox
	mv htmldoc $(DOCDIR)/Healpix_cxx
	rm libfftpack.tag cxxsupport.tag Healpix_cxx.tag
	cp index_cxx.html $(DOCDIR)/index.html

clean:
	rm -f *.aux *.out *.toc *.log *.tag
	rm -rf htmldoc

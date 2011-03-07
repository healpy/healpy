PKG:=docsrc

docsrc_idx: | $(DOCDIR)_mkdir
	cp $(SRCROOT)/docsrc/index_cxx.html $(DOCDIR)/index.html

docsrc_code_doc: docsrc_idx
	cd $(SRCROOT)/docsrc; \
	for i in libfftpack libpsht cxxsupport Healpix_cxx; do \
	  doxygen $${i}.dox; \
	  mv htmldoc $(DOCDIR)/$${i}; \
	done; \
	rm *.tag; \

docsrc_clean:
	cd $(SRCROOT)/docsrc; \
	rm -rf htmldoc

doc: docsrc_code_doc

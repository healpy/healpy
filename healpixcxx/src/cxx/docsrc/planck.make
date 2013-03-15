PKG:=docsrc

docsrc_idx: | $(DOCDIR)_mkdir
	cp $(SRCROOT)/docsrc/index_code.html $(DOCDIR)/index.html

docsrc_code_doc: docsrc_idx
	cd $(SRCROOT)/docsrc; \
	for i in c_utils libfftpack libsharp cxxsupport Healpix_cxx; do \
	  doxygen $${i}.dox; \
	  mv htmldoc $(DOCDIR)/$${i}; \
	done; \
	rm *.tag; \

docsrc_clean:
	cd $(SRCROOT)/docsrc; \
	rm -rf htmldoc

doc: docsrc_code_doc

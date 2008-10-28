TEMP1	= $(SRCROOT)/alice
VPATH	= $(TEMP1) $(LIBDIR)

BINARIES= generateTexture test alice2 testSoSSkyMap testMollweideSkyMap testOrthogonalSkyMap

include $(PARAMFILE)

TEST_OBJ= test.o PolarizationHolder.o TextureHolder.o
TESTSOS_OBJ = PolarizationHolder.o TextureHolder.o SoSSkyMap.o testSoSSkyMap.o
TESTMOL_OBJ = PolarizationHolder.o TextureHolder.o MollweideSkyMap.o testMollweideSkyMap.o
TESTORTH_OBJ = PolarizationHolder.o TextureHolder.o testOrthogonalSkyMap.o
ALICE2_OBJ= alice2.o PolarizationHolder.o TextureHolder.o color.o SoSSkyMap.o MollweideSkyMap.o

SoSSkyMap.o: SoSSkyMap.h
MollweideSkyMap.o: MollweideSkyMap.h
PolarizationHolder.o: PolarizationHolder.h
TextureHolder.o: TextureHolder.h
color.o: color.h
test.o: alice_utils.h
alice2.o: alice_usage.h OrthogonalSkyMap.h
testOrthogonalSkyMap.o: OrthogonalSkyMap.h

generateTexture.o: libhealpix_cxx.a libcxxsupport.a libcfitsio.a
generateTexture: generateTexture.o
	$(CXXL) $(CXXLFLAGS) -o $@ generateTexture.o  -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

$(ALICE2_OBJ): libhealpix_cxx.a libcxxsupport.a libcfitsio.a
alice2: $(ALICE2_OBJ) libhealpix_cxx.a libcxxsupport.a libcfitsio.a
	$(CXXL) $(CXXLFLAGS) -o $@ $(ALICE2_OBJ) -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)

test: $(TEST_OBJ) libhealpix_cxx.a libcxxsupport.a libcfitsio.a
	$(CXXL) $(CXXLFLAGS) -o $@ $(TEST_OBJ) -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)	

testSoSSkyMap: $(TESTSOS_OBJ) libhealpix_cxx.a libcxxsupport.a libcfitsio.a
	$(CXXL) $(CXXLFLAGS) -o $@ $(TESTSOS_OBJ) -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)		

testMollweideSkyMap: $(TESTMOL_OBJ) libhealpix_cxx.a libcxxsupport.a libcfitsio.a
	$(CXXL) $(CXXLFLAGS) -o $@ $(TESTMOL_OBJ) -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)	
	
testOrthogonalSkyMap: $(TESTORTH_OBJ) OrthogonalSkyMap.h  libhealpix_cxx.a libcxxsupport.a libcfitsio.a
	$(CXXL) $(CXXLFLAGS) -o $@ $(TESTORTH_OBJ) -lhealpix_cxx \
	-lcxxsupport -lcfitsio -lfftpack $(CXX_EXTRALIBS)		

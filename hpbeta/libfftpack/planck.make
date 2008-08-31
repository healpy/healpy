TEMP1	= $(SRCROOT)/libfftpack
VPATH	= $(TEMP1) $(INCDIR)

LIBRARIES= libfftpack.a

HEADERS= $(TEMP1)/ls_fft.h

include $(PARAMFILE)


FFTPACK_OBJ= fftpack.o bluestein.o ls_fft.o

fftpack.o: fftpack.h
bluestein.o: bluestein.h fftpack.h
ls_fft.o: ls_fft.h bluestein.h fftpack.h

libfftpack.a: $(FFTPACK_OBJ) $(HEADERS)
	$(ARCREATE) libfftpack.a $(FFTPACK_OBJ)

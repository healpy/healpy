#!/bin/sh

BINPATH=../$HEALPIX_TARGET/bin

time $BINPATH/syn_alm_cxx syn_alm.par && \
time $BINPATH/alm2map_cxx alm2map.par && \
time $BINPATH/map2tga test.fits test.tga -bar -title "Synthesized Map" && \
time $BINPATH/anafast_cxx anafast.par && \
time $BINPATH/alm2map_cxx alm2map2.par && \
time $BINPATH/map2tga test2.fits test2.tga -bar -title "Reconstructed Map" && \
time $BINPATH/udgrade_cxx udgrade.par && \
time $BINPATH/map2tga test3.fits test3.tga -bar -title "Downgraded Map" && \
time $BINPATH/map2tga test3.fits test4.tga -bar -interpol -title "Downgraded, Interpolated Map" && \
time $BINPATH/alm2map_cxx alm2map3.par && \
time $BINPATH/map2tga test4.fits test5.tga -bar -title "Synthesized Map (Nside=317)" && \
xv test.tga test2.tga test3.tga test4.tga test5.tga

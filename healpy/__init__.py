"""HealPy is a package to manipulate Healpix maps (ang2pix, pix2ang) and
compute spherical harmonics tranforms on them.
"""

from version import __version__
#try:
from pixelfunc import (ang2pix,pix2ang,pix2vec,vec2pix,
                       nside2npix,npix2nside,isnsideok,
                       ring2nest, nest2ring, get_neighbours,
                       get_interp_val)

from sphtfunc import (anafast,map2alm,
                      alm2map,Alm,synalm,synfast,
                      smoothing,smoothalm,almxfl,alm2cl,
                      pixwin,alm2signal,getylm)

from _healpy_sph_transform_lib import _alm2signal
from _healpy_pixel_lib import UNSEEN

#except ImportError:
#    print "Error: healpy package relies on numpy to work"

try:
    from visufunc import (mollview,graticule,delgraticules,gnomview,
                          projplot,projscatter, projtext)
##     from healpy_visu import (mollview,UNSEEN,graticule,
##                              remove_graticules,gnomview,blink,
##                              create_moll_image,create_gnom_image,
##                              mollzoom,get_cursor_coord)
except ImportError:
    print "Warning: Cannot import visualisation tools (needs matplotlib)"

try:
    from fitsfunc import write_map,read_map,mrdfits,mwrfits,read_alm
except:
    print "Warning: Cannot import fits i/o tools (needs pyfits)"


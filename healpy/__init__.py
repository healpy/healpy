"""HealPy is a package to manipulate Healpix maps (ang2pix, pix2ang) and
compute spherical harmonics tranforms on them.
"""

import warnings

from version import __version__
#try:
from pixelfunc import (ang2pix,pix2ang,pix2vec,vec2pix,
                       nside2npix,npix2nside,isnsideok,
                       ring2nest, nest2ring, get_neighbours,
                       get_interp_val,fit_dipole,fit_monopole,
                       remove_dipole,remove_monopole,
                       get_nside,maptype,ud_grade,reorder)

from sphtfunc import (anafast,map2alm,
                      alm2map,Alm,synalm,synfast,
                      smoothing,smoothalm,almxfl,alm2cl,
                      pixwin,alm2signal,getylm,alm2map_der1)

from query_disc_func import *

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
    if visufunc.matplotlib.__version__ == '0.98,3':
        warnings.warn("Bug in matplotlib 0.98.3 prevents mollview to work\n"+
                      "You should upgrade to matplotlib 0.98.4",
                      category=warnings.ImportWarning)
except ImportError:
    warnings.warn("Warning: Cannot import visualisation tools (needs matplotlib)",
                  category=warnings.ImportWarning)

try:
    from fitsfunc import write_map,read_map,mrdfits,mwrfits,read_alm
except:
    warnings.warn("Warning: Cannot import fits i/o tools (needs pyfits)",
                  category=warnings.ImportWarning)


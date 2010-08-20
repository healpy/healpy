# 
#  This file is part of Healpy.
# 
#  Healpy is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  Healpy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Healpy; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
#  For more information about Healpy, see http://code.google.com/p/healpy
# 
"""HealPy is a package to manipulate Healpix maps (ang2pix, pix2ang) and
compute spherical harmonics tranforms on them.
"""

import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass

from version import __version__

from pixelfunc import (ma,vec2ang,ang2vec,ang2pix,
                       pix2ang,pix2vec,vec2pix,ang2vec,
                       nside2npix,npix2nside,isnsideok,
                       ring2nest, nest2ring, get_neighbours,
                       get_all_neighbours,
                       get_interp_val,fit_dipole,fit_monopole,
                       remove_dipole,remove_monopole,
                       get_nside,maptype,ud_grade,reorder)

from sphtfunc import (anafast,map2alm,
                      alm2map,Alm,synalm,synfast,
                      smoothing,smoothalm,almxfl,alm2cl,
                      pixwin,alm2signal,getylm,alm2map_der1)

from query_disc_func import *

from zoomtool import mollzoom,set_g_clim

from rotator import Rotator

from _healpy_sph_transform_lib import _alm2signal
from _healpy_pixel_lib import UNSEEN

try:
    from visufunc import (mollview,graticule,delgraticules,gnomview,
                          projplot,projscatter, projtext, cartview)
    if visufunc.matplotlib.__version__ == '0.98,3':
        warnings.warn("Bug in matplotlib 0.98.3 prevents mollview from working\n"+
                      "You should upgrade to matplotlib 0.98.4 or above",
                      category=ImportWarning)
except ImportError:
    warnings.warn("Warning: Cannot import visualisation tools (needs matplotlib)",
                  category=ImportWarning)

try:
    from fitsfunc import write_map,read_map,mrdfits,mwrfits,read_alm
except:
    warnings.warn("Warning: Cannot import fits i/o tools (needs pyfits)",
                  category=ImportWarning)


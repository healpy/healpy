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

from pixelfunc import (ma, mask_good, mask_bad,
                       ang2pix, pix2ang,
                       pix2vec, vec2pix,
                       vec2ang, ang2vec,
                       nside2npix, npix2nside, 
                       isnsideok, isnpixok,
                       ring2nest, nest2ring, reorder,
                       get_neighbours, get_all_neighbours, max_pixrad, get_interp_val,
                       fit_dipole, fit_monopole,
                       remove_dipole, remove_monopole,
                       get_nside, maptype, ud_grade, nside2resol, nside2pixarea,
                       get_map_size)

from sphtfunc import (anafast, map2alm,
                      alm2map, Alm, synalm, synfast,
                      smoothing, smoothalm, almxfl, alm2cl,
                      pixwin, alm2map_der1, gauss_beam)

try:
    from _query_disc import query_disc, query_strip, query_polygon, boundaries
except ImportError:
    warnings.warn('Warning: cannot import query disc module')
try:
    from _pixelfunc import ringinfo, pix2ring
except ImportError:
    warnings.warn('Warning: cannot import pixelfunc module')

from zoomtool import mollzoom,set_g_clim

from rotator import Rotator, vec2dir, dir2vec

try:
    from _healpy_pixel_lib import UNSEEN
except ImportError:
    warnings.warn('Warning: cannot import pixel lib module')

try:
    from pshyt import job
    from pshyt import *
except ImportError:
    warnings.warn("Warning: Cannot import pshyt module)",
                  category=ImportWarning)

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
    from fitsfunc import write_map,read_map,mrdfits,mwrfits,read_alm,write_alm,write_cl,read_cl
except:
    warnings.warn("Warning: Cannot import fits i/o tools (needs pyfits)",
                  category=ImportWarning)


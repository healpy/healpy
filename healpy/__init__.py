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

from .version import __version__

from .pixelfunc import (ma, mask_good, mask_bad,
                       ang2pix, pix2ang,
                       xyf2pix, pix2xyf,
                       pix2vec, vec2pix,
                       vec2ang, ang2vec,
                       nside2npix, npix2nside,
                       nside2order, order2nside,
                       isnsideok, isnpixok,
                       ring2nest, nest2ring, reorder,
                       get_neighbours, get_all_neighbours, max_pixrad, get_interp_val,
                       get_interp_weights,
                       fit_dipole, fit_monopole,
                       remove_dipole, remove_monopole,
                       get_nside, maptype, ud_grade, nside2resol, nside2pixarea,
                       get_map_size)

from .sphtfunc import (anafast, map2alm,
                      alm2map, Alm, synalm, synfast,
                      smoothing, smoothalm, almxfl, alm2cl,
                      pixwin, alm2map_der1, gauss_beam)

from ._query_disc import query_disc, query_strip, query_polygon, boundaries
from ._pixelfunc import ringinfo, pix2ring

from ._sphtools import rotate_alm
from ._sphtools import alm2map_spin_healpy as alm2map_spin
from ._sphtools import map2alm_spin_healpy as map2alm_spin
from .rotator import Rotator, vec2dir, dir2vec
from ._healpy_pixel_lib import UNSEEN
from .visufunc import (mollview,graticule,delgraticules,gnomview,
                      projplot,projscatter, projtext, cartview, orthview, azeqview)
from .zoomtool import mollzoom,set_g_clim
from .fitsfunc import write_map,read_map,mrdfits,mwrfits,read_alm,write_alm,write_cl,read_cl


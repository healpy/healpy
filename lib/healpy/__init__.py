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
#  For more information about Healpy, see https://github.com/healpy/healpy
#
"""healpy is a package to manipulate HEALPix maps (ang2pix, pix2ang) and
compute spherical harmonics transforms on them.
"""

import os
import sys
from importlib.metadata import version

__version__ = version("healpy")


# Copied from https://mesonbuild.com/meson-python/how-to-guides/shared-libraries.html#internal-shared-libraries
def _append_to_sharedlib_load_path():
    """Ensure the shared libraries in this package can be loaded on Windows.

    Windows lacks a concept equivalent to RPATH: Python extension modules
    cannot find DLLs installed outside the DLL search path. This function
    ensures that the location of the shared libraries distributed inside this
    Python package is in the DLL search path of the process.

    The Windows DLL search path includes the path to the object attempting
    to load the DLL: it needs to be augmented only when the Python extension
    modules and the DLLs they require are installed in separate directories.
    Cygwin does not have the same default library search path: all locations
    where the shared libraries are installed need to be added to the search
    path.

    This function is very similar to the snippet inserted into the main
    ``__init__.py`` of a package by ``delvewheel`` when it vendors external
    shared libraries.

    .. note::

        `os.add_dll_directory` is only available for Python 3.8 and later, and
        in the Conda ``python`` packages it works as advertised only for
        version 3.10 and later. For older Python versions, pre-loading the DLLs
        with `ctypes.WinDLL` may be preferred.
    """
    basedir = os.path.dirname(__file__)
    subdir = os.path.join(basedir, '..', '.healpy.mesonpy.libs')
    if os.name == 'nt':
        os.add_dll_directory(subdir)
    elif sys.platform == 'cygwin':
        os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir, subdir))


_append_to_sharedlib_load_path()

from .pixelfunc import (
    ma,
    mask_good,
    mask_bad,
    ang2pix,
    pix2ang,
    xyf2pix,
    pix2xyf,
    pix2vec,
    vec2pix,
    vec2ang,
    ang2vec,
    nside2npix,
    npix2nside,
    nside2order,
    order2nside,
    order2npix,
    npix2order,
    isnsideok,
    isnpixok,
    ring2nest,
    nest2ring,
    reorder,
    get_all_neighbours,
    max_pixrad,
    get_interp_val,
    get_interp_weights,
    fit_dipole,
    fit_monopole,
    remove_dipole,
    remove_monopole,
    get_nside,
    maptype,
    ud_grade,
    nside2resol,
    nside2pixarea,
    get_map_size,
)

from .sphtfunc import (
    anafast,
    map2alm_lsq,
    map2alm,
    alm2map,
    harmonic_ud_grade,
    effective_resolution_fwhm,
    Alm,
    synalm,
    synfast,
    smoothing,
    smoothalm,
    almxfl,
    alm2cl,
    pixwin,
    alm2map_der1,
    gauss_beam,
    bl2beam,
    beam2bl,
    check_max_nside,
    map2alm_spin,
    alm2map_spin,
    blm_gauss,
    resize_alm
)

from ._pixelfunc import ringinfo, pix2ring
from ._query_disc import query_disc, query_strip, query_polygon, boundaries

from ._sphtools import rotate_alm
from .rotator import Rotator, vec2dir, dir2vec
from ._healpy_pixel_lib import UNSEEN

try:
    import matplotlib
except ImportError:
    pass
else:
    from .visufunc import (
        mollview,
        graticule,
        delgraticules,
        gnomview,
        projplot,
        projscatter,
        projtext,
        cartview,
        orthview,
        azeqview,
    )
    from .zoomtool import mollzoom, set_g_clim
    from .newvisufunc import projview, newprojplot

from .fitsfunc import write_map, read_map, read_alm, write_alm, write_cl, read_cl
from ._masktools import dist2holes_healpy as dist2holes
from ._hotspots import hotspots_healpy as hotspots
from ._line_integral_convolution import line_integral_convolution


from .utils.deprecation import deprecated


@deprecated("1.15.0")
def disable_warnings():
    """healpy uses logging now

    This function has no effect
    """
    pass


@deprecated("1.15.0")
def enable_warnings():
    """healpy uses logging now

    This function has no effect
    """
    pass

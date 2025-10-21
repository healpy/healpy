# Wrapper around the mask_tools functions from Healpix_cxx

import numpy as np
cimport numpy as np

import cython
from .pixelfunc import maptype, npix2nside

from _common cimport Healpix_Map, RING, ndarray2map

cdef extern from "mask_tools.h":
    cdef Healpix_Map[double] dist2holes(Healpix_Map[double] &mask,
                                        double max_distance)

def fill_small_holes(mask, nside, min_size=None, min_area_arcmin2=None):
    """
    Fill holes (regions of 0s) in a HEALPix mask that are smaller than min_size (pixels)
    or min_area_arcmin2 (arcmin^2). Returns a new mask with small holes filled.
    """
    import numpy as np
    from healpy.pixelfunc import get_all_neighbours
    npix = mask.size
    mask_out = mask.copy()
    visited = np.zeros(npix, dtype=bool)
    area_per_pix = 4 * np.pi / npix * (180 * 60 / np.pi) ** 2
    hole_id = 0
    for p in range(npix):
        if mask_out[p] == 0 and not visited[p]:
            # Start a new hole
            hole_id += 1
            stack = [p]
            hole_pixels = []
            while stack:
                q = stack.pop()
                if not visited[q] and mask_out[q] == 0:
                    visited[q] = True
                    hole_pixels.append(q)
                    n = get_all_neighbours(nside, q)
                    for nn in n:
                        if nn >= 0 and not visited[nn] and mask_out[nn] == 0:
                            stack.append(nn)
            size = len(hole_pixels)
            area = size * area_per_pix
            fill = False
            if min_size is not None and size < min_size:
                fill = True
            if min_area_arcmin2 is not None and area < min_area_arcmin2:
                fill = True
            if fill:
                mask_out[hole_pixels] = 1
    return mask_out

def dist2holes_healpy(m, maxdist=np.pi, hole_min_size=None, hole_min_surf_arcmin2=None):
    """Computes the distance (in radians) from pixel center to center of
    closest invalid pixel up to a maximal distance.

    Parameters
    ----------
    m : array-like, shape (Npix,)
      The input mask.

    maxdist : float
      The maximal distance in radians. Pixel farther from this distance are not
      taken into account (default: pi).

    hole_min_size : int, optional
      Minimum hole size (in pixels) to ignore. Holes smaller than this will be filled before distance calculation.

    hole_min_surf_arcmin2 : float, optional
      Minimum hole surface (in arcmin^2) to ignore. Holes smaller than this will be filled before distance calculation.

    Returns
    -------
    distances : out map of distances (in radians) in RING scheme as numpy arrays

    Example
    -------
    >>> import healpy as hp
    >>> import numpy as np
    >>> nside = 16
    >>> hp.dist2holes(np.random.randint(0, 2, 12*nside**2), hole_min_size=10)
    array([0.        , 0.        , 0.        , ..., 0.05831086, 0.05831086,
       0.        ])

    """
    info = maptype(m)
    if info == 0:
        mi = m.astype(np.float64, order='C', copy=True)
    elif info == 1:
        mi = m[0].astype(np.float64, order='C', copy=True)
    else:
        raise ValueError("Wrong input map (must be a valid healpix map)")

    # Optionally fill small holes before distance calculation
    if hole_min_size is not None or hole_min_surf_arcmin2 is not None:
        nside = npix2nside(mi.size)
        mi = fill_small_holes(mi, nside, hole_min_size, hole_min_surf_arcmin2)

    # View the ndarray as a Healpix_Map
    M = ndarray2map(mi, RING)

    # Prepare distance array
    npix = mi.size
    distances = np.empty(npix, dtype=np.float64)
    D = ndarray2map(distances, RING)

    D[0] = dist2holes(M[0], maxdist)

    del M, D
    return distances

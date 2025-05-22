# Wrapper around the mask_tools functions from Healpix_cxx

import numpy as np
cimport numpy as np

import cython
from .pixelfunc import maptype

from _common cimport Healpix_Map, RING, ndarray2map

cdef extern from "mask_tools.h":
    cdef Healpix_Map[double] dist2holes(Healpix_Map[double] &mask,
                                        double max_distance)

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
        try:
            from scipy.ndimage import label
        except ImportError:
            raise ImportError("scipy is required for hole filtering. Please install scipy.")
        nside = int(np.sqrt(mi.size // 12))
        # Identify holes (regions of 0s surrounded by 1s)
        mask = (mi == 0)
        structure = np.ones((3,), dtype=int)  # 1D connectivity for flat array
        labeled, num_features = label(mask, structure=structure)
        # Compute sizes
        sizes = np.bincount(labeled.ravel())
        # Area per pixel in arcmin^2
        area_per_pix = 4 * np.pi / mi.size * (180*60/np.pi)**2
        fill = np.zeros_like(mi, dtype=bool)
        for i in range(1, num_features+1):
            size = sizes[i]
            area = size * area_per_pix
            if (hole_min_size is not None and size < hole_min_size) or \
               (hole_min_surf_arcmin2 is not None and area < hole_min_surf_arcmin2):
                fill |= (labeled == i)
        mi[fill] = 1.0

    # View the ndarray as a Healpix_Map
    M = ndarray2map(mi, RING)

    # Prepare distance array
    npix = mi.size
    distances = np.empty(npix, dtype=np.float64)
    D = ndarray2map(distances, RING)

    D[0] = dist2holes(M[0], maxdist)

    del M, D
    return distances

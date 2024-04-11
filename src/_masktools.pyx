# Wrapper around the mask_tools functions from Healpix_cxx

import numpy as np
cimport numpy as np

import cython
from .pixelfunc import maptype

from _common cimport Healpix_Map, RING, ndarray2map

cdef extern from "mask_tools.h":
    cdef Healpix_Map[double] dist2holes(Healpix_Map[double] &mask,
                                        double max_distance)

def dist2holes_healpy(m, maxdist=np.pi):
    """Computes the distance (in radians) from pixel center to center of
    closest invalid pixel up to a maximal distance.

    Parameters
    ----------
    m : array-like, shape (Npix,)
      The input mask.

    maxdist : float
      The maximal distance in radians. Pixel farther from this distance are not
      taken into account (default: pi).

    Returns
    -------
    distances : out map of distances (in radians) in RING scheme as numpy arrays

    Example
    -------
    >>> import healpy as hp
    >>> import numpy as np
    >>> nside = 16
    >>> hp.dist2holes(np.random.randint(0, 2, 12*nside**2))
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

    # View the ndarray as a Healpix_Map
    M = ndarray2map(mi, RING)

    # Prepare distance array
    npix = mi.size
    distances = np.empty(npix, dtype=np.float64)
    D = ndarray2map(distances, RING)

    D[0] = dist2holes(M[0], maxdist)

    del M, D
    return distances

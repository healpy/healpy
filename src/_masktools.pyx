# Wrapper around the mask_tools functions from Healpix_cxx

import numpy as np
cimport numpy as np

import cython
from collections import deque
from operator import index as _operator_index

from .pixelfunc import get_all_neighbours, maptype, npix2nside

from _common cimport Healpix_Map, RING, ndarray2map

cdef extern from "mask_tools.h":
    cdef Healpix_Map[double] dist2holes(Healpix_Map[double] &mask,
                                        double max_distance)

def fill_small_holes(mask, min_size=None, min_area_arcmin2=None):
    """
    Fill holes (pixels marked as ``True``) in a binary HEALPix mask that are
    smaller than the requested thresholds.

    Parameters
    ----------
    mask : array_like
        Binary mask in RING ordering where ``True`` (or any non-zero value)
        marks masked pixels and ``False`` marks valid pixels. Masked arrays are
        converted using :func:`numpy.ma.getmaskarray`.
    min_size : int, optional
        Minimum number of connected pixels that should be preserved. Holes with
        fewer pixels are filled.
    min_area_arcmin2 : float, optional
        Minimum hole area, expressed in arcmin^2, that should be preserved.
        Area is computed assuming equal-area HEALPix pixels.

    Returns
    -------
    ndarray
        Boolean copy of ``mask`` with holes below the thresholds filled.

    Notes
    -----
    NSIDE is inferred from the mask length; a :class:`ValueError` is raised if
    the length does not match a valid HEALPix resolution. Connectivity uses the
    HEALPix neighbour graph in RING ordering. Pixels are considered part of a
    hole when they evaluate to ``True`` in the input mask.
    """
    if min_size is None and min_area_arcmin2 is None:
        if np.ma.isMaskedArray(mask):
            return np.ascontiguousarray(np.ma.getmaskarray(mask)).copy()
        return np.ascontiguousarray(np.asarray(mask, dtype=np.bool_)).copy()

    if min_size is not None:
        try:
            min_size = _operator_index(min_size)
        except TypeError as err:
            raise TypeError("min_size must be an integer") from err
        if min_size < 0:
            raise ValueError("min_size must be non-negative")

    if min_area_arcmin2 is not None:
        min_area_arcmin2 = float(min_area_arcmin2)
        if min_area_arcmin2 < 0:
            raise ValueError("min_area_arcmin2 must be non-negative")

    if np.ma.isMaskedArray(mask):
        mask_arr = np.ma.getmaskarray(mask)
    else:
        mask_arr = np.asarray(mask, dtype=np.bool_)
    mask_arr = np.ascontiguousarray(mask_arr, dtype=np.bool_)

    if mask_arr.ndim != 1:
        raise ValueError("mask must be one-dimensional")

    npix = mask_arr.size
    try:
        nside = npix2nside(npix)
    except ValueError as err:
        raise ValueError(
            "mask length must correspond to a valid HEALPix NSIDE"
        ) from err

    mask_out = mask_arr.copy()
    visited = np.zeros(npix, dtype=np.bool_)

    area_per_pix = None
    if min_area_arcmin2 is not None:
        area_per_pix = 4.0 * np.pi / npix * (180.0 * 60.0 / np.pi) ** 2

    for p in range(npix):
        if mask_out[p] and not visited[p]:
            stack = deque([p])
            hole_pixels = []
            while stack:
                q = stack.pop()
                if visited[q]:
                    continue
                visited[q] = True
                if not mask_out[q]:
                    continue
                hole_pixels.append(q)
                neighbours = get_all_neighbours(nside, q)
                for nn in neighbours:
                    if nn >= 0 and not visited[nn] and mask_out[nn]:
                        stack.append(nn)

            hole_size = len(hole_pixels)
            fill = False
            if min_size is not None and hole_size < min_size:
                fill = True
            if min_area_arcmin2 is not None:
                hole_area = hole_size * area_per_pix
                if hole_area < min_area_arcmin2:
                    fill = True
            if fill:
                mask_out[hole_pixels] = False

    return mask_out

def dist2holes_healpy(m, maxdist=np.pi, hole_min_size=None, hole_min_surf_arcmin2=None):
    """Computes the distance (in radians) from pixel center to center of
    closest invalid pixel up to a maximal distance.

    Parameters
    ----------
    m : array-like, shape (Npix,)
      Input mask in RING ordering where ``True`` (or non-zero entries) mark
      masked pixels and ``False`` (or zero) mark valid pixels.

    maxdist : float, optional
      The maximal distance in radians. Pixel farther from this distance are not
      taken into account (default: pi).

    hole_min_size : int, optional
      Minimum connected hole size, expressed as a number of pixels, to retain.
      Holes smaller than this threshold are filled before the distance
      calculation. Defaults to ``None`` (no size-based filtering).

    hole_min_surf_arcmin2 : float, optional
      Minimum hole area, expressed in square arcminutes, to retain. Holes with a
      smaller area than this threshold are filled before the distance
      calculation. Defaults to ``None`` (no area-based filtering).
      The filtering uses HEALPix neighbour connectivity in RING ordering.

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
    if hole_min_size is not None:
        try:
            hole_min_size = _operator_index(hole_min_size)
        except TypeError as err:
            raise TypeError("hole_min_size must be an integer") from err
        if hole_min_size < 0:
            raise ValueError("hole_min_size must be non-negative")

    if hole_min_surf_arcmin2 is not None:
        hole_min_surf_arcmin2 = float(hole_min_surf_arcmin2)
        if hole_min_surf_arcmin2 < 0:
            raise ValueError("hole_min_surf_arcmin2 must be non-negative")

    info = maptype(m)
    if info == 0:
        raw_mask = np.ma.getmaskarray(m) if np.ma.isMaskedArray(m) else np.asarray(m)
    elif info == 1:
        raw_mask = (
            np.ma.getmaskarray(m[0])
            if np.ma.isMaskedArray(m[0])
            else np.asarray(m[0])
        )
    else:
        raise ValueError("Wrong input map (must be a valid healpix map)")

    mask_bool = np.ascontiguousarray(np.asarray(raw_mask, dtype=np.bool_))
    if mask_bool.ndim != 1:
        raise ValueError("mask must be one-dimensional")

    # Optionally fill small holes before distance calculation
    if hole_min_size is not None or hole_min_surf_arcmin2 is not None:
        mask_bool = fill_small_holes(mask_bool, hole_min_size, hole_min_surf_arcmin2)

    mi = (~mask_bool).astype(np.float64, order='C', copy=True)

    # Optionally fill small holes before distance calculation
    # (already handled above)

    # View the ndarray as a Healpix_Map
    M = ndarray2map(mi, RING)

    # Prepare distance array
    npix = mi.size
    distances = np.empty(npix, dtype=np.float64)
    D = ndarray2map(distances, RING)

    D[0] = dist2holes(M[0], maxdist)

    del M, D
    return distances

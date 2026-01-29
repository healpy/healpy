# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
cimport cython

from _common cimport int64, pointing, rangeset, vec3, Healpix_Ordering_Scheme, RING, NEST, SET_NSIDE, T_Healpix_Base
import healpy.pixelfunc
from ._pixelfunc import isnsideok

@cython.boundscheck(False)
@cython.wraparound(False)
def query_disc(nside, vec, radius, inclusive = False, fact = 4, nest = False, np.ndarray[np.int64_t, ndim=1] buff=None, return_ranges=False):
    """Returns pixels whose centers lie within the disk defined by
    *vec* and *radius* (in radians) (if *inclusive* is False), or which
    overlap with this disk (if *inclusive* is True).

    Parameters
    ----------
    nside : int
      The nside of the Healpix map.
    vec : float, sequence of 3 elements
      The coordinates of unit vector defining the disk center.
    radius : float
      The radius (in radians) of the disk
    inclusive : bool, optional
      If False, return the exact set of pixels whose pixel centers lie 
      within the disk; if True, return all pixels that overlap with the disk,
      and maybe a few more. Default: False
    fact : int, optional
      Only used when inclusive=True. The overlapping test will be done at
      the resolution fact*nside. For NESTED ordering, fact must be a power of 2, less than 2**30,
      else it can be any positive integer. Default: 4.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    buff: int array, optional
      if provided, this numpy array is used to contain the return values and must be
      at least long enough to do so. Cannot be used with return_ranges=True.
    return_ranges: bool, optional
      if True, return a 2D array of pixel ranges with shape (num_ranges, 2), 
      where each row contains [start, end) of a range. This is more memory-efficient
      for queries returning large numbers of pixels. Default: False

    Returns
    -------
    ipix : int, array
      If return_ranges is False: 1D array of pixel indices.
      If return_ranges is True: 2D array of shape (num_ranges, 2) where each row
      is [start, end) representing a range of pixels.

    Note
    ----
    The performance characteristics of RING vs. NESTED schemes depend on the 
    query parameters. For simple disks with inclusive=False, RING scheme may 
    be faster. For inclusive=True, NESTED scheme often performs better due to 
    more efficient pixel indexing.
    
    For inclusive=True, the algorithm may return some pixels which don't overlap
    with the disk at all. The higher fact is chosen, the fewer false positives
    are returned, at the cost of increased run time.
    """
    # Validate parameters
    if return_ranges and buff is not None:
        raise ValueError("Cannot use both 'buff' and 'return_ranges' parameters")
    
    # Check Nside value
    if not isnsideok(nside, nest):
        raise ValueError('Wrong nside value, must be a power of 2, less than 2**30')
    cdef vec3 v = vec3(vec[0], vec[1], vec[2])
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    cdef rangeset[int64] pixset
    cdef int factor = fact
    cdef pointing ptg = pointing(v)
    ptg.normalize()
    if inclusive:
        factor = abs(fact)
        if nest and (factor == 0 or (factor & (factor - 1) != 0)):
            raise ValueError('fact must be a power of 2, less than 2**30 when '
                             'nest is True (fact=%d)' % (fact))
        hb.query_disc_inclusive(ptg, radius, pixset, factor)
    else:
        hb.query_disc(ptg, radius, pixset)

    if return_ranges:
        return rangeset_to_array(pixset)
    else:
        return pixset_to_array(pixset, buff)


def query_polygon(nside, vertices, inclusive = False, fact = 4, nest = False, np.ndarray[np.int64_t, ndim=1] buff=None, return_ranges=False):
    """ Returns the pixels whose centers lie within the convex polygon 
    defined by the *vertices* array (if *inclusive* is False), or which 
    overlap with this polygon (if *inclusive* is True).

    Parameters
    ----------
    nside : int
      The nside of the Healpix map.
    vertices : float, array-like
      Vertex array containing the vertices of the polygon, shape (N, 3).
    inclusive : bool, optional
      If False, return the exact set of pixels whose pixel centers lie
      within the polygon; if True, return all pixels that overlap with the
      polygon, and maybe a few more. Default: False.
    fact : int, optional
      Only used when inclusive=True. The overlapping test will be done at
      the resolution fact*nside. For NESTED ordering, fact must be a power of 2, less than 2**30,
      else it can be any positive integer. Default: 4.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    buff: int array, optional
      if provided, this numpy array is used to contain the return values and must be
      at least long enough to do so. Cannot be used with return_ranges=True.
    return_ranges: bool, optional
      if True, return a 2D array of pixel ranges with shape (num_ranges, 2), 
      where each row contains [start, end) of a range. This is more memory-efficient
      for queries returning large numbers of pixels. Default: False
      
    Returns
    -------
    ipix : int, array
      If return_ranges is False: 1D array of pixel indices.
      If return_ranges is True: 2D array of shape (num_ranges, 2) where each row
      is [start, end) representing a range of pixels.

    Note
    ----
    The performance characteristics of RING vs. NESTED schemes depend on the 
    query parameters. For simple polygons with inclusive=False, RING scheme may 
    be faster. For inclusive=True, NESTED scheme often performs better due to 
    more efficient pixel indexing, especially for complex polygons.
    
    For inclusive=True, the algorithm may return some pixels which don't overlap
    with the polygon at all. The higher fact is chosen, the fewer false positives
    are returned, at the cost of increased run time.
    """
    # Validate parameters
    if return_ranges and buff is not None:
        raise ValueError("Cannot use both 'buff' and 'return_ranges' parameters")
    
    # Check Nside value
    if not isnsideok(nside, nest):
        raise ValueError('Wrong nside value, must be a power of 2, less than 2**30')
    # Create vector of vertices
    cdef vector[pointing] vert
    cdef pointing ptg
    for v in vertices:
        ptg = pointing(vec3(v[0], v[1], v[2]))
        ptg.normalize()
        vert.push_back(ptg)
    # Create the Healpix_Base2 structure
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    # Call query_polygon
    cdef rangeset[int64] pixset
    cdef int factor
    if inclusive:
        factor = abs(fact)
        if nest and (factor == 0 or (factor & (factor - 1) != 0)):
            raise ValueError('fact must be a power of 2, less than 2**30 when '
                             'nest is True (fact=%d)' % (fact))
        hb.query_polygon_inclusive(vert, pixset, factor)
    else:
        hb.query_polygon(vert, pixset)

    if return_ranges:
        return rangeset_to_array(pixset)
    else:
        return pixset_to_array(pixset, buff)

def query_strip(nside, theta1, theta2, inclusive = False, nest = False, np.ndarray[np.int64_t, ndim=1] buff=None, return_ranges=False):
    """Returns pixels whose centers lie within the colatitude range
    defined by *theta1* and *theta2* (if inclusive is False), or which 
    overlap with this region (if *inclusive* is True). If theta1<theta2, the
    region between both angles is considered, otherwise the regions 
    0<theta<theta2 and theta1<theta<pi.
    
    Parameters
    ----------
    nside : int
      The nside of the Healpix map.
    theta1 : float
      First colatitude (radians)
    theta2 : float
      Second colatitude (radians)
    inclusive : bool
      If False, return the exact set of pixels whose pixels centers lie 
      within the region; if True, return all pixels that overlap with the
      region.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    buff: int array, optional
      if provided, this numpy array is used to contain the return values and must be
      at least long enough to do so. Cannot be used with return_ranges=True.
    return_ranges: bool, optional
      if True, return a 2D array of pixel ranges with shape (num_ranges, 2), 
      where each row contains [start, end) of a range. This is more memory-efficient
      for queries returning large numbers of pixels. Default: False
      
    Returns
    -------
    ipix : int, array
      If return_ranges is False: 1D array of pixel indices.
      If return_ranges is True: 2D array of shape (num_ranges, 2) where each row
      is [start, end) representing a range of pixels.
      
    Note
    ----
    When using nest=True with return_ranges=True, the implementation must materialize
    all pixels in memory to convert from RING to NESTED ordering before creating ranges.
    This defeats the memory efficiency purpose of return_ranges. For large queries with
    NESTED ordering, consider using RING ordering or the standard return format.
    """
    # Validate parameters
    if return_ranges and buff is not None:
        raise ValueError("Cannot use both 'buff' and 'return_ranges' parameters")
    
    # Check Nside value
    if not isnsideok(nside, nest):
        raise ValueError('Wrong nside value, must be a power of 2, less than 2**30')
    # Create the Healpix_Base2 structure
    cdef Healpix_Ordering_Scheme scheme_internal
    cdef T_Healpix_Base[int64] hb
    cdef rangeset[int64] pixset
    cdef np.ndarray[np.int64_t, ndim=1] result_pixels
    cdef np.ndarray[np.int64_t, ndim=2] result_ranges

    # Always call with RING ordering internally, as query_strip_internal does not support NESTED
    scheme_internal = RING
    hb = T_Healpix_Base[int64](nside, scheme_internal, SET_NSIDE)
    hb.query_strip(theta1, theta2, inclusive, pixset)
    
    if return_ranges:
        result_ranges = rangeset_to_array(pixset)
        if nest:
            # If original request was for NESTED, convert the ranges
            # Each pixel in the range needs to be converted
            result_pixels = healpy.pixelfunc.ring2nest(nside, 
                np.concatenate([np.arange(result_ranges[i, 0], result_ranges[i, 1]) 
                               for i in range(result_ranges.shape[0])]))
            # Reconstruct ranges from the converted pixels
            # This is inefficient but necessary for NESTED ordering with query_strip
            # Sort the pixels first
            result_pixels.sort()
            # Find consecutive ranges
            ranges_list = []
            if len(result_pixels) > 0:
                start = result_pixels[0]
                prev = result_pixels[0]
                for pix in result_pixels[1:]:
                    if pix != prev + 1:
                        ranges_list.append([start, prev + 1])
                        start = pix
                    prev = pix
                ranges_list.append([start, prev + 1])
            if len(ranges_list) > 0:
                result_ranges = np.array(ranges_list, dtype=np.int64)
            else:
                result_ranges = np.empty((0, 2), dtype=np.int64)
        return result_ranges
    else:
        result_pixels = pixset_to_array(pixset, buff)
        if nest:
            # If original request was for NESTED, convert the result
            result_pixels = healpy.pixelfunc.ring2nest(nside, result_pixels)
        return result_pixels


def _boundaries_single(nside, pix, step=1, nest=False):
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    cdef vector[vec3] bounds
    if pix >= (12*nside*nside):
        raise ValueError('Pixel identifier is too large')
    hb.boundaries(pix, step, bounds)
    cdef size_t n = bounds.size()
    cdef np.ndarray[double, ndim = 2] out = np.empty((3, n), dtype=np.float64)
    for i in range(n):
        out[0,i] = bounds[i].x
        out[1,i] = bounds[i].y
        out[2,i] = bounds[i].z
    return out

def _boundaries_multiple(nside, pix, step=1, nest=False):
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    cdef size_t npix = len(pix)
    cdef size_t n = step * 4
    cdef size_t maxnpix = 12*nside*nside
    cdef np.ndarray[double, ndim = 3] out = np.empty((npix, 3, n), dtype=np.float64)
    cdef vector[vec3] bounds
 
    for j in range(npix):
        if pix[j] >= maxnpix:
            raise ValueError('Pixel identifier is too large')

        hb.boundaries(pix[j], step, bounds)
        for i in range(n):
            out[j, 0, i] = bounds[i].x
            out[j, 1, i] = bounds[i].y
            out[j, 2, i] = bounds[i].z
    return out

def boundaries(nside, pix, step=1, nest=False):
    """Returns an array containing vectors to the boundary of
    the nominated pixel.

    The returned array has shape (3, 4*step), the elements of
    which are the x,y,z positions on the unit sphere of the
    pixel boundary.  In order to get vector positions for just
    the corners, specify step=1.

    Parameters
    ----------
    nside : int
      The nside of the Healpix map.
    pix : int
      Pixel identifier
    step : int, optional
      Number of elements for each side of the pixel.
    nest : bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    Returns
    -------
    boundary : float, array
      x,y,z for positions on the boundary of the pixel

    Example
    -------
    # Get corners of HEALPixel 5 in nside=2, RINGS ordering.

    >>> import healpy as hp
    >>> import numpy as np
    >>> nside = 2
    >>> corners = hp.boundaries(nside, 5)

    # Now convert to phi,theta representation:
    >>> phi_theta = hp.vec2ang(np.transpose(corners))

    >>> corners = hp.boundaries(nside, np.array([5,5]))

    # doctest moved to test_query_disc.py
    """

    if not isnsideok(nside, nest):
        raise ValueError('Wrong nside value, must be a power of 2, less than 2**30')
    if np.isscalar(pix):
        if not np.can_cast(type(pix), int):
            raise ValueError('Array of pixels has to be integers')
        return _boundaries_single(nside, pix, step=step, nest=nest)
    else:
        pix = np.asarray(pix)
        if np.can_cast(pix.dtype, int):
            if pix.ndim!=1:
                raise ValueError('Array has to be one dimensional')
            return _boundaries_multiple(nside, pix.astype(int), step=step, nest=nest)
        else:
            raise ValueError('Array of pixels has to be integers')

# Try to implement pix2ang
### @cython.boundscheck(False)
### @cython.wraparound(False)
### def pix2ang(nside, ipix, nest = False):
###     # Check Nside value
###     if not isnsideok(nside, nest):
###         raise ValueError('Wrong nside value, must be a power of 2, less than 2**30')
###     # Create the Healpix_Base2 structure
###     cdef Healpix_Ordering_Scheme scheme
###     if nest:
###         scheme = NEST
###     else:
###         scheme = RING
###     cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
###     cdef pointing p
###     ipix_ = np.asarray(ipix, dtype = np.int64)
###     out_shape = (2,) + ipix_.shape
###     out = np.empty(out_shape, dtype = np.float64)
###     cdef np.ndarray[int64, ndim = 1] ipix_r = ipix_.reshape(-1)
###     cdef np.ndarray[double, ndim = 2] out_r = out.reshape(2, -1)
###     cdef int i, n
###     n = ipix_r.size
###     for i in range(n):
###         p = hb.pix2ang(ipix_r[i])
###         out_r[0, i] = p.theta
###         out_r[1, i] = p.phi
###     del out_r, ipix_r, ipix_
###     return out
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef rangeset_to_array(rangeset[int64] &pixset):
    """Convert rangeset to a 2D array of ranges.
    
    Returns a (num_ranges, 2) array where each row is [start, end) of a range.
    """
    cdef int64 i, n
    n = pixset.size()
    cdef np.ndarray[np.int64_t, ndim=2] ranges = np.empty((n, 2), dtype=np.int64)
    
    for i in range(n):
        ranges[i, 0] = pixset.ivbegin(i)
        ranges[i, 1] = pixset.ivend(i)
    
    return ranges

@cython.boundscheck(False)
@cython.wraparound(False)
cdef pixset_to_array(rangeset[int64] &pixset, buff=None):
    cdef int64 i, n
    n = pixset.nval()
    cdef np.ndarray[np.int64_t, ndim=1] ipix 
   		
    if buff is None :
       ipix = np.empty(n, dtype=np.int64)
    else :
       if n>=len(buff) :
           raise ValueError("Buffer too small to contain return value")
       ipix = buff[:n] 		
    
    cdef int64 a, b, ii, ip
    ii = 0
    n = pixset.size()
    for i in range(n):
        a = pixset.ivbegin(i)
        b = pixset.ivend(i)
        for ip in range(a, b):
            ipix[ii] = ip
            ii += 1
    return ipix




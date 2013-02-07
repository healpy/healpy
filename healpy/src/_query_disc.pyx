# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
cimport cython

from _common cimport int64, pointing, rangeset, vec3, Healpix_Ordering_Scheme, RING, NEST, SET_NSIDE, T_Healpix_Base

@cython.boundscheck(False)
@cython.wraparound(False)
def query_disc(nside, vec, radius, inclusive = False, fact = 4, nest = False):
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
      the resolution fact*nside. For NESTED ordering, fact must be a power of 2,
      else it can be any positive integer. Default: 4.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    Returns
    -------
    ipix : int, array
      The pixels which lie within the given disk.

    Note
    ----
    This method is more efficient in the RING scheme.
    For inclusive=True, the algorithm may return some pixels which don't overlap
    with the disk at all. The higher fact is chosen, the fewer false positives
    are returned, at the cost of increased run time.
    """
    # Check Nside value
    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef vec3 v = vec3(vec[0], vec[1], vec[2])
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    cdef rangeset[int64] pixset
    cdef int factor = fact
    if inclusive:
        factor = abs(fact)
        if nest and (factor == 0 or (factor & (factor - 1) != 0)):
            raise ValueError('fact must be a power of 2 when '
                             'nest is True (fact=%d)' % (fact))
        hb.query_disc_inclusive(pointing(v), radius, pixset, factor)
    else:
        hb.query_disc(pointing(v), radius, pixset)

    return pixset_to_array(pixset)


def query_polygon(nside, vertices, inclusive = False, fact = 4, nest = False):
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
      the resolution fact*nside. For NESTED ordering, fact must be a power of 2,
      else it can be any positive integer. Default: 4.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    
    Returns
    -------
    ipix : int, array
      The pixels which lie within the given polygon.

    Note
    ----
    This method is more efficient in the RING scheme.
    For inclusive=True, the algorithm may return some pixels which don't overlap
    with the disk at all. The higher fact is chosen, the fewer false positives
    are returned, at the cost of increased run time.
    """
    # Check Nside value
    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    # Create vector of vertices
    cdef vector[pointing] vert
    for v in vertices:
        vert.push_back(pointing(vec3(v[0], v[1], v[2])))
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
            raise ValueError('fact must be a power of 2 when '
                             'nest is True (fact=%d)' % (fact))
        hb.query_polygon_inclusive(vert, pixset, factor)
    else:
        hb.query_polygon(vert, pixset)

    return pixset_to_array(pixset)

def query_strip(nside, theta1, theta2, inclusive = False, nest = False):
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
      First colatitude
    theta2 : float
      Second colatitude
    inclusive ; bool
      If False, return the exact set of pixels whose pixels centers lie 
      within the region; if True, return all pixels that overlap with the
      region.
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    Returns
    -------
    ipix : int, array
      The pixels which lie within the given strip.
    """
    # Check Nside value
    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    # Create the Healpix_Base2 structure
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    # Call query_polygon
    cdef rangeset[int64] pixset
    hb.query_strip(theta1, theta2, inclusive, pixset)

    return pixset_to_array(pixset)


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
    >>> corners = hp.boundary(nside, 5)
    >>> print corners
    [[  2.44716655e-17   5.27046277e-01   3.60797400e-01   4.56398915e-17]
     [  3.99652627e-01   5.27046277e-01   8.71041977e-01   7.45355992e-01]
     [  9.16666667e-01   6.66666667e-01   3.33333333e-01   6.66666667e-01]]

    # Now convert to phi,theta representation:

    >>> hp.vec2ang(np.transpose(corners))
    (array([ 0.41113786,  0.84106867,  1.23095942,  0.84106867]), array([ 1.57079633,  0.78539816,  1.17809725,  1.57079633]))
    """


    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    cdef vector[vec3] bounds
    hb.boundaries(pix, step, bounds)
    cdef size_t n = bounds.size()
    cdef np.ndarray[double, ndim = 2] out = np.empty((3, n), dtype=np.float)
    for i in range(n):
        out[0,i] = bounds[i].x
        out[1,i] = bounds[i].y
        out[2,i] = bounds[i].z
    return out


# Try to implement pix2ang
### @cython.boundscheck(False)
### @cython.wraparound(False)
### def pix2ang(nside, ipix, nest = False):
###     # Check Nside value
###     if not isnsideok(nside):
###         raise ValueError('Wrong nside value, must be a power of 2')
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
    

cdef pixset_to_array(rangeset[int64] &pixset):
    cdef int64 i, n
    n = pixset.nval()
    cdef np.ndarray[np.int64_t, ndim=1] ipix = np.empty(n, dtype=np.int64)
    n = pixset.size()
    cdef int64 a, b, ii, ip
    ii = 0
    for i in range(n):
        a = pixset.ivbegin(i)
        b = pixset.ivend(i)
        for ip in range(a, b):
            ipix[ii] = ip
            ii += 1
    return ipix

cdef bool isnsideok(int nside):
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        return False
    else:
        return True



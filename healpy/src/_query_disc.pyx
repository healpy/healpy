# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
cimport cython

cdef extern from "healpix_base.h":
    ctypedef long int64

cdef extern from "rangeset.h":
    cdef cppclass interval[T]:
        T a()
        T b()
    cdef cppclass rangeset[T]:
        int64 size()
        int64 nval()
        interval[T] operator[](int i)

cdef extern from "vec3.h":
    cdef cppclass vec3:
        double x, y, z
        vec3()
        vec3(double xc, double yc, double zc)

cdef extern from "pointing.h":
    cdef cppclass pointing:
        pointing(vec3 inp)
        double theta
        double phi

cdef extern from "healpix_base.h":
    cdef enum Healpix_Ordering_Scheme:
        RING, NEST
    cdef cppclass nside_dummy:
        pass
    cdef nside_dummy SET_NSIDE
    cdef cppclass Healpix_Base2:
       Healpix_Base2(int64 nside, Healpix_Ordering_Scheme scheme,
                     nside_dummy)
       void query_disc (pointing ptg, double radius, bool inclusive,
                        rangeset[int64] pixset)
       void query_polygon(vector[pointing] vert, bool inclusive, 
                          rangeset[int64] pixset)
       void query_strip(double theta1, double theta2, bool inclusive,
                        rangeset[int64] pixset)

@cython.boundscheck(False)
@cython.wraparound(False)
def query_disc(nside, vec, radius, inclusive = False, nest = False):
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
    inclusive : bool
      If False, return the exact set of pixels whose pixel centers lie 
      within the disk; if True, return all pixels that overlap with the disk,
      and maybe a few more.

    Returns
    -------
    ipix : int, array
      The pixels which lie within the given disk.

    Note
    ----
    This method is more efficient in the RING scheme, but the algorithm
    used for inclusive==True returns fewer false positives in the NEST scheme.
    """
    # Check Nside value
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef vec3 v = vec3(vec[0], vec[1], vec[2])
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef Healpix_Base2 *hb = new Healpix_Base2(nside, scheme, SET_NSIDE)
    cdef rangeset[int64] pixset
    hb.query_disc(pointing(v), radius, inclusive, pixset)
    del hb

    return pixset_to_array(pixset)


def query_polygon(nside, vertices, inclusive = False, nest = False):
    """ Returns the pixels whose centers lie within the convex polygon 
    defined by the *vertices* array (if *inclusive* is False), or which 
    overlap with this polygon (if *inclusive* is True).

    Parameters
    ----------
    nside : int
      The nside of the Healpix map.
    vertices : float, array-like
      Vertex array containing the vertices of the polygon, shape (N, 3).
    inclusive : bool
      If False, return the exact set of pixels whose pixel centers lie
      within the polygon; if True, return all pixels that overlap with the
      polygon, and maybe a few more.
    
    Note
    ----
    This method is more efficient in the RING scheme, but the algorithm used
    for inclusive==True returns fewer false positives in the NEST scheme.
    """
    # Check Nside value
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
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
    cdef Healpix_Base2 *hb = new Healpix_Base2(nside, scheme, SET_NSIDE)
    # Call query_polygon
    cdef rangeset[int64] pixset
    hb.query_polygon(vert, inclusive, pixset)
    del hb

    return pixset_to_array(pixset)

def query_strip(nside, theta1, theta2, inclusive = False, nest = False):
    """Returns a range set of pixels whose centers lie within the colatitude
    range defined by *theta1* and *theta2* (if inclusive is False), or which 
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
    """
    # Check Nside value
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        raise ValueError('Wrong nside value, must be a power of 2')
    # Create the Healpix_Base2 structure
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef Healpix_Base2 *hb = new Healpix_Base2(nside, scheme, SET_NSIDE)
    # Call query_polygon
    cdef rangeset[int64] pixset
    hb.query_strip(theta1, theta2, inclusive, pixset)
    del hb

    return pixset_to_array(pixset)


cdef pixset_to_array(rangeset[int64] &pixset):
    cdef int64 i, n
    n = pixset.nval()
    cdef np.ndarray[np.int64_t, ndim=1] ipix = np.empty(n, dtype=np.int64)
    n = pixset.size()
    cdef int64 a, b, ii, ip
    ii = 0
    for i in range(n):
        a = pixset[i].a()
        b = pixset[i].b()
        for ip in range(a, b):
            ipix[ii] = ip
            ii += 1
    return ipix

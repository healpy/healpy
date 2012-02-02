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

cdef extern from "arr.h":
    cdef cppclass arr[T]:
        pass
    cdef cppclass fix_arr[T, int]:
        pass

cdef extern from "vec3.h":
    cdef cppclass vec3:
        double x, y, z
        vec3()
        vec3(double xc, double yc, double zc)

cdef extern from "pointing.h":
    cdef cppclass pointing:
        pointing()
        pointing(vec3 inp)
        double theta
        double phi

cdef extern from "healpix_base.h":
    ctypedef int val_4 "4"
    ctypedef int val_8 "8"
    cdef enum Healpix_Ordering_Scheme:
        RING, NEST
    cdef cppclass nside_dummy:
        pass
    cdef nside_dummy SET_NSIDE

    cdef cppclass T_Healpix_Base[I]:
       int nside2order(I nside)
       I npix2nside(I npix)
       T_Healpix_Base()
       T_Healpix_Base(int order, Healpix_Ordering_Scheme scheme)
       T_Healpix_Base(I nside, Healpix_Ordering_Scheme scheme,
                     nside_dummy)
       void query_disc (pointing ptg, double radius, bool inclusive,
                        rangeset[I]& pixset)
       void query_polygon(vector[pointing] vert, bool inclusive, 
                          rangeset[I]& pixset)
       void query_strip(double theta1, double theta2, bool inclusive,
                        rangeset[I]& pixset)
       void Set(int order, Healpix_Ordering_Scheme scheme)
       void SetNside(I nside, Healpix_Ordering_Scheme scheme)
       double ring2z(I ring)
       I pix2ring(I pix)
       I xyf2pix(int ix, int iy, int face_num)
       void pix2xyf(I pix, int &ix, int &iy, int &face_num)
       I nest2ring (I pix)
       I ring2nest (I pix)
       I nest2peano (I pix)
       I peano2nest (I pix)
       I zphi2pix (double z, double phi)
       I ang2pix (pointing &ang)
       I vec2pix (vec3 &vec)
       void pix2zphi (I pix, double &z, double &phi)
       pointing pix2ang (I pix)
       vec3 pix2vec (I pix)
       void get_ring_info (I ring, I &startpix, I &ringpix,
                           double &costheta, double &sintheta, bool &shifted)
       void get_ring_info2 (I ring, I &startpix, I &ringpix,
                            double &theta, bool &shifted)

       void get_ring_info_small (I ring, I &startpix, I &ringpix,
                                 bool &shifted)
       void neighbors (I pix, fix_arr[I,val_8] &result)
       void get_interpol (pointing &ptg, fix_arr[I,val_4] &pix,
                          fix_arr[double,val_4] &wgt)

       int Order()
       I Nside()
       I Npix()
       Healpix_Ordering_Scheme Scheme()
       bool conformable (T_Healpix_Base &other)
       void swap (T_Healpix_Base &other)
       double max_pixrad()
       double max_pixrad(I ring)
       arr[int] swap_cycles()


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
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

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
    cdef bool inc = inclusive
    hb.query_disc(pointing(v), radius, inc, pixset)

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
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    
    Returns
    -------
    ipix : int, array
      The pixels which lie within the given polygon.

    Note
    ----
    This method is more efficient in the RING scheme, but the algorithm used
    for inclusive==True returns fewer false positives in the NEST scheme.
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
    hb.query_polygon(vert, inclusive, pixset)

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
        a = pixset[i].a()
        b = pixset[i].b()
        for ip in range(a, b):
            ipix[ii] = ip
            ii += 1
    return ipix

cdef bool isnsideok(int nside):
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        return False
    else:
        return True



# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
from libcpp cimport bool
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

@cython.boundscheck(False)
@cython.wraparound(False)
def query_disc(nside, vec, radius, inclusive = False, nest = False):
    """query_disc
    """
    cdef vec3 *v
    v = new vec3(vec[0], vec[1], vec[2])
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    # Check Nside value
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef Healpix_Base2 *hb = new Healpix_Base2(nside, scheme, SET_NSIDE)
    cdef rangeset[int64] pixset
    hb.query_disc(pointing(v[0]), radius, inclusive, pixset)
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
    
    del hb
    del v
    return ipix

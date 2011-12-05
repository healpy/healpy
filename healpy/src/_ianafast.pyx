# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
#from libcpp cimport bool
#from libcpp.vector cimport vector
cimport cython

cdef typedef size_t tsize

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr(T *ptr, tsize sz)

cdef extern from "xcomplex.h":
    cdef cppclass xcomplex[T]:
       T re, im

cdef extern from "alm.h":
    cdef cppclass Alm[T]:
        Alm(int lmax_=0, int mmax_=0)

cdef extern from "healpix_map.h":
    cdef enum Healpix_Ordering_Scheme:
        RING, NEST
    cdef cppclass Healpix_Map[T]:
        Healpix_Map()
        Set(arr[T] &data, Healpix_Ordering_Scheme scheme)

cdef extern from "alm_healpix_tools.h":
    void map2alm_iter(Healpix_Map[T] &m,
                      Alm[xcomplex[T]] &alm,
                      int num_iter, arr[double] &weight)
    
def ianafast(m, lmax = -1, mmax = -1):
    pass


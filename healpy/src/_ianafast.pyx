# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
#from libcpp cimport bool
#from libcpp.vector cimport vector
cimport cython

cdef typedef size_t tsize

cdef extern from "healpix_tables.h":
    ctypedef enum Healpix_Ordering_Scheme {RING, NESTED}

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr(T *ptr, tsize sz)

cdef extern from "healpix_map.h":
    cdef cppclass Healpix_Map[T]:
        Healpix_Map()
        Set(arr[T] &data, Healpix_Ordering_Scheme scheme)

cdef extern from "alm_healpix_tools.h":
    pass

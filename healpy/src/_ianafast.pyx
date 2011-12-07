# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
#from libcpp cimport bool
#from libcpp.vector cimport vector
cimport cython
from healpy import npix2nside
ctypedef unsigned size_t
ctypedef size_t tsize

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr()
        arr(T *ptr, tsize sz)
        void allocAndFill (tsize sz, T &inival)

cdef extern from "xcomplex.h":
    cdef cppclass xcomplex[T]:
       T re, im

cdef extern from "alm.h":
    cdef cppclass Alm[T]:
        Alm(int lmax_=0, int mmax_=0)
        void Set (int lmax_, int mmax_)
        void Set (arr[T] &data, int lmax_, int mmax_)
        tsize Num_Alms (int l, int m)

cdef extern from "healpix_map.h":
    cdef enum Healpix_Ordering_Scheme:
        RING, NEST
    cdef cppclass Healpix_Map[T]:
        Healpix_Map()
        void Set(arr[T] &data, Healpix_Ordering_Scheme scheme)
        T average()
        void Add(T x)

cdef extern from "alm_healpix_tools.h":
    cdef void map2alm_iter(Healpix_Map[double] &m,
                           Alm[xcomplex[double]] &alm,
                           int num_iter, arr[double] &weight)


cdef extern from "hack.h":
    cdef xcomplex[double]* cast_to_ptr_xcomplex_d(char *)


cdef Healpix_Map[double] as_map(np.ndarray[double, ndim=1] m):
    cdef arr[double] m_arr = arr[double](<double*>m.data, m.size)
    cdef Healpix_Map[double] m_hpx
    m_hpx.Set(m_arr, RING)
    return m_hpx

cdef Alm[xcomplex[double]] as_alm(np.ndarray[double, ndim=1] a, int lmax, int mmax):
    cdef Alm[xcomplex[double]] alm_hpx
    cdef tsize n_alm = Num_Alms(lmax, mmax)
    if a.size != n_alm:
        raise ValueError("inconsistent array size, lmax, mmax")
    cdef arr[xcomplex[double]] alm_arr = arr[xcomplex[double]](cast_to_ptr_xcomplex_d(a.data), a.size)
    alm_hpx.Set(alm_arr, lmax, mmax)
    return alm_hpx

cdef Num_Alms(int l, int m):
    if not m <= l:
        raise ValueError("mmax must be <= lmax")
    return ((m+1)*(m+2))/2 + (m+1)*(l-m)


def ianafast(m, lmax = -1, mmax = -1, niter = 3, regression = True):

    # Get the map as a contiguous ndarray object if it isn't
    cdef np.ndarray[np.float64_t, ndim=1] m_
    m_ = np.ascontiguousarray(m, dtype = np.float64)

    # Adjust lmax and mmax
    cdef int lmax_ = lmax, mmax_ = mmax, nside, npix
    npix = m_.size
    nside = npix2nside(npix)
    if lmax_ < 0:
        lmax_ = 3 * nside - 1
    if mmax_ < 0 or mmax_ > lmax_:
        mmax_ = lmax_

    # Create an Healpix_Map object pointing to data from ndarray map
    m_hpx = as_map(m_)
    cdef double avg = 0.0

    if regression:
        avg = m_hpx.average()
        m_hpx.Add(-avg)

    # Create an ndarray object that will contain the alm for output (to be returned)
    cdef np.ndarray alm = np.empty(Num_Alms(lmax_, mmax_), dtype = np.complex128)
    alm_hpx = as_alm(alm, lmax_, mmax_)

    # ring weights
    cdef arr[double] w_arr
    w_arr.allocAndFill(2 * nside, 1.)
    
    map2alm_iter(m_hpx, alm_hpx, niter, w_arr)
 
    if regression:
        m_hpx.Add(avg)
        alm[0] += avg * np.sqrt(4 * np.pi)
   
    return alm



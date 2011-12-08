# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
#from libcpp cimport bool
#from libcpp.vector cimport vector
from libc.math cimport sqrt
cimport libc
from healpy import npix2nside
from healpy.pixelfunc import maptype
ctypedef unsigned size_t
ctypedef size_t tsize

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr()
        arr(T *ptr, tsize sz)
        void allocAndFill (tsize sz, T &inival)
        tsize size()

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
    cdef void map2alm_pol_iter(Healpix_Map[double] &mapT,
                               Healpix_Map[double] &mapQ,
                               Healpix_Map[double] &mapU,
                               Alm[xcomplex[double]] &almT,
                               Alm[xcomplex[double]] &almG,
                               Alm[xcomplex[double]] &almC,
                               int num_iter,
                               arr[double] &weight)

cdef extern from "hack.h":
    cdef xcomplex[double]* cast_to_ptr_xcomplex_d(char *)

cdef Num_Alms(int l, int m):
    if not m <= l:
        raise ValueError("mmax must be <= lmax")
    return ((m+1)*(m+2))/2 + (m+1)*(l-m)

cdef class WrapMap(object):
    """This class provides a wrapper to a ndarray so it can be sent as Map to healpix_cxx functions.
    """
    cdef Healpix_Map[double] * h
    cdef arr[double] * a
    cdef object m
    cdef int is_init

    def __init__(self, np.ndarray[double] m):
        if self.is_init == 1:
            raise Exception('Already init...')
        self.is_init = 1
        self.m = np.ascontiguousarray(m)
        self.a = new arr[double](<double*>(<np.ndarray>self.m).data, (<np.ndarray>self.m).size)
        self.h = new Healpix_Map[double]()
        self.h.Set(self.a[0], RING)

    def __dealloc__(self):
        if self.is_init == 1:
            #print "deallocating map wrapper..."
            del self.a, self.h

cdef class WrapAlm(object):
    """This class provides a wrapper to a ndarray so it can be sent as Alm to healpix_cxx functions.
    """
    cdef Alm[xcomplex[double]] * h
    cdef arr[xcomplex[double]] * a
    cdef object m
    cdef int is_init

    def __init__(self, np.ndarray[np.complex128_t] m, int lmax, int mmax):
        if self.is_init == 1:
            raise Exception('Already init...')
        self.is_init = 1
        self.m = np.ascontiguousarray(m)
        self.h = new Alm[xcomplex[double]]()
        self.a = new arr[xcomplex[double]](cast_to_ptr_xcomplex_d((<np.ndarray>self.m).data),
                                           (<np.ndarray>self.m).size)
        self.h.Set(self.a[0], lmax, mmax)

    def __dealloc__(self):
        if self.is_init == 1:
            #print "deallocating alm wrapper..."
            del self.a, self.h

def ianafast(m, lmax = -1, mmax = -1, niter = 3, regression = True):

    # Check if the input map is polarized or not
    info = maptype(m)
    if info == 0:
        polarization = False
        mmi = m
    elif info == 1:
        polarization = False
        mmi = m[0]
    elif info == 3:
        polarization = True
        mmi = m[0]
        mmq = m[1]
        mmu = m[2]
    else:
        raise ValueError("Wrong input map (must be a valid healpix map or a sequence of 1 or 3 maps)")
        
    # Get the map as a contiguous ndarray object if it isn't
    cdef np.ndarray[np.float64_t, ndim=1] mi, mq, mu
    mi = np.ascontiguousarray(mmi, dtype = np.float64)
    if polarization:
        mq = np.ascontiguousarray(mmq, dtype = np.float64)
        mu = np.ascontiguousarray(mmu, dtype = np.float64)
    
    # Adjust lmax and mmax
    cdef int lmax_ = lmax, mmax_ = mmax, nside, npix
    npix = mi.size
    nside = npix2nside(npix)
    if lmax_ < 0:
        lmax_ = 3 * nside - 1
    if mmax_ < 0 or mmax_ > lmax_:
        mmax_ = lmax_
    
    # Wrap the map into an Healpix_Map
    MI = WrapMap(mi)
    if polarization:
        MQ = WrapMap(mq)
        MU = WrapMap(mu)

    # if regression is True, remove average of the intensity map before computing alm
    cdef double avg = 0.0
    if regression:
        avg = MI.h.average()
        MI.h.Add(-avg)

    # Create an ndarray object that will contain the alm for output (to be returned)
    cdef np.ndarray almI, almQ, almC
    n_alm = Num_Alms(lmax_, mmax_)
    almI = np.empty(n_alm, dtype = np.complex128)
    if polarization:
        almG = np.empty(n_alm, dtype = np.complex128)
        almC = np.empty(n_alm, dtype = np.complex128)
    
    # Wrap it into an healpix Alm object
    AI = WrapAlm(almI, lmax_, mmax_)
    if polarization:
        AG = WrapAlm(almG, lmax_, mmax_)
        AC = WrapAlm(almC, lmax_, mmax_)

    # ring weights
    cdef arr[double] * w_arr = new arr[double]()
    w_arr.allocAndFill(2 * nside, 1.)
    
    if polarization:
        map2alm_pol_iter(MI.h[0], MQ.h[0], MU.h[0], AI.h[0], AG.h[0], AC.h[0],
                         niter, w_arr[0])
    else:
        map2alm_iter(MI.h[0], AI.h[0], niter, w_arr[0])
    
    if regression:
        MI.h.Add(avg)
        almI[0] += avg * sqrt(4 * np.pi)

    del w_arr
    if polarization:
        return almI, almG, almC
    else:
        return almI

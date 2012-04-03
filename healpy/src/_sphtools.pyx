# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
#from libcpp cimport bool
#from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.math cimport sqrt, floor, fabs
cimport libc
from healpy import npix2nside
from healpy.pixelfunc import maptype
ctypedef unsigned size_t
ctypedef size_t tsize
import os
import cython

cdef double UNSEEN = -1.6375e30
cdef double rtol_UNSEEN = 1.e-7 * 1.6375e30

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr()
        arr(T *ptr, tsize sz)
        void allocAndFill (tsize sz, T &inival)
        tsize size()
        T &operator[] (tsize n)

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

cdef extern from "healpix_data_io.h":
    cdef void read_weight_ring (string &dir, int nside, arr[double] &weight)

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

DATAPATH = None

def get_datapath():
    global DATAPATH
    if DATAPATH is None:
        DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    return DATAPATH

def map2alm(m, lmax = None, mmax = None, niter = 3, use_weights = False, 
            regression = True, datapath = None):
    """Computes the alm of an Healpix map.

    Parameters
    ----------
    m : array-like, shape (Npix,) or (3, Npix)
      The input map or a list of 3 input maps (polariztion).
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    iter : int, scalar, optional
      Number of iteration (default: 1)
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    regression: bool, scalar, optional
      If True, subtract map average before computing alm. Default: True.
    
    Returns
    -------
    alm : array or tuple of arrays
      alm or a tuple of 3 alm (almT, almE, almB) if polarized input.
    """
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
        raise ValueError("Wrong input map (must be a valid healpix map "
                         "or a sequence of 1 or 3 maps)")
    
    # Get the map as a contiguous ndarray object if it isn't
    cdef np.ndarray[np.float64_t, ndim=1] mi, mq, mu
    mi = np.ascontiguousarray(mmi, dtype = np.float64)
    # create UNSEEN mask for I map
    mask_mi = False if count_bad(mi) == 0 else mkmask(mi)
    # same for polarization maps if needed
    if polarization:
        mq = np.ascontiguousarray(mmq, dtype = np.float64)
        mask_mq = False if count_bad(mq) == 0 else mkmask(mq)
        mu = np.ascontiguousarray(mmu, dtype = np.float64)
        mask_mu = False if count_bad(mu) == 0 else mkmask(mu)

    # replace UNSEEN pixels with zeros
    if mask_mi is not False:
        mi[mask_mi] = 0.0
    if polarization:
        if mask_mq is not False:
            mq[mask_mq] = 0.0
        if mask_mu is not False:
            mu[mask_mu] = 0.0

    # Adjust lmax and mmax
    cdef int lmax_, mmax_, nside, npix
    npix = mi.size
    nside = npix2nside(npix)
    if lmax is None:
        lmax_ = 3 * nside - 1
    else:
        lmax_ = lmax
    if mmax is None:
        mmax_ = lmax_
    else:
        mmax_ = mmax

    # Check all maps have same npix
    if polarization:
        if mq.size != npix or mu.size != npix:
            raise ValueError("Input maps must have same size")
    
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
    cdef int i
    cdef char *c_datapath
    if use_weights:
        if datapath is None:
            datapath = get_datapath()
        c_datapath = datapath
        weightfile = 'weight_ring_n%05d.fits' % (nside)
        if not os.path.isfile(os.path.join(datapath, weightfile)):
            raise IOError('Weight file not found in %s' % (datapath))
        read_weight_ring(string(c_datapath), nside, w_arr[0])
        for i in range(w_arr.size()):
            w_arr[0][i] += 1
    else:
        w_arr.allocAndFill(2 * nside, 1.)
    
    if polarization:
        map2alm_pol_iter(MI.h[0], MQ.h[0], MU.h[0], AI.h[0], AG.h[0], AC.h[0],
                         niter, w_arr[0])
    else:
        map2alm_iter(MI.h[0], AI.h[0], niter, w_arr[0])
    
    # restore input map with UNSEEN pixels
    if mask_mi is not False:
        mi[mask_mi] = UNSEEN
    if polarization:
        if mask_mq is not False:
            mq[mask_mq] = UNSEEN
        if mask_mu is not False:
            mu[mask_mu] = UNSEEN
    
    if regression:
        MI.h.Add(avg)
        almI[0] += avg * sqrt(4 * np.pi)

    del w_arr
    if polarization:
        return almI, almG, almC
    else:
        return almI


def alm2cl(alms, alms2 = None, lmax = None, mmax = None, lmax_out = None):
    """Computes (cross-)spectra from alm(s). If alm2 is given, cross-spectra between
    alm and alm2 are computed. If alm (and alm2 if provided) contains n alm,
    then n(n+1)/2 auto and cross-spectra are returned.

    Parameters
    ----------
    alms : complex, array or sequence of arrays
      The alm from which to compute the power spectrum. If n>=2 arrays are given,
      computes both auto- and cross-spectra.
    alms2 : complex, array or sequence of 3 arrays, optional
      If provided, computes cross-spectra between alm and alm2.
      Default: alm2=alm, so auto-spectra are computed.
    lmax : None or int, optional
      The maximum l of the input alm. Default: computed from size of alm
      and mmax_in
    mmax : None or int, optional
      The maximum m of the input alm. Default: assume mmax_in = lmax_in
    lmax_out : None or int, optional
      The maximum l of the returned spectra. By default: the lmax of the given
      alm(s).

    Returns
    -------
    cl : array or tuple of n(n+1)/2 arrays
      the spectrum <*alm* x *alm2*> if *alm* (and *alm2*) is one alm, or 
      the auto- and cross-spectra <*alm*[i] x *alm2*[j]> if alm (and alm2)
      contains more than one spectra.
      If more than one spectrum is returned, they are ordered by diagonal.
      For example, if *alm* is almT, almE, almB, then the returned spectra are:
      TT, EE, BB, TE, EB, TB.
    """
    #############################
    # Check alm and number of spectra
    #
    cdef int Nspec, Nspec2
    if not hasattr(alms, '__len__'):
        raise ValueError('alms must be an array or a sequence of arrays')
    if not hasattr(alms[0], '__len__'):
        alms_lonely = True
        alms = [alms]
    else:
        alms_lonely = False

    Nspec = len(alms)

    if alms2 is None:
        alms2 = alms
    
    if not hasattr(alms2, '__len__'):
        raise ValueError('alms2 must be an array or a sequence of arrays')
    if not hasattr(alms2[0], '__len__'):
        alms2 = [alms2]
    Nspec2 = len(alms2)
    
    if Nspec != Nspec2:
        raise ValueError('alms and alms2 must have same number of spectra')
    
    ##############################################
    # Check sizes of alm's and lmax/mmax/lmax_out
    #
    cdef int almsize
    almsize = alms[0].size
    for i in xrange(Nspec):
        if alms[i].size != almsize or alms2[i].size != almsize:
            raise ValueError('all alms must have same size')

    if lmax is None:
        if mmax is None:
            lmax = alm_getlmax(almsize)
            mmax = lmax
        else:
            lmax = alm_getlmax2(almsize, mmax)

    if mmax is None:
        mmax = lmax

    if lmax_out is None:
        lmax_out = lmax


    #######################
    # Computing the spectra
    #
    cdef int j, l, m, limit
    cdef int lmax_ = lmax, mmax_ = mmax
    cdef int lmax_out_ = lmax_out

    cdef np.ndarray[double, ndim=1] powspec_
    cdef np.ndarray[np.complex128_t, ndim=1] alm1_
    cdef np.ndarray[np.complex128_t, ndim=1] alm2_

    spectra = []
    for n in xrange(Nspec): # diagonal rank
        for m in xrange(0, Nspec - n): # position in the diagonal
            powspec_ = np.zeros(lmax + 1)
            alm1_ = alms[m]
            alm2_ = alms2[m + n]
            # compute cross-spectrum alm1[n] x alm2[n+m]
            # and place result in result list
            for l in range(lmax_ + 1):
                j = alm_getidx(lmax_, l, 0)
                powspec_[l] = alm1_[j].real * alm2_[j].real
                limit = l if l <= mmax else mmax
                for m in range(1, limit + 1):
                    j = alm_getidx(lmax_, l, m)
                    powspec_[l] += 2 * (alm1_[j].real * alm2_[j].real +
                                        alm1_[j].imag * alm2_[j].imag)
                powspec_[l] /= (2 * l + 1)
            spectra.append(powspec_)
    
    # if only one alm was given, returns only cl and not a list with one cl
    if alms_lonely:
        spectra = spectra[0]

    return spectra


@cython.wraparound(False)
@cython.boundscheck(False)
def almxfl(alm, fl, mmax = None, inplace = False):
    """Multiply an a_lm by a vector b_l.

    Parameters
    ----------
    alm : array, double
      The array representing the spherical harmonics coefficients
    fl : array, double
      The array giving the factor f_l by which to multiply a_lm
    mmax : None or int, optional
      The maximum m of the input alm
    inplace : bool, optional
      If True, performs the computation in-place if possible (input alm
      is modified if it is a 1d-array of type float64). Otherwise,
      a copy of alm is done.

    Returns
    -------
    alm : array, double
      The result of a_lm * f_l. If *inplace* is True, returns the input
      alm modified
    """
    cdef np.ndarray[np.complex128_t, ndim=1] alm_
    cdef np.ndarray[np.complex128_t, ndim=1] fl_

    if inplace:
        alm_ = np.ascontiguousarray(alm, dtype = np.complex128)
    else:
        alm_ = np.array(alm, dtype = np.complex128, copy = True)
    
    fl_ = np.ascontiguousarray(fl, dtype = np.complex128)

    cdef int lmax_, mmax_
    cdef int l, m
    if mmax is None:
        lmax_ = alm_getlmax(alm_.size)
        mmax_ = lmax_
    else:
        lmax_ = alm_getlmax2(alm_.size, mmax)
        mmax_ = mmax
    
    cdef np.complex128_t f
    cdef int maxm, i
    cdef int flsize = fl_.size
    for l in xrange(lmax_ + 1):
        f = fl_[l] if l < flsize else 0.
        maxm = l if l <= mmax_ else mmax_
        for m in xrange(maxm + 1):
            i = alm_getidx(lmax_, l, m)
            alm_[i] *= f

    return alm_


@cython.cdivision(True)
cdef inline int alm_getidx(int lmax, int l, int m):
    return m*(2*lmax+1-m)/2+l


@cython.cdivision(True)
cdef inline int alm_getlmax(int s):
    cdef double x
    x=(-3+np.sqrt(1+8*s))/2
    if x != floor(x):
        return -1
    else:
        return <int>floor(x)

@cython.cdivision(True)
cdef inline int alm_getlmax2(int s, int mmax):
    cdef double x
    x = (2 * s + mmax ** 2 - mmax - 2.) / (2 * mmax + 2.)
    if x != floor(x):
        return -1
    else:
        return <int>floor(x)

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef mkmask(np.ndarray[double, ndim=1] m):
    cdef int nbad
    cdef int size = m.size
    cdef int i
    # first, count number of bad pixels, to see if allocating a mask is needed
    nbad = count_bad(m)
    cdef np.ndarray[np.int8_t, ndim=1] mask
    #cdef np.ndarray[double, ndim=1] m_
    if nbad == 0:
        return False
    else:
        mask = np.zeros(size, dtype = np.int8)
        #m_ = m
        for i in range(size):
            if fabs(m[i] - UNSEEN) < rtol_UNSEEN:
                mask[i] = 1
    mask.dtype = bool
    return mask
            
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef int count_bad(np.ndarray[double, ndim=1] m):
    cdef int i
    cdef int nbad = 0
    cdef size = m.size
    for i in xrange(m.size):
        if fabs(m[i] - UNSEEN) < rtol_UNSEEN:
            nbad += 1
    return nbad


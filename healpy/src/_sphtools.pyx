# coding: utf-8
# Wrapper around the query_disc method of Healpix_base class

import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libc.math cimport sqrt, floor, fabs
cimport libc
from healpy import npix2nside, nside2npix
from healpy.pixelfunc import maptype
import os
import cython
from libcpp cimport bool as cbool

from _common cimport tsize, arr, xcomplex, Healpix_Ordering_Scheme, RING, NEST, Healpix_Map, Alm, ndarray2map, ndarray2alm

cdef double UNSEEN = -1.6375e30
cdef double rtol_UNSEEN = 1.e-7 * 1.6375e30

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
    cdef void map2alm_spin(    Healpix_Map[double] &map1, 
                               Healpix_Map[double] &map2,
                               Alm[xcomplex[double]] &alm1, 
                               Alm[xcomplex[double]] &alm2, 
                               int spin, 
                               arr[double] &weight, 
                               cbool add_alm)
    cdef void alm2map_spin(    Alm[xcomplex[double]] &alm1, 
                               Alm[xcomplex[double]] &alm2, 
                               Healpix_Map[double] &map1, 
                               Healpix_Map[double] &map2,
                               int spin) 

cdef extern from "alm_powspec_tools.h":
    cdef void c_rotate_alm "rotate_alm" (Alm[xcomplex[double]] &alm, double psi, double theta, double phi)
    cdef void c_rotate_alm "rotate_alm" (Alm[xcomplex[double]] &ai, Alm[xcomplex[double]] &ag, Alm[xcomplex[double]] &ac, double psi, double theta, double phi)

cdef extern from "healpix_data_io.h":
    cdef void read_weight_ring (string &dir, int nside, arr[double] &weight)

DATAPATH = None

def get_datapath():
    global DATAPATH
    if DATAPATH is None:
        DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    return DATAPATH

def map2alm_spin_healpy(maps, spin, lmax = None, mmax = None):
    """Computes the spinned alm of a 2 Healpix maps.

    Parameters
    ----------
    m : list of 2 arrays
        list of 2 input maps as numpy arrays
    spin : int
        spin of the alms (either 1, 2 or 3)
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    
    Returns
    -------
    alms : list of 2 arrays
      list of 2 alms
    """
    maps_c = [np.ascontiguousarray(m, dtype=np.float64) for m in maps]

    # create UNSEEN mask for map
    masks = [False if count_bad(m) == 0 else mkmask(m) for m in maps_c]

    # Adjust lmax and mmax
    cdef int lmax_, mmax_, nside, npix
    npix = maps_c[0].size
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
    if maps_c[1].size != npix:
        raise ValueError("Input maps must have same size")
    
    # View the ndarray as a Healpix_Map
    M1 = ndarray2map(maps_c[0], RING)
    M2 = ndarray2map(maps_c[1], RING)

    # replace UNSEEN pixels with zeros
    for m, mask in zip(maps_c, masks):
        if mask:
            m[mask] = 0.0

    # Create an ndarray object that will contain the alm for output (to be returned)
    n_alm = alm_getn(lmax_, mmax_)
    alms = [np.empty(n_alm, dtype=np.complex128) for m in maps]

    # View the ndarray as an Alm
    # Alms = [ndarray2alm(alm, lmax_, mmax_) for alm in alms]
    A1 = ndarray2alm(alms[0], lmax_, mmax_) 
    A2 = ndarray2alm(alms[1], lmax_, mmax_) 
    
    # ring weights
    cdef arr[double] * w_arr = new arr[double]()
    cdef int i
    cdef char *c_datapath
    w_arr.allocAndFill(2 * nside, 1.)
    
    map2alm_spin(M1[0], M2[0], A1[0], A2[0], spin, w_arr[0], False)
    
    # restore input map with UNSEEN pixels
    for m, mask in zip(maps_c, masks):
        if mask:
            m[mask] = UNSEEN

    del w_arr
    del M1, M2, A1, A2
    return alms

def alm2map_spin_healpy(alms, nside, spin, lmax, mmax=None):
    """Computes maps from a set of 2 spinned alm

    Parameters
    ----------
    alms : list of 2 arrays
      list of 2 alms
    nside : int
        requested nside of the output map 
    spin : int
        spin of the alms (either 1, 2 or 3)
    lmax : int, scalar
      Maximum l of the power spectrum.
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    
    Returns
    -------
    m : list of 2 arrays
        list of 2 out maps in RING scheme as numpy arrays
    """
    alms_c = [np.ascontiguousarray(alm, dtype=np.complex128) for alm in alms]

    npix = nside2npix(nside)
    maps = [np.zeros(npix, dtype=np.float64) for alm in alms]

    # View the ndarray as a Healpix_Map
    M1 = ndarray2map(maps[0], RING)
    M2 = ndarray2map(maps[1], RING)

    if not mmax:
        mmax = lmax

    # View the ndarray as an Alm
    A1 = ndarray2alm(alms_c[0], lmax, mmax) 
    A2 = ndarray2alm(alms_c[1], lmax, mmax) 
    
    alm2map_spin(A1[0], A2[0], M1[0], M2[0], spin)
    
    del M1, M2, A1, A2
    return maps

def map2alm(m, lmax = None, mmax = None, niter = 3, use_weights = False, 
            datapath = None):
    """Computes the alm of a Healpix map.

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
    
    Returns
    -------
    alm : array or tuple of arrays
      alm or a tuple of 3 alm (almT, almE, almB) if polarized input.
    """
    # Check if the input map is polarized or not
    info = maptype(m)
    if info == 0:
        polarization = False
        mi = np.ascontiguousarray(m, dtype=np.float64)
    elif info == 1:
        polarization = False
        mi = np.ascontiguousarray(m[0], dtype=np.float64)
    elif info == 3:
        polarization = True
        mi = np.ascontiguousarray(m[0], dtype=np.float64)
        mq = np.ascontiguousarray(m[1], dtype=np.float64)
        mu = np.ascontiguousarray(m[2], dtype=np.float64)
    else:
        raise ValueError("Wrong input map (must be a valid healpix map "
                         "or a sequence of 1 or 3 maps)")

    # create UNSEEN mask for I map
    mask_mi = False if count_bad(mi) == 0 else mkmask(mi)
    # same for polarization maps if needed
    if polarization:
        mask_mq = False if count_bad(mq) == 0 else mkmask(mq)
        mask_mu = False if count_bad(mu) == 0 else mkmask(mu)

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
    
    # View the ndarray as a Healpix_Map
    MI = ndarray2map(mi, RING)
    if polarization:
        MQ = ndarray2map(mq, RING)
        MU = ndarray2map(mu, RING)

    # replace UNSEEN pixels with zeros
    if mask_mi is not False:
        mi[mask_mi] = 0.0
    if polarization:
        if mask_mq is not False:
            mq[mask_mq] = 0.0
        if mask_mu is not False:
            mu[mask_mu] = 0.0


    # Create an ndarray object that will contain the alm for output (to be returned)
    n_alm = alm_getn(lmax_, mmax_)
    almI = np.empty(n_alm, dtype=np.complex128)
    if polarization:
        almG = np.empty(n_alm, dtype=np.complex128)
        almC = np.empty(n_alm, dtype=np.complex128)

    # View the ndarray as an Alm
    AI = ndarray2alm(almI, lmax_, mmax_)
    if polarization:
        AG = ndarray2alm(almG, lmax_, mmax_)
        AC = ndarray2alm(almC, lmax_, mmax_)
    
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
        map2alm_pol_iter(MI[0], MQ[0], MU[0], AI[0], AG[0], AC[0], niter, w_arr[0])
    else:
        map2alm_iter(MI[0], AI[0], niter, w_arr[0])
    
    # restore input map with UNSEEN pixels
    if mask_mi is not False:
        mi[mask_mi] = UNSEEN
    if polarization:
        if mask_mq is not False:
            mq[mask_mq] = UNSEEN
        if mask_mu is not False:
            mu[mask_mu] = UNSEEN
    
    del w_arr
    if polarization:
        del MI, MQ, MU, AI, AG, AC
        return almI, almG, almC
    else:
        del MI, AI
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

    lmax, mmax = alm_getlmmax(alms[0], lmax, mmax)

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
    lmax_, mmax_ = alm_getlmmax(alm_, None, mmax)
    
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


def rotate_alm(alm not None, double psi, double theta, double phi, lmax=None,
               mmax=None):
    """
    This routine transforms the scalar (and tensor) a_lm coefficients
    to emulate the effect of an arbitrary rotation of the underlying
    map. The rotation is done directly on the a_lm using the Wigner
    rotation matrices, computed by recursion. To rotate the a_lm for
    l ≤ l_max the number of operations scales like l_max^3.

    Parameters
    ----------
    alm : array-like of shape (n,) or (k,n), or list of arrays
        Complex a_lm values before and after rotation of the coordinate system.
    psi : float
        First rotation: angle ψ about the z-axis. All angles are in radians
        and should lie in [-2pi,2pi], the rotations are active and the
        referential system is assumed to be right handed. The routine
        coordsys2euler zyz can be used to generate the Euler angles ψ, θ, φ
        for rotation between standard astronomical coordinate systems.
    theta : float
        Second rotation: angle θ about the original (unrotated) y-axis
    phi : float.
        Third rotation: angle φ about the original (unrotated) z-axis.
    lmax : int
        Maximum multipole order l of the data set.
    mmax : int
        Maximum degree m of data set.

    """
    if isinstance(alm, np.ndarray) and alm.ndim == 1:
        alm = [alm]

    if not isinstance(alm, (list, tuple, np.ndarray)) or len(alm) == 0:
        raise ValueError('Invalid input.')

    # C++ rotate_alm only handles 1 or 3 maps. The function handling 3 maps
    # is faster than running 3 times the 1-map function, but gives identical
    # results.
    if len(alm) not in (1, 3):
        for a in alm:
            rotate_alm(a, psi, theta, phi)
        return

    lmax, mmax = alm_getlmmax(alm[0], lmax, mmax)
    ai = np.ascontiguousarray(alm[0], dtype=np.complex128)
    AI = ndarray2alm(ai, lmax, mmax)
    if len(alm) == 1:
        c_rotate_alm(AI[0], psi, theta, phi)
        del AI
    else:
        ag = np.ascontiguousarray(alm[1], dtype=np.complex128)
        ac = np.ascontiguousarray(alm[2], dtype=np.complex128)
        AG = ndarray2alm(ag, lmax, mmax)
        AC = ndarray2alm(ac, lmax, mmax)
        c_rotate_alm(AI[0], AG[0], AC[0], psi, theta, phi)
        del AI, AG, AC


cdef int alm_getn(int l, int m):
    if not m <= l:
        raise ValueError("mmax must be <= lmax")
    return ((m+1)*(m+2))/2 + (m+1)*(l-m)


def alm_getlmmax(a, lmax, mmax):
    if lmax is None:
        if mmax is None:
            lmax = alm_getlmax(a.size)
            mmax = lmax
        else:
            lmax = alm_getlmax2(a.size, mmax)
    elif mmax is None:
        mmax = lmax
    return lmax, mmax


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


@cython.cdivision(True)
cdef inline int alm_getidx(int lmax, int l, int m):
    return m*(2*lmax+1-m)/2+l


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

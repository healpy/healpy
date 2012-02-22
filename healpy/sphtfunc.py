# 
#  This file is part of Healpy.
# 
#  Healpy is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  Healpy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Healpy; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
#  For more information about Healpy, see http://code.google.com/p/healpy
# 
import warnings
import exceptions
import numpy as np
pi = np.pi

import _healpy_sph_transform_lib as sphtlib
import _healpy_fitsio_lib as hfitslib
import healpy._ianafast as _ianafast

import os.path
import pixelfunc

from pixelfunc import mask_bad, maptype, UNSEEN

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Spherical harmonics transformation
def anafast(map1, map2 = None, lmax = None, mmax = None, 
            iter = 3, alm = False, use_weights = False, regression = True, 
            datapath = None, nspec = None):
    """Computes the power spectrum of an Healpix map, or the cross-spectrum
    between two maps if *map2* is given.

    Parameters
    ----------
    map1 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps
    map2 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps
    nspec : None or int, optional
      The number of spectra to return. If None, returns all, otherwise
      returns cls[:nspec]
    lmax : int, scalar, optional
      Maximum l of the power spectrum (default: 3*nside-1)
    mmax : int, scalar, optional
      Maximum m of the alm (default: lmax)
    iter : int, scalar, optional
      Number of iteration (default: 3)
    alm : bool, scalar, optional
      If True, returns both cl and alm, otherwise only cl is returned
    regression : bool, scalar, optional
      If True, map average is removed before computing alm. Default: True.
    datapath : None or str, optional
      If given, the directory where to find the weights data.

    Returns
    -------
    res : array or sequence of arrays
      If *alm* is False, returns cl or a list of cl's (TT, EE, BB, TE, EB, TB for
      polarized input map)
      Otherwise, returns a tuple (cl, alm), where cl is as above and
      alm is the spherical harmonic transform or a list of almT, almE, almB
      for polarized input
    """
    alms1 = map2alm(map1, lmax = lmax, mmax = mmax, niter = iter, 
                    use_weights = use_weights, regression = regression, 
                    datapath = datapath)
    if map2 is not None:
        alms2 = map2alm(map2, lmax = lmax, mmax = mmax, niter = iter, 
                    use_weights = use_weights, regression = regression, 
                    datapath = datapath)
    else:
        alms2 = None
    
    cls = alm2cl(alms1, alm2 = alms2, lmax = lmax, mmax = mmax,
                 lmax_out = lmax, nspec = nspec)

    if alm:
        if map2 is not None:
            return (cls, alms1, alms2)
        else:
            return (cls, alms1)
    else:
        return cls

def map2alm(m, lmax = None, mmax = None, iter = 3, use_weights = False, 
            regression = True, datapath = None):
    """Computes the alm of an Healpix map.

    Parameters
    ----------
    m : array-like, shape (Npix,) or (3, Npix)
      The input map or a list of 3 input maps (polarization).
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    iter : int, scalar, optional
      Number of iteration (default: 3)
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    regression: bool, scalar, optional
      If True, subtract map average before computing alm. Default: True.
    
    Returns
    -------
    alm : array or tuple of array
      alm or a tuple of 3 alm (almT, almE, almB) if polarized input.
    """
    alm = _ianafast.map2alm(m, niter = iter, regression = regression, 
                            datapath = datapath, use_weights = use_weights,
                            lmax = lmax, mmax = mmax)
    return alm

def alm2map(alm, nside, lmax = None, mmax = None, pixwin = False,
            fwhm = 0.0, sigma = None):
    """Computes an Healpix map given the alm.

    The alm are given as a complex array. You can specify lmax
    and mmax, or they will be computed from array size (assuming
    lmax==mmax).

    Parameters
    ----------
    alm : complex, array
      A complex array of alm. Size must be of the form mmax*(lmax-mmax+1)/2+lmax
    nside : int, scalar
      The nside of the output map.
    lmax : None or int, scalar, optional
      Explicitly define lmax (needed if mmax!=lmax)
    mmax : None or int, scalar, optional
      Explicitly define mmax (needed if mmax!=lmax)
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
      [in radians]

    Returns
    -------
    map : array or list of arrays
      An Healpix map in RING scheme at nside or a list of T,Q,U maps (if
      polarized input)
    """
    smoothalm(alm, fwhm = fwhm, sigma = sigma)
    if pixwin:
        pw=globals()['pixwin'](nside,True)
        if type(alm[0]) is np.ndarray:
            if len(alm) != 3:
                raise TypeError("alm must be a sequence of 3 ndarray "
                                "or a 1D ndarray")
            alm[0] = almxfl(alm[0],pw[0],inplace=True)
            alm[1] = almxfl(alm[1],pw[1],inplace=True)
            alm[2] = almxfl(alm[2],pw[1],inplace=True)
        else:
            alm = almxfl(alm,pw[0],inplace=True)
    if lmax is None:
        lmax = -1
    if mmax is None:
        mmax = -1
    return sphtlib._alm2map(alm, nside, lmax = lmax, mmax = mmax)

def synalm(cls, lmax = None, mmax = None, new = False):
    """Generate a set of alm given cl.
    The cl are given as a float array. Corresponding alm are generated.
    If lmax is None, it is assumed lmax=cl.size-1
    If mmax is None, it is assumed mmax=lmax.

    Parameters
    ----------
    cls : float, array or tuple of arrays
      Either one cl (1D array) or a tuple of either 4 cl (TT,TE,EE,BB)
      or of n*(n+1)/2 cl. Some of the cl may be None, implying no
      cross-correlation. For example, (TT,TE,TB,EE,EB,BB).
      NOTE: this order differs from the alm2cl function !
    lmax : int, scalar, optional
      The lmax (if None or <0, the largest size-1 of cls)
    mmax : int, scalar, optional
      The mmax (if None or <0, =lmax)

    Returns
    -------
    alms : array or list of arrays
      the generated alm if one spectrum is given, or a list of n alms 
      (with n(n+1)/2 the number of input cl, or n=3 if there are 4 input cl).
    """
    if not hasattr(cls, '__len__'):
        raise TypeError('cls must be an array or a sequence of arrays')

    if not hasattr(cls[0], '__len__'):
        if lmax is None or lmax < 0:
            lmax = cls.size-1
        if mmax is None or mmax < 0:
            mmax = lmax
        cls_list = [cls]
        szalm = Alm.getsize(lmax,mmax)
        alm = np.zeros(szalm,'D')
        alm.real = np.random.standard_normal(szalm)
        alm.imag = np.random.standard_normal(szalm)
        alms_list=[alm]
        sphtlib._synalm(cls_list,alms_list,lmax,mmax)
        return alm

    cls_list = list(cls)
    maxsize = max([len(c) for c in cls])

    if lmax is None or lmax < 0:
        lmax = maxsize-1
    if mmax is None or mmax < 0:
        mmax = lmax

    Nspec = sphtlib._getn(len(cls_list))

    if Nspec <= 0:
        if len(cls_list) == 4:
            if new: ## new input order: TT EE BB TE -> TT EE BB TE 0 0
                cls_list = [cls[0], cls[1], cls[2], cls[3], None, None]
            else: ## old input order: TT TE EE BB -> TT TE 0 EE 0 BB
                cls_list = [cls[0], cls[1], None, cls[2], None, cls[3]]
            Nspec = 3
        else:
            raise TypeError("The sequence of arrays must have either 4 elements "
                            "or n(n+1)/2 elements (some may be None)")
    
    print 'Nspec=', Nspec
    
    szalm = Alm.getsize(lmax,mmax)
    alms_list = []
    for i in xrange(Nspec):
        alm = np.zeros(szalm,'D')
        alm.real = np.random.standard_normal(szalm)
        alm.imag = np.random.standard_normal(szalm)
        alms_list.append(alm)
    if new: # new input order: input given by diagonal, should be given by row
        cls_list = new_to_old_spectra_order(cls_list)
    sphtlib._synalm(cls_list, alms_list, lmax, mmax)
    return alms_list

def synfast(cls, nside, lmax = None, mmax = None, alm = False,
            pixwin = False,fwhm = 0.0,sigma = None, new = False):
    """Create a map(s) from cl(s).

    Parameters
    ----------
    cls : array or tuple of array
      A cl or a list of cl (either 4 or 6, see :func:`synalm`)
    nside : int, scalar
      The nside of the output map(s)
    lmax : int, scalar, optional
      Maximum l for alm. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m for alm. Default: 3*nside-1
    alm : bool, scalar, optional
      If True, return also alm(s). Default: False.
    pixwin : bool, scalar, optional
      If True, convolve the alm by the pixel window function. Default: False.
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
      [in radians]

    Returns
    -------
    map : array or tuple of arrays
      The output map (possibly list of maps if polarized input).
      or, if alm is True, a tuple of (map,alm)
      (alm possibly a list of alm if polarized input)
    """
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    if lmax is None or lmax < 0:
        lmax = 3*nside-1
    alms = synalm(cls, lmax = lmax, mmax = mmax, new = new)
    maps = alm2map(alms, nside, lmax, mmax, pixwin = pixwin,
                   fwhm = fwhm, sigma = sigma, new = new)
    if alm:
        return maps,alms
    else:
        return maps
    
class Alm(object):
    """This class provides some static methods for alm index computation.

    Methods
    -------
    getlm
    getidx
    getsize
    getlmax
    """
    def __init__(self):
        pass

    @staticmethod
    def getlm(lmax,i=None):
        """Get the l and m from index and lmax.
        
        Parameters
        ----------
        lmax : int
          The maximum l defining the alm layout
        i : int or None
          The index for which to compute the l and m.
          If None, the function return l and m for i=0..Alm.getsize(lmax)
        """
        if i is None:
            i=np.arange(Alm.getsize(lmax))
        m=(np.ceil(((2*lmax+1)-np.sqrt((2*lmax+1)**2-8*(i-lmax)))/2)).astype(int)
        l = i-m*(2*lmax+1-m)/2
        return (l,m)

    @staticmethod
    def getidx(lmax,l,m):
        """Returns index corresponding to (l,m) in an array describing alm up to lmax.
        
        Parameters
        ----------
        lmax : int
          The maximum l, defines the alm layout
        l : int
          The l for which to get the index
        m : int
          The m for which to get the index

        Returns
        -------
        idx : int
          The index corresponding to (l,m)
        """
        return m*(2*lmax+1-m)/2+l

    @staticmethod
    def getsize(lmax,mmax = None):
        """Returns the size of the array needed to store alm up to *lmax* and *mmax*

        Parameters
        ----------
        lmax : int
          The maximum l, defines the alm layout
        mmax : int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        size : int
          The size of the array needed to store alm up to lmax, mmax.
        """
        if mmax is None or mmax < 0 or mmax > lmax:
            mmax = lmax
        return mmax * (2 * lmax + 1 - mmax) / 2 + lmax + 1

    @staticmethod
    def getlmax(s, mmax = None):
        """Returns the lmax corresponding to a given array size.
        
        Parameters
        ----------
        s : int
          Size of the array
        mmax : None or int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        lmax : int
          The maximum l of the array, or -1 if it is not a valid size.
        """
        if mmax is not None and mmax >= 0:
            x = (2 * s + mmax ** 2 - mmax - 2) / (2 * mmax + 2)
        else:
            x = (-3 + np.sqrt(1 + 8 * s)) / 2
        if x != np.floor(x):
            return -1
        else:
            return int(x)


def alm2cl(alms1, alms2 = None, lmax = None, mmax = None,
           lmax_out = None, nspec = None):
    """Computes (cross-)spectra from alm(s). If alm2 is given, cross-spectra between
    alm and alm2 are computed. If alm (and alm2 if provided) contains n alm,
    then n(n+1)/2 auto and cross-spectra are returned.

    Parameters
    ----------
    alm : complex, array or sequence of arrays
      The alm from which to compute the power spectrum. If n>=2 arrays are given,
      computes both auto- and cross-spectra.
    alm2 : complex, array or sequence of 3 arrays, optional
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
    nspec : None or int, optional
      The number of spectra to return. None means all, otherwise returns cl[:nspec]

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
    cls = _ianafast.alm2cl(alms1, alms2 = alms2, lmax = lmax,
                           mmax = mmax, lmax_out = lmax_out)
    if nspec is None:
        return cls
    else:
        return cls[:nspec]
        
    
def almxfl(alm, fl, mmax = None, inplace = False):
    """Multiply alm by a function of l. The function is assumed
    to be zero where not defined.

    Parameters
    ----------
    alm : array
      The alm to multiply
    fl : array
      The function (at l=0..fl.size-1) by which alm must be multiplied.
    mmax : None or int, optional
      The maximum m defining the alm layout. Default: lmax.
    inplace : bool, optional
      If True, modify the given alm, otherwise make a copy before multiplying.

    Returns
    -------
    alm : array
      The modified alm, either a new array or a reference to input alm, if inplace is True.
    """
    # this is the expected lmax, given mmax
    lmax = Alm.getlmax(alm.size, mmax)
    if lmax < 0:
        raise TypeError('Wrong alm size for the given mmax.')
    if mmax is None or mmax<0:
        mmax=lmax
    fl = np.array(fl)
    if inplace:
        almout = alm
    else:
        almout = alm.copy()
    for l in xrange(lmax+1):
        if l < fl.size:
            a=fl[l]
        else:
            a=0
        i=Alm.getidx(lmax, l, np.arange(min(mmax,l) + 1))
        almout[i] *= a
    return almout

def smoothalm(alm, fwhm = 0.0, sigma = None, mmax = -1,
              verbose = False):
    """Smooth alm with a Gaussian symmetric beam function in place.

    Parameters
    ----------
    alm : array or sequence of 3 arrays
      Either an array representing one alm, or a sequence of
      3 arrays representing 3 alm
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
      [in radians]
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.
      [in radians]
    mmax : None or int, optional
      The maximum m for alm. Default: mmax=lmax
    verbose : bool, optional
      If True prints diagnostic information. Default: False

    Returns
    -------
    None
    """
    if sigma is None:
        sigma = fwhm / (2.*np.sqrt(2.*np.log(2.)))
    if verbose:
        print "Sigma is %f arcmin (%f rad) " %  (sigma*60*180/pi,sigma)
        print "-> fwhm is %f arcmin" % (sigma*60*180/pi*(2.*np.sqrt(2.*np.log(2.))))
    if type(alm[0]) == np.ndarray:
        if len(alm) != 3:
            raise ValueError("alm must be en array or a sequence of 3 arrays")
        retval = []
        for a in alm:
            lmax = Alm.getlmax(a.size, mmax)
            if lmax < 0:
                raise TypeError('Wrong alm size for the given '
                                'mmax (alms[%d]).'%(a.size))
            if mmax is None or mmax < 0:
                mmax=lmax
            ell = np.arange(lmax+1)
            fact = np.exp(-0.5*ell*(ell+1)*sigma**2)
            almxfl(a,fact,mmax,inplace=True)
        return None
    else:
        lmax = Alm.getlmax(alm.size,mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given '
                            'mmax (alms[%d]).'%(a.size))
        if mmax is None or mmax<0:
            mmax=lmax
        ell = np.arange(lmax+1)
        fact = np.exp(-0.5*ell*(ell+1)*sigma**2)
        almxfl(alm,fact,mmax,inplace=True)
        return None

def smoothing(m, fwhm = 0.0, sigma = None):
    """Smooth a map with a Gaussian symmetric beam.

    Parameters
    ----------
    map : array or sequence of 3 arrays
      Either an array representing one map, or a sequence of
      3 arrays representing 3 maps
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.

    Returns
    -------
    map_smo : array or tuple of 3 arrays
      The smoothed map(s)
    """
    if type(m[0]) is np.ndarray:
        if len(m) != 3:
            raise TypeError("map must be en array or a list of 3 arrays")
        nside = pixelfunc.npix2nside(m[0].size)
        if (pixelfunc.npix2nside(m[1].size) != nside
            or pixelfunc.npix2nside(m[2].size) != nside):
            raise TypeError("all maps in the array must have identical nside")
    elif type(m) == np.ndarray:
        nside=pixelfunc.npix2nside(m.size)
    else:
        raise TypeError("map must be en array or a list of 3 arrays")
    # Replace UNSEEN pixels with zeros
    mask = mask_bad(m)
    m[mask] = 0
    alm = map2alm(m)
    return alm2map(alm,nside,fwhm=fwhm,sigma=sigma)


def pixwin(nside, pol = False):
    """Return the pixel window function for the given nside.

    Parameters
    ----------
    nside : int
      The nside for which to return the pixel window function
    pol : bool, optional
      If True, return also the polar pixel window. Default: False

    Returns
    -------
    pw or pwT,pwP : array or tuple of 2 arrays
      The temperature pixel window function, or a tuple with both
      temperature and polarisation pixel window functions.
    """
    datapath = DATAPATH
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    fname = os.path.join(datapath, 'pixel_window_n%04d.fits'%nside)
    if not os.path.isfile(fname):
        raise ValueError("No pixel window for this nside "
                         "or data files missing")
    # return hfitslib._pixwin(nside,datapath,pol)  ## BROKEN -> seg fault...
    try:
        import pyfits
    except ImportError:
        print "*********************************************************"
        print "**   You need to install pyfits to use this function   **"
        print "*********************************************************"
        raise
    pw = pyfits.getdata(fname)
    pw_temp, pw_pol = pw.field(0), pw.field(1)
    if pol:
        return pw_temp, pw_pol
    else:
        return pw_temp

def alm2map_der1(alm, nside, lmax = None, mmax = None):
   """Computes an Healpix map and its first derivatives given the alm.

   The alm are given as a complex array. You can specify lmax
   and mmax, or they will be computed from array size (assuming
   lmax==mmax).

   Parameters
   ----------
   alm : array, complex
     A complex array of alm. Size must be of the form mmax(lmax-mmax+1)/2+lmax
   nside : int
     The nside of the output map.
   lmax : None or int, optional
     Explicitly define lmax (needed if mmax!=lmax)
   mmax : None or int, optional
     Explicitly define mmax (needed if mmax!=lmax)

   Returns
   -------
   m, d_theta, d_phi : tuple of arrays
     The maps correponding to alm, and its derivatives with respect to
     theta and phi.
   """
   return sphtlib._alm2map_der1(alm,nside,lmax=lmax,mmax=mmax)

# Helper function : get nside from m, an array or a sequence
# of arrays
def _get_nside(m):
    if hasattr(m,'__len__'):
        if len(m) == 0:
            raise TypeError('Empty sequence !')
        if hasattr(m[0],'__len__'):
            nside=pixelfunc.npix2nside(len(m[0]))
        else:
            nside=pixelfunc.npix2nside(len(m))
    else:
        raise TypeError('You must give an array or a tuple of arrays '
                        'as input')
    return nside

def new_to_old_spectra_order(cls_new_order):
    """Reorder the cls from new order (by diagonal) to old order (by row).
    For example : TT, EE, BB, TE, EB, BB => TT, TE, TB, EE, EB, BB
    """
    Nspec = sphtlib._getn(len(cls_new_order))
    cls_old_order = []
    for i in xrange(Nspec):
        for j in xrange(i, Nspec):
            p = j - i
            q = i
            idx_new = p * (2 * Nspec + 1 - p) / 2 + q
            cls_old_order.append(cls_new_order[idx_new])
    return cls_old_order

def load_sample_spectra():
    """Read a sample power spectra for testing and demo purpose.
    Based on LambdaCDM. Gives TT, EE, BB, TE.

    Returns
    -------
    ell, f, cls : arrays
      ell is the array of ell values (from 0 to lmax)
      f is the factor ell*(ell+1)/2pi (in general, plots show f * cl)
      cls is a sequence of the power spectra TT, EE, BB and TE
    """
    cls = np.loadtxt(os.path.join(DATAPATH, 'totcls.dat'), unpack = True)
    ell = cls[0]
    f = ell * (ell + 1) / 2 / np.pi
    cls[1:, 1:] /= f[1:]
    return ell, f, cls[1:]


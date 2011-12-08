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
import numpy as npy
pi = npy.pi

import _healpy_sph_transform_lib as sphtlib
import _healpy_fitsio_lib as hfitslib

import os.path
import pixelfunc

from pixelfunc import mask_bad, maptype, UNSEEN

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Spherical harmonics transformation
def anafast(m,lmax=None,mmax=None,iter=1,alm=False, use_weights=False, regression=True):
    """Computes the power spectrum of an Healpix map.

    Parameters
    ----------
    m : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps
    lmax : int, scalar, optional
      Maximum l of the power spectrum (default: 3*nside-1)
    mmax : int, scalar, optional
      Maximum m of the alm (default: lmax)
    iter : int, scalar, optional
      Number of iteration (default: 1)
    alm : bool, scalar, optional
      If True, returns both cl and alm, otherwise only cl is returned
    regression : bool, scalar, optional
      If True, map average is removed before computing alm. Default: True.
    
    Returns
    -------
    res : array or sequence of arrays
      If *alm* is False, returns cl or a list of cl (TT, EE, BB, TE for
      polarized input map)
      Otherwise, returns a tuple (cl, alm), where cl is as above and
      alm is the spherical harmonic transform or a list of almT, almE, almB
      for polarized input
    """
    datapath = DATAPATH #os.path.dirname(__file__)+'/data'
    nside = _get_nside(m)
    if lmax is None:
        lmax = 3*nside-1
    if mmax is None or mmax < 0 or mmax > lmax:
        mmax = lmax
    # Check the presence of weights file
    if use_weights:
        weightfile = 'weight_ring_n%05d.fits' % (nside)
        if not os.path.isfile(datapath+'/'+weightfile):
            raise IOError('File not found : '+datapath+'/'+weightfile)
    # Replace UNSEEN pixels with zeros
    info = maptype(m)
    if info == 0:
        m[m == UNSEEN] = 0
    elif info == 3:
        mi, mq, mu = m
        mask = mask_bad(mi)
        mask |= mask_bad(mq)
        mask |= mask_bad(mu)
        mi[mask] = 0
        mq[mask] = 0
        mu[mask] = 0
    else:
        raise TypeError("Input map must be an array or a sequence of 3 arrays (for polarization)")
    #m[m == UNSEEN] = 0
    clout,almout = sphtlib._map2alm(m,lmax=lmax,mmax=mmax,iter=iter,cl=True,
                                    use_weights=use_weights,data_path=datapath,
                                    regression=regression)
    if alm:
        return (clout,almout)
    else:
        return clout

def map2alm(m,lmax=None,mmax=None,iter=1,use_weights=False,regression=True):
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
    alm : array or tuple of array
      alm or a tuple of 3 alm (almT, almE, almB) if polarized input.
    """
    datapath = DATAPATH #os.path.dirname(__file__)+'/data'
    nside = _get_nside(m)
    if lmax is None:
        lmax = 3*nside-1
    if mmax is None or mmax < 0 or mmax > lmax:
        mmax = lmax
    # Replace UNSEEN pixels with zeros
    m[m == UNSEEN] = 0
    # Check the presence of weights file
    if use_weights:
        weightfile = 'weight_ring_n%05d.fits' % (nside)
        if not os.path.isfile(datapath+'/'+weightfile):
            raise IOError('File not found : '+datapath+'/'+weightfile)
    alm = sphtlib._map2alm(m,lmax=lmax,mmax=mmax,cl=False,
                           iter=iter,
                           use_weights=use_weights,data_path=datapath,
                           regression=regression)
    return alm


def alm2map(alm, nside, lmax=-1, mmax=-1,pixwin=False,
            fwhm=0.0,sigma=None,degree=False,arcmin=False):
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
    lmax : int, scalar, optional
      Explicitly define lmax (needed if mmax!=lmax)
    mmax : int, scalar, optional
      Explicitly define mmax (needed if mmax!=lmax)
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
    degree : bool, scalar, optional
      If True, unit of sigma or fwhm is degree, otherwise it is radian
    arcmin : bool, scalar, optional
      If True, unit of sigma or fwhm is arcmin, otherwise it is radian

    Returns
    -------
    map : array or list of arrays
      An Healpix map in RING scheme at nside or a list of T,Q,U maps (if
      polarized input)
    """
    smoothalm(alm,fwhm=fwhm,sigma=sigma,degree=degree,arcmin=arcmin)
    if pixwin:
        pw=globals()['pixwin'](nside,True)
        if type(alm[0]) is npy.ndarray:
            if len(alm) != 3:
                raise TypeError("alm must be a sequence of 3 ndarray "
                                "or a 1D ndarray")
            alm[0] = almxfl(alm[0],pw[0],inplace=True)
            alm[1] = almxfl(alm[1],pw[1],inplace=True)
            alm[2] = almxfl(alm[2],pw[1],inplace=True)
        else:
            alm = almxfl(alm,pw[0],inplace=True)
    return sphtlib._alm2map(alm,nside,lmax=lmax,mmax=mmax)

def synalm(cls, lmax=-1, mmax=-1):
    """Generate a set of alm given cl.
    The cl are given as a float array. Corresponding alm are generated.
    If lmax is not given or negative, it is assumed lmax=cl.size-1
    If mmax is not given or negative, it is assumed mmax=lmax.

    Parameters
    ----------
    cls : float, array or tuple of arrays
      Either one cl (1D array) or a tuple of either 4 cl (TT,TE,EE,BB)
      or of n*(n+1)/2 cl. Some of the cl may be None, implying no
      cross-correlation. For example, (TT,TE,TB,EE,EB,BB).
      NOTE: this order differs from the alm2cl function !
    lmax : int, scalar, optional
      The lmax (if <0, the largest size-1 of cls)
    mmax : int, scalar, optional
      The mmax (if <0, =lmax)

    Returns
    -------
    : a list of n alms (with n(n+1)/2 the number of input cl,
            or n=3 if there are 4 input cl).
    """
    if not isinstance(cls[0], npy.ndarray):
        if lmax < 0: lmax = cls.size-1
        if mmax < 0: mmax = lmax
        cls_list = [cls]
        szalm = Alm.getsize(lmax,mmax)
        alm = npy.zeros(szalm,'D')
        alm.real = npy.random.standard_normal(szalm)
        alm.imag = npy.random.standard_normal(szalm)
        alms_list=[alm]
        sphtlib._synalm(cls_list,alms_list,lmax,mmax)
        return alm
    # otherwise, cls must be a sequence of arrays
    try:
        cls_list = list(cls)
        maxsize = 0
        for c in cls_list:
            if c is not None:
                if c.size > maxsize: maxsize=c.size
    except:
        raise TypeError("First argument must be an array or a sequence of arrays.")
    if lmax < 0: lmax = maxsize-1
    if mmax < 0: mmax = lmax
    if sphtlib._getn(len(cls_list)) <= 0:
        if len(cls_list) == 4:  # if 4 cls are given, assume TT, TE, EE, BB
            cls_list = [cls[0], cls[1], None, cls[2], None, cls[3]]
        else:
            raise TypeError("The sequence of arrays must have either 4 elements "
                            "(TT,TE,EE,BB)\n"
                            "or n(n+1)/2 elements (some may be None)")
    szalm = Alm.getsize(lmax,mmax)
    alms_list = []
    for i in xrange(sphtlib._getn(len(cls_list))):
        alm = npy.zeros(szalm,'D')
        alm.real = npy.random.standard_normal(szalm)
        alm.imag = npy.random.standard_normal(szalm)
        alms_list.append(alm)
    sphtlib._synalm(cls_list, alms_list, lmax, mmax)
    if len(alms_list) > 1:
        return alms_list
    else:
        return alms_list[0]

def synfast(cls,nside,lmax=-1,mmax=-1,alm=False,
            pixwin=False,fwhm=0.0,sigma=None,degree=False,
            arcmin=False):
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
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
    degree : bool, scalar, optional
      If True, unit of sigma or fwhm is degree, otherwise it is radian
    arcmin : bool, scalar, optional
      If True, unit of sigma or fwhm is arcmin, otherwise it is radian

    Returns
    -------
    map : array or tuple of arrays
      The output map (possibly list of maps if polarized input).
      or, if alm is True, a tuple of (map,alm)
      (alm possibly a list of alm if polarized input)
    """
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    if lmax < 0:
        lmax = 3*nside-1
    alms = synalm(cls,lmax,mmax)
    maps = alm2map(alms,nside,lmax,mmax,pixwin=pixwin,
                   fwhm=fwhm,sigma=sigma,degree=degree,
                   arcmin=arcmin)
    if alm: return maps,alms
    else: return maps
    
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
            i=npy.arange(Alm.getsize(lmax))
        m=(npy.ceil(((2*lmax+1)-npy.sqrt((2*lmax+1)**2-8*(i-lmax)))/2)).astype(int)
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
    def getsize(lmax,mmax=-1):
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
        if mmax<0 or mmax > lmax:
            mmax=lmax
        return mmax*(2*lmax+1-mmax)/2+lmax+1

    @staticmethod
    def getlmax(s,mmax=-1):
        """Returns the lmax corresponding to a given array size.

        Parameters
        ----------
        s : int
          Size of the array
        mmax : int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        lmax : int
          The maximum l of the array, or -1 if it is not a valid size.
        """
        if mmax >= 0:
            x=(2*s+mmax**2-mmax-2)/(2*mmax+2)
        else:
            x=(-3+npy.sqrt(1+8*s))/2
        if x != npy.floor(x):
            return -1
        else:
            return int(x)



def alm2cl(alm,mmax=-1,nspec=4):
    """Compute the auto- and cross- spectra of the given alm's
    (either 1 alm or 3 alm).

    Parameters
    ----------
    alm : array or sequence of arrays
      one array or a sequence of 3 arrays of identical size
    mmax : int, optional
      The maximum m for alm(s)
    nspec : int, optional
      The number of spectra to return if 3 alms were given:
      nspec==[0-6], in order TT,EE,BB,TE,TB,EB. 

    Returns
    -------
    cl : array or sequence of arrays
      the power spectrum estimated from alm.
    """
    # this is the expected lmax, given mmax
    if isinstance(alm,npy.ndarray) and alm.ndim == 1:
        lmax = Alm.getlmax(alm.size,mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given mmax.')
        if mmax<0:
            mmax=lmax
        cl_est = npy.zeros(lmax+1)
        for l in range(lmax+1):
            i=Alm.getidx(lmax,l,npy.arange(min(mmax,l)+1))
            cl_est[l] = (npy.abs(alm[i[0]])**2
                         +2.*npy.sum(npy.abs(alm[i[1:]])**2))/(2*l+1)
        return cl_est
    else:
        almT,almE,almB=tuple(alm)
        lmax = Alm.getlmax(almT.size,mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given mmax.')
        if mmax<0:
            mmax=lmax
        ctt_est = npy.zeros(lmax+1)
        cee_est = npy.zeros(lmax+1)
        cbb_est = npy.zeros(lmax+1)
        cte_est = npy.zeros(lmax+1)
        ctb_est = npy.zeros(lmax+1)
        ceb_est = npy.zeros(lmax+1)
        for l in range(lmax+1):
            i=Alm.getidx(lmax,l,npy.arange(min(mmax,l)+1))
            ctt_est[l] = (npy.abs(almT[i[0]])**2
                          +2.*npy.sum(npy.abs(almT[i[1:]])**2))/(2*l+1)
            cee_est[l] = (npy.abs(almE[i[0]])**2
                          +2.*npy.sum(npy.abs(almE[i[1:]])**2))/(2*l+1)
            cbb_est[l] = (npy.abs(almB[i[0]])**2
                          +2.*npy.sum(npy.abs(almB[i[1:]])**2))/(2*l+1)
            cte_est[l] = (almT[i[0]]*almE[i[0]].conj()
                          +2.*npy.sum(almT[i[1:]]*almE[i[1:]].conj()))/(2*l+1)
            ctb_est[l] = (almT[i[0]]*almB[i[0]].conj()
                          +2.*npy.sum(almT[i[1:]]*almB[i[1:]].conj()))/(2*l+1)
            ceb_est[l] = (almE[i[0]]*almB[i[0]].conj()
                          +2.*npy.sum(almE[i[1:]]*almB[i[1:]].conj()))/(2*l+1)
        return (ctt_est,cee_est,cbb_est,cte_est,ctb_est,ceb_est)[:nspec]
        
    
def almxfl(alm,fl,mmax=-1,inplace=False):
    """Multiply alm by a function of l. The function is assumed
    to be zero where not defined.

    Parameters
    ----------
    alm : array
      The alm to multiply
    fl : array
      The function (at l=0..fl.size-1) by which alm must be multiplied.
    mmax : int, optional
      The maximum m defining the alm layout. Default: lmax.
    inplace : bool, optional
      If True, modify the given alm, otherwise make a copy before multiplying.

    Returns
    -------
    alm : array
      The modified alm, either a new array or a reference to input alm, if inplace is True.
    """
    # this is the expected lmax, given mmax
    lmax = Alm.getlmax(alm.size,mmax)
    if lmax < 0:
        raise TypeError('Wrong alm size for the given mmax.')
    if mmax<0:
        mmax=lmax
    fl = npy.array(fl)
    if inplace:
        almout = alm
    else:
        almout = alm.copy()
    for l in xrange(lmax+1):
        if l < fl.size:
            a=fl[l]
        else:
            a=0
        i=Alm.getidx(lmax,l,npy.arange(min(mmax,l)+1))
        almout[i] *= a
    return almout

def smoothalm(alm,fwhm=0.0,sigma=None,degree=False,
              arcmin=False,mmax=-1,verbose=False):
    """Smooth alm with a Gaussian symmetric beam function in place.

    Parameters
    ----------
    alm : array or sequence of 3 arrays
      Either an array representing one alm, or a sequence of
      3 arrays representing 3 alm
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.
    degree : bool, optional
      If True, parameter given in degree. Override arcmin. Default: False
    arcmin : bool, optional
      If True, parameter given in arcmin. Default: False
    mmax : int, optional
      The maximum m for alm. Default: mmax=lmax
    verbose : bool, optional
      If True prints diagnostic information. Default: False

    Returns
    -------
    None
    """
    if sigma is None:
        sigma = fwhm / (2.*npy.sqrt(2.*npy.log(2.)))
    if degree:
        sigma *= (pi/180.)
    elif arcmin:
        sigma *= (pi/180./60.)
    if verbose:
        print "Sigma is %f arcmin (%f rad) " %  (sigma*60*180/pi,sigma)
        print "-> fwhm is %f arcmin" % (sigma*60*180/pi*(2.*npy.sqrt(2.*npy.log(2.))))
    if type(alm[0]) == npy.ndarray:
        if len(alm) != 3:
            raise ValueError("alm must be en array or a sequence of 3 arrays")
        retval = []
        for a in alm:
            lmax = Alm.getlmax(a.size,mmax)
            if lmax < 0:
                raise TypeError('Wrong alm size for the given '
                                'mmax (alms[%d]).'%(a.size))
            if mmax<0:
                mmax=lmax
            ell = npy.arange(lmax+1)
            fact = npy.exp(-0.5*ell*(ell+1)*sigma**2)
            almxfl(a,fact,mmax,inplace=True)
        return None
    else:
        lmax = Alm.getlmax(alm.size,mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given '
                            'mmax (alms[%d]).'%(a.size))
        if mmax<0:
            mmax=lmax
        ell = npy.arange(lmax+1)
        fact = npy.exp(-0.5*ell*(ell+1)*sigma**2)
        almxfl(alm,fact,mmax,inplace=True)
        return None

def smoothing(m,fwhm=0.0,sigma=None,degree=False,
              arcmin=False):
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
    degree : bool, optional
      If True, parameter given in degree. Override arcmin. Default: False
    arcmin : bool, optional
      If True, parameter given in arcmin. Default: False
    mmax : int, optional
      The maximum m for alm. Default: mmax=lmax

    Returns
    -------
    map_smo : array or tuple of 3 arrays
      The smoothed map(s)
    """
    if type(m[0]) is npy.ndarray:
        if len(m) != 3:
            raise TypeError("map must be en array or a list of 3 arrays")
        nside = pixelfunc.npix2nside(m[0].size)
        if (pixelfunc.npix2nside(m[1].size) != nside
            or pixelfunc.npix2nside(m[2].size) != nside):
            raise TypeError("all maps in the array must have identical nside")
    elif type(m) == npy.ndarray:
        nside=pixelfunc.npix2nside(m.size)
    else:
        raise TypeError("map must be en array or a list of 3 arrays")
    # Replace UNSEEN pixels with zeros
    m[m == UNSEEN] = 0
    alm = map2alm(m)
    return alm2map(alm,nside,fwhm=fwhm,sigma=sigma,
                   degree=degree,arcmin=arcmin)

def pixwin(nside,pol=False):
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

def alm2map_der1(alm, nside, lmax=-1, mmax=-1):
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
   lmax : int, optional
     Explicitly define lmax (needed if mmax!=lmax)
   mmax : int, optional
     Explicitly define mmax (needed if mmax!=lmax)

   Returns
   -------
   m, d_theta, d_phi : tuple of arrays
     The maps correponding to alm, and its derivatives with respect to theta and phi.
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

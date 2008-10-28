import numpy as npy
import _healpy_sph_transform_lib as sphtlib
import _healpy_fitsio_lib as hfitslib
import os.path
import pixelfunc

pi = npy.pi

# Spherical harmonics transformation
def anafast(m,lmax=None,mmax=None,iter=3,alm=False, use_weights=False):
    """Computes the power spectrum of an Healpix map.

    Input:
      - m : either an array representing a map, or a list of 3 arrays
            representing I, Q, U maps
    Parameters:
      - lmax : maximum l of the power spectrum (default: 3*nside-1)
      - mmax : maximum m of the alm (default: lmax)
      - iter : number of iteration (default: 3)
      - alm : (boolean) whether to return alm or not (if True, both are
               returned in a tuple)
    Return:
      - if alm==False: return cl or a list of cl's (TT,TE,EE,BB)
      - if alm==True: return a tuple with cl or a list of cl's and alm
                      or a list of alm's
    """
    datapath=os.path.dirname(__file__)+'/data'
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
    clout,almout = sphtlib._map2alm(m,lmax=lmax,mmax=mmax,iter=iter,cl=True,
                                    use_weights=use_weights,data_path=datapath)
    if alm:
        return (clout,almout)
    else:
        return clout

def map2alm(m,lmax=None,mmax=None,iter=3,use_weights=False):
    """Computes the alm of an Healpix map.

    Input:
      - m: a ndarray (not polarised) or a list of 3 ndarray (polarised)
    Parameters:
      - lmax : maximum l of the power spectrum. Default: 3*nside-1
      - mmax : maximum m of the alm. Default: lmax
      - iter : number of iteration (default: 3)
      - use_weights: whether to use ring weights or not. Default: False.
    Return:
      - alm as one ndarray or a tuple of 3 ndarrays
    """
    datapath=os.path.dirname(__file__)+'/data'
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
    alm = sphtlib._map2alm(m,lmax=lmax,mmax=mmax,cl=False,
                           iter=iter,
                           use_weights=use_weights,data_path=datapath)
    return alm


def alm2map(alm, nside, lmax=-1, mmax=-1,pixwin=False,
            fwhm=0.0,sigma=None,degree=False,arcmin=False):
    """Computes an Healpix map given the alm.

    The alm are given as a complex array. You can specitify lmax
    and mmax, or they will be computed from array size (assuming
    lmax==mmax).

    Parameters:
    - alm: a complex array of alm. Size must be of the form
           size=mmax(lmax-mmax+1)/2+lmax
    - nside: the nside of the output map.
    - lmax: explicitly define lmax (needed if mmax!=lmax)
    - mmax: explicitly define mmax (needed if mmax!=lmax)
    - fwhm, sigma, degree and arcmin (as in smoothalm): smooth by a gaussian
      symmetric beam

    Return: an Healpix map in RING scheme at nside.
    """
    smoothalm(alm,fwhm=fwhm,sigma=sigma,degree=degree,arcmin=arcmin)
    if pixwin:
        pw=globals()['pixwin'](nside,True)
        if type(alm[0]) is ndarray:
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

    Call: alms = synalm(cls, lmax=-1, mmax=-1)

    Input:
      - cls: either one cl (1D array) or a tuple of either 4 cl (TT,TE,EE,BB)
             or of n(n+1)/2 cl. Some of the cl may be None, implying no
             cross-correlation. For example, (TT,TE,TB,EE,EB,BB).
             NOTE: this order differs from the alm2cl function !
    Parameters:
      - lmax: the lmax (if <0, the largest size-1 of cls)
      - mmax: the mmax (if <0, =lmax)

    Return: a list of n alms (with n(n+1)/2 the number of input cl,
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

    Input:
      - cls: an array of cl or a list of cls (either 4 or 6, see synalm)
      - nside: the nside of the output map(s)
    Parameters:
      - lmax, mmax: maximum l and m for alm. Default: 3*nside-1
      - alm : if True, return also alm(s). Default: False.
      - pixwin: convolve the alm by the pixel window function. Default: False.
      - fwhm,sigma,degree,arcmin: see smoothalm. Convolve the map(s)
        by a symmetric gaussian beam
    Output:
      - if alm==False: return a map or a tuple of maps
      - if alm==True: return a tuple of map(s) and alm(s)
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
    * getlm(lmax,i=None)
    * getidx(lmax,l,m)
    * getsize(lmax,mmax=-1)
    * getlmax(s,mmax=-1)
    """
    def __init__(self,lmax):
        pass

    @staticmethod
    def getlm(lmax,i=None):
        """Get the l and m from index and lmax.
        
        Parameters:
        - lmax
        - i: the index. If not given, the function return l and m
             for i=0..Alm.getsize(lmax)
        """
        if i is None:
            i=npy.arange(Alm.getsize(lmax))
        m=(npy.ceil(((2*lmax+1)-npy.sqrt((2*lmax+1)**2-8*(i-lmax)))/2)).astype(int)
        l = i-m*(2*lmax+1-m)/2
        return (l,m)

    @staticmethod
    def getidx(lmax,l,m):
        """Get index from lmax, l and m.
        
        Parameters:
        - lmax
        - l
        - m
        """
        return m*(2*lmax+1-m)/2+l

    @staticmethod
    def getsize(lmax,mmax=-1):
        if mmax<0 or mmax > lmax:
            mmax=lmax
        return mmax*(2*lmax+1-mmax)/2+lmax+1

    @staticmethod
    def getlmax(s,mmax=-1):
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

    Input:
      - alm: one array or a sequence of 3 arrays of identical size
    Parameters:
      - mmax: maximum m for alm(s)
      - nspec: number of spectra to return if 3 alms were given:
        nspec==[0-6], in order TT,EE,BB,TE,TB,EB. 
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
            i=Alm.getidx(lmax,l,arange(min(mmax,l)+1))
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
    If inplace is True, the operation is done in place. Always return
    the alm array (either a new array or the modified input array).
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
    for l in range(lmax+1):
        if l < fl.size:
            a=fl[l]
        else:
            a=0
        i=Alm.getidx(lmax,l,npy.arange(min(mmax,l)+1))
        almout[i] *= a
    return almout

def smoothalm(alm,fwhm=0.0,sigma=None,degree=False,
              arcmin=False,mmax=-1):
    """Smooth alm with a gaussian symmetric beam function in place.

    Input:
      - alm: either an array representing one alm, or a sequence of
             3 arrays representing 3 alm
    Parameters:
      - fwhm: the full width half max parameter of the gaussian. Default:0.0
      - sigma: the sigma of the Gaussian. Override fwhm.
      - degree: if True, parameter given in degree. Override arcmin.
                Default: False
      - arcmin: if True, parameter given in arcmin. Default: False
      - mmax: the maximum m for alm. Default: mmax=lmax
    Return:
      None
    """
    if sigma is None:
        sigma = fwhm / (2.*npy.sqrt(2.*npy.log(2.)))
    if degree:
        sigma *= (pi/180.)
    elif arcmin:
        sigma *= (pi/180./60.)
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
                                'mmax (alms[%d]).'%i)
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
                            'mmax (alms[%d]).'%i)
        if mmax<0:
            mmax=lmax
        ell = npy.arange(lmax+1)
        fact = npy.exp(-0.5*ell*(ell+1)*sigma**2)
        almxfl(alm,fact,mmax,inplace=True)
        return None

def smoothing(m,fwhm=0.0,sigma=None,degree=False,
              arcmin=False):
    """Smooth a map with a gaussian symmetric beam.

    Input:
      - map: either an array representing one map, or a sequence of
             3 arrays representing 3 maps
    Parameters:
      - fwhm: the full width half max parameter of the gaussian. Default:0.0
      - sigma: the sigma of the Gaussian. Override fwhm.
      - degree: if True, parameter given in degree. Override arcmin.
                Default: False
      - arcmin: if True, parameter given in arcmin. Default: False
      - mmax: the maximum m for alm. Default: mmax=lmax
    Return:
      - the smoothed map(s)
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
    alm = map2alm(m)
    return alm2map(alm,nside,fwhm=fwhm,sigma=sigma,
                   degree=degree,arcmin=arcmin)

def pixwin(nside,pol=False):
    """Return the pixel window function for the given nside.

    Input:
      - nside
    Parameters:
      - pol: if True, return also the polar pixel window. Default: False
    Return:
      - the temperature pixel window function, or a tuple with both
        temperature and polarisation pixel window functions.
    """
    datapath=os.path.dirname(__file__)+'/data'
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    if not os.path.isfile(datapath+'/pixel_window_n%04d.fits'%nside):
        raise ValueError("No pixel window for this nside "
                         "or data files missing")
    return hfitslib._pixwin(nside,datapath,pol)

def alm2signal(alm, theta, phi, lmax=-1, mmax=-1):
    """This function use alm (temperature only) to compute the signal at
    an arbitrary direction on the sky, ie computes sum alm Ylm(theta,phi).

    Input:
     - alm: a vector containing alm
     - theta, phi : a direction on the sky
     - lmax, mmax (optionnal): the lmax, mmax of the alm. If not given, assume
                               mmax=lmax and find lmax from alm size 
    Return:
     - a scalar (double) = sum alm Ylm
    """
    sig = sphtlib._alm2signal(npy.asarray(alm,
                                          dtype=npy.comple128,
                                          order='C').ravel(),
                              theta, phi, lmax, mmax)
    return sig

def getylm(lmax, m, theta):
    """Computes Ylm for given theta and phi=0, for l=0..lmax and given m.
    Values for l<|m| are set to zero.

    Input:
     - lmax: maximum l value for which to compute Ylm
     - m: value m for which to compute Ylm
     - theta: value for which Ylm(theta,0) is computed
    Return:
     - ylm : a vector of size lmax+1. ylm[0..m-1] = 0.0
    """
    ylm = sphtlib._getylm(lmax,m,m,theta)
    return ylm

def alm2map_der1(alm, nside, lmax=-1, mmax=-1):
   """Computes an Healpix map and its first derivatives given the alm.

   The alm are given as a complex array. You can specitify lmax
   and mmax, or they will be computed from array size (assuming
   lmax==mmax).

   Parameters:
   - alm: a complex array of alm. Size must be of the form
          size=mmax(lmax-mmax+1)/2+lmax
   - nside: the nside of the output map.
   - lmax: explicitly define lmax (needed if mmax!=lmax)
   - mmax: explicitly define mmax (needed if mmax!=lmax)

   Return: an Healpix map in RING scheme at nside.
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

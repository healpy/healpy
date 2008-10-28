import numpy as npy
import _healpy_pixel_lib as pixlib

def ang2pix(nside,theta,phi,nest=False):
    """ang2pix : nside,theta,phi,nest=False -> ipix (default:RING)
    """
    if nest:
        return pixlib._ang2pix_nest(nside,theta,phi)
    else:
        return pixlib._ang2pix_ring(nside,theta,phi)

def pix2ang(nside,ipix,nest=False):
    """pix2ang : nside,ipix,nest=False -> theta,phi (default RING)
    """
    if nest:
        return pixlib._pix2ang_nest(nside, ipix)
    else:
        return pixlib._pix2ang_ring(nside,ipix)

def vec2pix(nside,x,y,z,nest=False):
    """vec2pix : nside,x,y,z,nest=False -> ipix (default:RING)
    """
    if nest:
        return pixlib._vec2pix_nest(nside,x,y,z)
    else:
        return pixlib._vec2pix_ring(nside,x,y,z)

def pix2vec(nside,ipix,nest=False):
    """pix2vec : nside,ipix,nest=False -> x,y,z (default RING)
    """
    if nest:
        return pixlib._pix2vec_nest(nside,ipix)
    else:
        return pixlib._pix2vec_ring(nside,ipix)

def ring2nest(nside, ipix):
    """Convert pixel number from ring scheme number to nest scheme number.

    Input:
      - nside: the nside to work with
      - ipix: the pixel number (can be an array) in ring scheme
    Return:
      - a pixel number or an array of pixel numbers in nest scheme
    """
    return pixlib._ring2nest(nside, ipix)

def nest2ring(nside, ipix):
    """Convert pixel number from nest scheme number to ring scheme number.

    Input:
      - nside: the nside to work with
      - ipix: the pixel number (can be an array) in nest scheme
    Return:
      - a pixel number or an array of pixel numbers in ring scheme
    """
    return pixlib._nest2ring(nside, ipix)

def nside2npix(nside):
    """Give the number of pixel for the given nside.

    Input:
      - nside: nside paramater
    Return:
      - npix: corresponding number of pixels
    Raise a ValueError exception if nside is not valid.
    """
    if not isnsideok(nside):
        raise ValueError("Given number is not a valid mside parameter "
                         "(must be a power of 2)")
    return 12*nside**2

def npix2nside(npix):
    """Give the nside parameter for the given number of pixels.

    Input:
      - npix: number of pixel
    Return:
      - nside: the nside parameter.
    Raise a ValueError exception if number of pixel does not correspond to
    the number of pixel of an healpix map.
    """
    nside = npy.sqrt(npix/12.)
    if nside != npy.floor(nside):
        raise ValueError("Wrong pixel number (it is not 12*n**2)")
    nside=int(npy.floor(nside))
    if not isnsideok(nside):
        raise ValueError("Wrong pixel number (it is not 12*n**2)")
    return nside

def isnsideok(nside):
    """Return True if nside is a valid nside parameter, False otherwise.
    """
    if( npy.log2(nside) == npy.floor(npy.log2(nside)) ):
        return True
    else:
        return False

def isnpixok(npix):
    if hasattr(npix,'__len__'):
        nside = npy.sqrt(npy.asarray(npix)/12.)
        return (nside == npy.floor(nside))
    else:
        nside = npy.sqrt(npix/12.)
        return (nside == npy.floor(nside))

def get_interp_val(m,theta,phi,nest=False):
    """Return the bi-linear interpolation value of a map at given direction.

    Input:
      - m: a map (an ndarray)
      - theta, phi : the direction (either scalar or arrays of same size)
    Parameters:
      - nest: if True, the map is in NEST scheme. Default: False (ie RING)
    Return:
      - the interpolated value(s)
    """
    m2=m.ravel()
    nside=npix2nside(m2.size)
    if nest:
        r=pixlib._get_interpol_nest(nside,theta,phi)
    else:
        r=pixlib._get_interpol_ring(nside,theta,phi)
    p=npy.array(r[0:4])
    w=npy.array(r[4:8])
    del r
    return npy.sum(m2[p]*w,0)

def get_neighbours(nside,theta,phi=None,nest=False):
    """Return the 4 nearest pixels and the corresponding weights for
    bi-linear interpolation for the given direction.

    Input:
      - nside: the nside to work with
      - theta, phi: if phi is not given, theta is actually a pixel number
                    if phi is given, theta,phi is a direction
    Parameters:
      - nest: if True, NEST scheme. Default: False (RING)
    Return:
      - tuple of pixels and weights
    """
    if not isnsideok(nside):
        raise ValueError('Wrong nside value. Must be a power of 2.')
    if phi == None:
        theta,phi = pix2ang(nside,theta)    
    if nest:
        r=pixlib._get_interpol_nest(nside,theta,phi)
    else:
        r=pixlib._get_interpol_ring(nside,theta,phi)
    p=npy.array(r[0:4])
    w=npy.array(r[4:8])
    return (p,w)

def reorder(map_in, inp=None, out=None, r2n=None, n2r=None):
    """Reorder an healpix map from RING/NESTED ordering to NESTED/RING

    Input:
      - map_in: the input map (can be a sequence of maps)
    Parameters:
      - inp: the ordering of the input map 'RING' or 'NESTED'
      - out: the ordering of the output map 'RING' or 'NRSTED'
    Output:
      - map_out: the reordered map
    """
    pass
    
def fit_dipole(m,nest=False,bad=pixlib.UNSEEN,gal_cut=0):
    """Fit a dipole and a monopole to the map, excluding unseen pixels.
    Input:
      - m: the map from which a dipole is fitted and subtracted
      - nest: False if m is in RING scheme, True if it is NESTED
      - bad: bad values of pixel, default to UNSEEN.
      - gal_cut: latitude below which pixel are not taken into account
    Return:
      - the monopole value  and the dipole vector
    """
    m=npy.asarray(m)
    npix = m.size
    nside = npix2nside(npix)
    if nside>128:
        bunchsize = npix/24
    else:
        bunchsize = npix
    aa = npy.zeros((4,4),dtype=npy.float64)
    v = npy.zeros(4,dtype=npy.float64)
    for ibunch in range(npix/bunchsize):
        ipix = npy.arange(ibunch*bunchsize,
                          (ibunch+1)*bunchsize)
        ipix = ipix[m.flat[ipix]!=bad]
        x,y,z = pix2vec(nside, ipix, nest)
        if gal_cut>0:
            w = (npy.abs(z)>=npy.sin(gal_cut*npy.pi/180))
            ipix=ipix[w]
            x=x[w]
            y=y[w]
            z=z[w]
            del w
        aa[0,0] += ipix.size
        aa[1,0] += x.sum()
        aa[2,0] += y.sum()
        aa[3,0] += z.sum()
        aa[1,1] += (x**2).sum()
        aa[2,1] += (x*y).sum()
        aa[3,1] += (x*z).sum()
        aa[2,2] += (y**2).sum()
        aa[3,2] += (y*z).sum()
        aa[3,3] += (z**2).sum()
        v[0] += m.flat[ipix].sum()
        v[1] += (m.flat[ipix]*x).sum()
        v[2] += (m.flat[ipix]*y).sum()
        v[3] += (m.flat[ipix]*z).sum()
    aa[0,1] = aa[1,0]
    aa[0,2] = aa[2,0]
    aa[0,3] = aa[3,0]
    aa[1,2] = aa[2,1]
    aa[1,3] = aa[3,1]
    aa[2,3] = aa[3,2]
    res = npy.dot(npy.linalg.inv(aa),v)
    mono = res[0]
    dipole = res[1:4]
    return mono,dipole

def remove_dipole(m,nest=False,bad=pixlib.UNSEEN,gal_cut=0,fitval=False,
                  copy=True,verbose=False):
    """Fit and subtract the dipole and the monopole from the given map m.
    Input:
      - m: the map from which a dipole is fitted and subtracted
      - nest: False if m is in RING scheme, True if it is NESTED
      - bad: bad values of pixel, default to UNSEEN.
      - gal_cut: latitude below which pixel are not taken into account
      - fitval: whether to return or not the fitted values of monopole and dipole
    Return:
      if fitval is False:
      - the map with monopole and dipole subtracted
      if fitval is True:
      - the map, the monopole value  and the dipole vector
    """
    m=npy.array(m,copy=copy)
    npix = m.size
    nside = npix2nside(npix)
    if nside>128:
        bunchsize = npix/24
    else:
        bunchsize = npix
    mono,dipole = fit_dipole(m,nest=nest,bad=bad,gal_cut=gal_cut)
    for ibunch in range(npix/bunchsize):
        ipix = npy.arange(ibunch*bunchsize,
                          (ibunch+1)*bunchsize)
        ipix = ipix[m.flat[ipix]!=bad]
        x,y,z = pix2vec(nside, ipix, nest)
        m.flat[ipix] -= (dipole[0]*x)
        m.flat[ipix] -= dipole[1]*y
        m.flat[ipix] -= dipole[2]*z
        m.flat[ipix] -= mono
    if verbose:
        import rotator as R
        lon,lat = R.vec2dir(dipole,lonlat=True)
        amp = npy.sqrt((dipole**2).sum())
        print 'monopole: %g  dipole: lon: %g, lat: %g, amp: %g'%(mono,
                                                                 lon,
                                                                 lat,
                                                                 amp)
    if fitval:
        return m,mono,dipole
    else:
        return m

def fit_monopole(m,nest=False,bad=pixlib.UNSEEN,gal_cut=0):
    """Fit a monopole to the map, excluding unseen pixels.
    Input:
      - m: the map from which a dipole is fitted and subtracted
      - nest: False if m is in RING scheme, True if it is NESTED
      - bad: bad values of pixel, default to UNSEEN.
      - gal_cut: latitude below which pixel are not taken into account
    Return:
      - the monopole value
    """
    m=npy.asarray(m)
    npix=m.size
    nside = npix2nside(npix)
    if nside>128:
        bunchsize=npix/24
    else:
        bunchsize=npix
    for ibunch in range(npix/bunchsize):
        ipix = npy.arange(ibunch*bunchsize,
                          (ibunch+1)*bunchsize)
        ipix = ipix[m.flat[ipix]!=bad]
        x,y,z = pix2vec(nside, ipix, nest)
        aa = v = 0.0
        if gal_cut>0:
            w = (npy.abs(z)>=npy.sin(gal_cut*npy.pi/180))
            ipix=ipix[w]
            x=x[w]
            y=y[w]
            z=z[w]
            del w
        aa += ipix.size
        v += m.flat[ipix].sum()
    mono = v/aa
    return mono

def remove_monopole(m,nest=False,bad=pixlib.UNSEEN,gal_cut=0,fitval=False,
                    copy=True,verbose=False):
    """Fit and subtract the monopole from the given map m.
    Input:
      - m: the map from which a dipole is fitted and subtracted
      - nest: False if m is in RING scheme, True if it is NESTED
      - bad: bad values of pixel, default to UNSEEN.
      - gal_cut: latitude below which pixel are not taken into account
      - fitval: whether to return or not the fitted values of monopole
    Return:
      if fitval is False:
      - the map with monopole subtracted
      if fitval is True:
      - the map with monopole subtracted and the monopole value
    """
    m=npy.array(m,copy=copy)
    npix = m.size
    nside = npix2nside(npix)
    if nside>128:
        bunchsize = npix/24
    else:
        bunchsize = npix
    mono = fit_monopole(m,nest=nest,bad=bad,gal_cut=gal_cut)
    for ibunch in range(npix/bunchsize):
        ipix = npy.arange(ibunch*bunchsize,
                          (ibunch+1)*bunchsize)
        ipix = ipix[m.flat[ipix]!=bad]
        x,y,z = pix2vec(nside, ipix, nest)
        m.flat[ipix] -= mono
    if verbose:
        print 'monopole: %g'%mono
    if fitval:
        return m,mono
    else:
        return m

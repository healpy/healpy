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


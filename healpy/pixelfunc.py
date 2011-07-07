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
"""
=====================================================
pixelfunc.py : Healpix pixelization related functions
=====================================================

This module provides functions related to Healpix pixelization scheme.

conversion from/to sky coordinates
----------------------------------

- :func:`pix2ang` converts pixel number to angular coordinates
- :func:`pix2vec` converts pixel number to unit 3-vector direction
- :func:`ang2pix` converts angular coordinates to pixel number
- :func:`vec2pix` converts 3-vector to pixel number 
- :func:`vec2ang` converts 3-vector to angular coordinates
- :func:`ang2vec` converts angular coordinates to unit 3-vector
- :func:`get_neighbours` returns the 4 nearest pixels for given
  angular coordinates
- :func:`get_all_neighbours` return the 8 nearest pixels for given
  angular coordinates

conversion between NESTED and RING schemes
------------------------------------------

- :func:`nest2ring` converts NESTED scheme pixel numbers to RING
  scheme pixel number
- :func:`ring2nest` converts RING scheme pixel number to NESTED
  scheme pixel number
- :func:`reorder` reorders a healpix map pixels from one scheme to another

nside/npix/resolution
---------------------

- :func:`nside2npix` converts healpix nside parameter to number of pixel
- :func:`npix2nside` converts number of pixel to healpix nside parameter
- :func:`nside2resol` converts nside to mean angular resolution
- :func:`nside2pixarea` converts nside to pixel area
- :func:`isnsideok` checks the validity of nside
- :func:`isnpixok` checks the validity of npix
- :func:`get_map_size` gives the number of pixel of a map
- :func:`get_min_valid_nside` gives the minimum nside possible for a given
  number of pixel
- :func:`get_nside` returns the nside of a map
- :func:`maptype` checks the type of a map (one map or sequence of maps)
- :func:`ud_grade` upgrades or degrades the resolution (nside) of a map

Masking pixels
--------------

- :const:`UNSEEN` is a constant value interpreted as a masked pixel
- :func:`mask_bad` returns a map with ``True`` where map is :const:`UNSEEN`
- :func:`mask_good` returns a map with ``False`` where map is :const:`UNSEEN`

Map data manipulation
---------------------

- :func:`fit_dipole` fits a monopole+dipole on the map
- :func:`fit_monopole` fits a monopole on the map
- :func:`remove_dipole` fits and removes a monopole+dipole from the map
- :func:`remove_monopole` fits and remove a monopole from the map
- :func:`get_interp_val` computes a bilinear interpolation of the map
  at given angular coordinates, using 4 nearest neighbours
"""

import numpy as npy
import _healpy_pixel_lib as pixlib
from _healpy_pixel_lib import UNSEEN
import exceptions

def mask_bad(m, badval = UNSEEN, rtol = 1.e-5, atol = 1.e-8):
    """Returns a bool array with ``True`` where m is close to badval.

    Parameters
    ----------
    m : a map (may be a sequence of maps)
    badval: float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol: float, optional
        The relative tolerance
    atol: float, optional
        The absolute tolerance

    Returns
    -------
    a bool array with the same shape as the input map, ``True`` where input map is
    close to badval, and ``False`` elsewhere.

    See Also
    --------
    :func:`mask_good`, :func:`ma`

    Examples
    --------
    >>> m = np.arange(12.)
    >>> m[3] = hpy.UNSEEN
    >>> hpy.mask_bad(m)
    array([False, False, False,  True, False, False, False, False, False,
           False, False, False], dtype=bool)
    """
    atol = npy.absolute(atol)
    rtol = npy.absolute(rtol)
    return npy.absolute(m - badval) <= atol + rtol * npy.absolute(badval)

def mask_good(m, badval = UNSEEN, rtol = 1.e-5, atol = 1.e-8):
    """Returns a bool array with ``False`` where m is close to badval.

    Parameters
    ----------
    m : a map (may be a sequence of maps)
    badval: float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol: float, optional
        The relative tolerance
    atol: float, optional
        The absolute tolerance

    Returns
    -------
    a bool array with the same shape as the input map, ``False`` where input map is
    close to badval, and ``True`` elsewhere.

    See Also
    --------
    :func:`mask_bad`, :func:`ma`

    Examples
    --------
    >>> m = np.arange(12.)
    >>> m[3] = hpy.UNSEEN
    >>> hpy.mask_good(m)
    array([ True,  True,  True, False,  True,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)
    """
    atol = npy.absolute(atol)
    rtol = npy.absolute(rtol)
    return npy.absolute(m - badval) > atol + rtol * npy.absolute(badval)

def ma(m, badval = UNSEEN, rtol = 1e-5, atol = 1e-8, copy = True):
    """Return map as a masked array, with ``badval`` pixels masked.

    Parameters
    ----------
    m : a map (may be a sequence of maps)
    badval: float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol: float, optional
        The relative tolerance
    atol: float, optional
        The absolute tolerance
    copy: bool, optional
        If ``True``, a copy of the input map is made.

    Returns
    -------
    a masked array with the same shape as the input map, masked where input map is
    close to badval.

    See Also
    --------
    :func:`mask_good`, :func:`mask_bad`, :func:`numpy.ma.masked_values`

    Examples
    --------
    >>> m = np.arange(12.)
    >>> m[3] = hpy.UNSEEN
    >>> hpy.ma(m)
    masked_array(data = [0.0 1.0 2.0 -- 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0],
                 mask = [False False False  True False False False False False False False False],
           fill_value = -1.6375e+30)
    <BLANKLINE>
    """
    return npy.ma.masked_values(m, badval, rtol = rtol, atol = atol, copy = copy)

def ang2pix(nside,theta,phi,nest=False):
    """ang2pix : nside,theta[rad],phi[rad],nest=False -> ipix (default:RING)

    Parameters
    ----------
    nside: int
      The healpix nside parameter, must be a power of 2
    theta, phi: array-like
      Angular coordinates of a point on the sphere
    nest: bool, optional
      if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    Returns
    -------
    pix: int or array of int
      The healpix pixel numbers. Scalar if all input are scalar, array otherwise.
      Usual numpy broadcasting rules apply.

    Examples
    --------
    >>> hpy.ang2pix(16, np.pi/2, 0)
    1440
    >>> hpy.ang2pix(16, [np.pi/2, np.pi/4, np.pi/2, 0, np.pi],
                    [0., np.pi/4, np.pi/2, 0, 0])
    array([1440,  427, 1520,    0, 3068])
    """
    if nest:
        return pixlib._ang2pix_nest(nside,theta,phi)
    else:
        return pixlib._ang2pix_ring(nside,theta,phi)

def pix2ang(nside,ipix,nest=False):
    """pix2ang : nside,ipix,nest=False -> theta[rad],phi[rad] (default RING)
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
        raise ValueError("Given number is not a valid nside parameter "
                         "(must be a power of 2)")
    return 12*nside**2

def nside2resol(nside, arcmin=False):
    """Give approximate resolution for nside, resolution is just the square root of the pixel area, which is a gross approximation given the different pixel shapes

    Input:
      - nside: nside paramater
    Return:
      - resol: approximate pixel size in radians (or arcmin)
    Raise a ValueError exception if nside is not valid.
    """
    if not isnsideok(nside):
        raise ValueError("Given number is not a valid nside parameter "
                         "(must be a power of 2)")
    
    resol = npy.sqrt(nside2pixarea(nside))

    if arcmin:
        resol = npy.rad2deg(resol) * 60
        
    return resol


def nside2pixarea(nside, degrees=False):
    """Give pixel area given nside

    Input:
      - nside: nside paramater
    Return:
      - pixarea: pixel area in radians
    Raise a ValueError exception if nside is not valid.
    """
    if not isnsideok(nside):
        raise ValueError("Given number is not a valid nside parameter "
                         "(must be a power of 2)")
    
    pixarea = 4*npy.pi/nside2npix(nside)

    if degrees:
        pixarea = npy.rad2deg(npy.rad2deg(pixarea))
        
    return pixarea

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
        raise ValueError("Wrong pixel number (it is not 12*nside**2)")
    nside=int(npy.floor(nside))
    if not isnsideok(nside):
        raise ValueError("Wrong nside value (it is not 2**N)")
    return nside

def isnsideok(nside):
    """Return ``True`` if nside is a valid nside parameter, False otherwise.
    Accept sequence as input, in this case return a bool array.
    """
    if hasattr(nside, '__len__'):
        return nside == 2**npy.int32(npy.around(npy.ma.log2(nside).filled(0)))
    elif nside <= 0:
        return False
    else:
        return nside == 2**int(round(npy.log2(nside)))

def isnpixok(npix):
    """Return ``True`` if npix is a valid value for healpix map size, False otherwise.
    Accept sequence as input, in this case return a bool array.
    """
    if hasattr(npix,'__len__'):
        nside = npy.sqrt(npy.asarray(npix)/12.)
        return isnsideok(nside)
    else:
        nside = npy.sqrt(npix/12.)
        return isnsideok(nside)

def get_map_size(map):
    """Try to figure out the size of the given map :
     - if map is a dict type (explicit pixel) : use nside key if present, or
       use nside attribute if present, otherwise use the smallest valid
       npix given the maximum key value
     - otherwise, return len(map)
    """
    if isinstance(map, dict):
        if 'nside' in map:
            return nside2npix(map['nside'])
        elif hasattr(map, 'nside'):
            return nside2npix(map.nside)
        else:
            nside = get_min_valid_nside(max(map.keys())+1)
            return nside2npix(nside)
    else:
        return len(map)

def get_min_valid_nside(npix):
    """Return the minimum acceptable nside so that npix <= nside2npix(nside)
    """
    order = 0.5 * npy.log2(npix / 12.)
    return 2**int(npy.ceil(order))

def get_interp_val(m,theta,phi,nest=False):
    """Return the bi-linear interpolation value of a map at given direction.

    Input:
      - m: a map (an ndarray)
      - theta, phi : the direction [rad] (either scalar or arrays of same size)
    Parameters:
      - nest: if ``True``, the map is in NEST scheme. Default: False (ie RING)
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
                    if phi is given, theta[rad],phi[rad] is a direction
    Parameters:
      - nest: if ``True``, NEST scheme. Default: False (RING)
    Return:
      - tuple of pixels and weights
    """
    if not isnsideok(nside):
        raise ValueError('Wrong nside value. Must be a power of 2.')
    if phi == None:
        theta,phi = pix2ang(nside,theta,nest=nest)
    if nest:
        r=pixlib._get_interpol_nest(nside,theta,phi)
    else:
        r=pixlib._get_interpol_ring(nside,theta,phi)
    p=npy.array(r[0:4])
    w=npy.array(r[4:8])
    return (p,w)

def get_all_neighbours(nside, theta, phi=None, nest=False):
    """Return the 8 nearest pixels
    Input:
      - nside: the nside to work with
      - theta, phi: if phi is not given, theta is actually a pixel number
                    if phi is given, theta[rad],phi[rad] is a direction

    Parameters:
      - nest: if ``True``, NEST scheme. Default: False (RING)
    Return:
      - tuple of pixels
    """
    if not isnsideok(nside):
        raise ValueError('Wrong nside value. Must be a power of 2.')
    if not (phi is None):
        theta = ang2pix(nside,theta, phi,nest=nest)
    if nest:
        r=pixlib._get_neighbors_nest(nside,theta)
    else:
        r=pixlib._get_neighbors_ring(nside,theta)
    res=npy.array(r[0:8])
    return res

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
    typ = maptype(map_in)
    if typ < 0:
        raise TypeError('map_in is not a map nor a sequence of maps')
    if typ == 0:
        npix = len(map_in)
    else:
        npix = len(map_in[0])
    nside = npix2nside(npix)
    if nside>128:
        bunchsize = npix/24
    else:
        bunchsize = npix
    if r2n:
        inp='RING'
        out='NEST'
    if n2r:
        inp='NEST'
        out='RING'
    inp = str(inp).upper()[0:4]
    out = str(out).upper()[0:4]
    if inp not in ['RING','NEST'] or out not in ['RING','NEST']:
        raise ValueError('inp and out must be either RING or NEST')
    if typ == 0:
        mapin = [map_in]
    else:
        mapin = map_in
    mapout = []
    for m_in in mapin:
        if inp == out:
            mapout.append(m_in)
        elif inp == 'RING':
            m_out = npy.zeros(npix,dtype=type(m_in[0]))
            for ibunch in range(npix/bunchsize):
                ipix_n = npy.arange(ibunch*bunchsize,
                                    (ibunch+1)*bunchsize)
                ipix_r = nest2ring(nside, ipix_n)
                m_out[ipix_n] = m_in[ipix_r]
            mapout.append(m_out)
        elif inp == 'NEST':
            m_out = npy.zeros(npix,dtype=type(m_in[0]))
            for ibunch in range(npix/bunchsize):
                ipix_r = npy.arange(ibunch*bunchsize,
                                    (ibunch+1)*bunchsize)
                ipix_n = ring2nest(nside, ipix_r)
                m_out[ipix_r] = m_in[ipix_n]
            mapout.append(m_out)
    if typ == 0:
        return mapout[0]
    else:
        return mapout

def fit_dipole(m,nest=False,bad=pixlib.UNSEEN,gal_cut=0):
    """Fit a dipole and a monopole to the map, excluding unseen pixels.
    Input:
      - m: the map from which a dipole is fitted and subtracted
      - nest: False if m is in RING scheme, ``True`` if it is NESTED
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
        ipix = ipix[(m.flat[ipix]!=bad) & (npy.isfinite(m.flat[ipix]))]
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
      - nest: False if m is in RING scheme, ``True`` if it is NESTED
      - bad: bad values of pixel, default to UNSEEN.
      - gal_cut: latitude below which pixel are not taken into account
      - fitval: whether to return or not the fitted values of monopole and dipole
    Return:
      if fitval is False:
      - the map with monopole and dipole subtracted
      if fitval is ``True``:
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
        ipix = ipix[(m.flat[ipix]!=bad) & (npy.isfinite(m.flat[ipix]))]
        x,y,z = pix2vec(nside, ipix, nest)
        m.flat[ipix] -= (dipole[0]*x)
        m.flat[ipix] -= dipole[1]*y
        m.flat[ipix] -= dipole[2]*z
        m.flat[ipix] -= mono
    if verbose:
        import rotator as R
        lon,lat = R.vec2ang(dipole,lonlat=True)
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
      - m: the map from which a monopole is fitted and subtracted
      - nest: if ``True``, input map is NESTED (default False)
      - bad: bad pixel value (default UNSEEN)
      - gal_cut: latitude cut in degrees (default 0.0)
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
    aa = v = 0.0
    for ibunch in range(npix/bunchsize):
        ipix = npy.arange(ibunch*bunchsize,
                          (ibunch+1)*bunchsize)
        ipix = ipix[(m.flat[ipix]!=bad) & (npy.isfinite(m.flat[ipix]))]
        x,y,z = pix2vec(nside, ipix, nest)
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
      - nest: if ``True``, input map is NESTED (default False)
      - bad: bad pixel value (default UNSEEN)
      - gal_cut: latitude cut in degrees (default 0.0)
      - fitval: return fitted monopole with map (default False)
      - copy: if False, input map is changed directly (default ``True``)
      - verbose: print the fitted monopole (default False)
    Return:
      if fitval is False:
      - the map with monopole subtracted
      if fitval is ``True``:
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
        ipix = ipix[(m.flat[ipix]!=bad) & (npy.isfinite(m.flat[ipix]))]
        x,y,z = pix2vec(nside, ipix, nest)
        m.flat[ipix] -= mono
    if verbose:
        print 'monopole: %g'%mono
    if fitval:
        return m,mono
    else:
        return m

def maptype(m):
    """Return -1 if the given object is not a map.
    Return 0 if it is a map.
    Return k>0 if it is a list of map (k: number of maps in the sequence)
    """
    if not hasattr(m, '__len__'):
        return -1
    if len(m) == 0:
        return -1
    if hasattr(m[0], '__len__'):
        npix=len(m[0])
        for mm in m[1:]:
            if len(mm) != npix:
                return -1
        if isnpixok(len(m[0])):
            return len(m)
    else:
        if isnpixok(len(m)):
            return 0
        else:
            return -1

def get_nside(m):
    """Return the nside of the given map.
    Can be a single map or a list of maps.
    """
    typ = maptype(m)
    if typ < 0:
        raise TypeError('m is not a map nor a sequence of maps of same size')
    if typ == 0:
        return npix2nside(len(m))
    else:
        return npix2nside(len(m[0]))

def _ud_grade_core(m,nside_out,pess=False,power=None, dtype=None):
    """Internal routine used by ud_grade. It assumes that the map is NESTED
    and single (not a list of maps)
    """
    nside_in = get_nside(m)
    if dtype:
        type_out = dtype
    else:
        type_out = type(m[0])
    if not isnsideok(nside_out):
        raise ValueError('invalid nside_out value')
    npix_in = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    if power:
        power = float(power)
        ratio = (float(nside_out)/float(nside_in))**power
    else:
        ratio = 1
    
    if nside_out > nside_in:
        rat2 = npix_out/npix_in
        fact = npy.ones(rat2, dtype=type_out)*ratio
        map_out = npy.outer(m,fact).reshape(npix_out)
    elif nside_out < nside_in:
        rat2 = npix_in/npix_out
        bads = (mask_bad(m) | (~npy.isfinite(m)))
        hit = npy.ones(npix_in,dtype=npy.int16)
        hit[bads] = 0
        m[bads] = 0
        mr = m.reshape(npix_out,rat2)
        hit = hit.reshape(npix_out,rat2)
        map_out = mr.sum(axis=1).astype(type_out)
        nhit = hit.sum(axis=1)
        if pess:
            badout = npy.where(nhit != rat2)
        else:
            badout = npy.where(nhit == 0)
        if power: nhit /= ratio
        map_out /= nhit
        try:
            map_out[badout] = UNSEEN
            m[bads] = UNSEEN
        except OverflowError:
            pass
    else:
        map_out = m
    return map_out.astype(type_out)

def ud_grade(map_in,nside_out,pess=False,order_in='RING',order_out=None,
             power=None, dtype=None):
    """Upgrade or degrade resolution of a map (or list of maps).

    Input:
     - map_in: the input map(s)
     - nside_out: the desired nside of the output
    Parameters:
     - pess: if ``True``, pessismistic, in degrading, reject pixels which contains
             a bad sub_pixel. Otherwise, estimate average with other pixels
    Output:
     - the upgraded or degraded map(s)
    """
    if not isnsideok(nside_out):
        raise ValueError('Invalid nside for output')
    typ = maptype(map_in)
    if typ<0:
        raise TypeError('Invalid map')
    if typ == 0:
        m_in = [map_in]
    else:
        m_in = map_in
    mapout = []
    if order_out is None: order_out = order_in
    for m in m_in:
        if str(order_in).upper()[0:4] == 'RING':
            m = reorder(m,r2n=True)
        mout = _ud_grade_core(m,nside_out,pess=pess, dtype=dtype)
        if str(order_out).upper()[0:4] == 'RING':
            mout = reorder(mout,n2r=True)
        mapout.append(mout)
    if typ == 0:
        return mapout[0]
    else:
        return mapout

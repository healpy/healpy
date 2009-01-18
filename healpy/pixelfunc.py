import numpy as npy
import _healpy_pixel_lib as pixlib
from _healpy_pixel_lib import UNSEEN

def ang2pix(nside,theta,phi,nest=False):
    """ang2pix : nside,theta[rad],phi[rad],nest=False -> ipix (default:RING)
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
      - theta, phi : the direction [rad] (either scalar or arrays of same size)
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
                    if phi is given, theta[rad],phi[rad] is a direction
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

def _ud_grade_core(m,nside_out,pess=False,power=None):
    """Internal routine used by ud_grade. It assumes that the map is NESTED
    and single (not a list of maps)
    """
    nside_in = get_nside(m)
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
        fact = ones(rat2,dtype=type(m[0]))*ratio
        map_out = outer(m,fact).reshape(npix_out)
    elif nside_out < nside_in:
        try:
            bad_data_val = type(m[0])(UNSEEN)
        except OverflowError:
            bad_data_present = False
        else:
            bad_data_present = True
        rat2 = npix_in/npix_out
        bads = npy.where(m==UNSEEN)
        hit = npy.ones(npix_in,dtype=npy.int16)
        hit[bads] = 0
        m[bads] = 0
        mr = m.reshape(npix_out,rat2)
        hit = hit.reshape(npix_out,rat2)
        map_out = mr.sum(axis=1)
        nhit = hit.sum(axis=1)
        if pess:
            badout = npy.where(nhit != rat2)
        else:
            badout = npy.where(nhit == 0)
        if power: nhit /= ratio
        map_out /= nhit
        if bad_data_present:
            map_out[badout] = UNSEEN
            m[bads] = UNSEEN
    else:
        map_out = m
    return map_out

def ud_grade(map_in,nside_out,pess=False,order_in='RING',order_out=None,
             power=None):
    """Upgrade or degrade resolution of a map (or list of maps).

    Input:
     - map_in: the input map(s)
     - nside_out: the desired nside of the output
    Parameters:
     - pess: if True, pessismistic, in degrading, reject pixels which contains
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
        mout = _ud_grade_core(m,nside_out,pess=pess)
        if str(order_out).upper()[0:4] == 'RING':
            mout = reorder(mout,n2r=True)
        mapout.append(mout)
    if typ == 0:
        return mapout[0]
    else:
        return mapout

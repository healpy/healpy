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
"""Provides input and output functions for Healpix maps, alm, and cl.
"""

import pyfits as pf
import numpy as np
import pixelfunc
from sphtfunc import Alm
import warnings
from pixelfunc import UNSEEN

class HealpixFitsWarning(Warning):
    pass

def read_cl(filename, dtype=np.float64, h=False):
    """Reads Cl from an healpix file, as IDL fits2cl.

    Parameters
    ----------
    filename : str
      the fits file name
    dtype : data type, optional
      the data type of the returned array

    Returns
    -------
    cl : array
      the cl array
    """
    hdulist=pf.open(filename)
    cl = [hdulist[1].data.field(n) for n in range(len(hdulist[1].data.columns))]
    hdulist.close()
    if len(cl) == 1:
        return cl[0]
    else:
        return cl

def write_cl(filename, cl, dtype=np.float64):
    """Writes Cl into an healpix file, as IDL cl2fits.

    Parameters
    ----------
    filename : str
      the fits file name
    cl : array
      the cl array to write to file, currently TT only
    """
    # check the dtype and convert it
    fitsformat = getformat(dtype)
    column_names = ['TEMPERATURE','GRADIENT','CURL','G-T','C-T','C-G']
    if isinstance(cl, list):
        cols = [pf.Column(name=column_name,
                               format='%s'%fitsformat,
                               array=column_cl) for column_name, column_cl in zip(column_names[:len(cl)], cl)]
    else: # we write only one TT
        cols = [pf.Column(name='TEMPERATURE',
                               format='%s'%fitsformat,
                               array=cl)]
            
    tbhdu = pf.new_table(cols)
    # add needed keywords
    tbhdu.header.update('CREATOR','healpy')
    tbhdu.writeto(filename,clobber=True)

def write_map(filename,m,nest=False,dtype=np.float32,fits_IDL=True,coord=None):
    """Writes an healpix map into an healpix file.

    Parameters
    ----------
    filename : str
      the fits file name
    m : array or sequence of 3 arrays
      the map to write. Possibly a sequence of 3 maps of same size.
      They will be considered as I, Q, U maps. 
      Supports masked maps, see the `ma` function.
    nest : bool, optional
      If False, ordering scheme is NESTED, otherwise, it is RING. Default: RING.
    fits_IDL : bool, optional
      If True, reshapes columns in rows of 1024, otherwise all the data will 
      go in one column. Default: True
    coord : str
      The coordinate system, typically 'E' for Ecliptic, 'G' for Galactic or 'Q' for Equatorial  
    """
    if not hasattr(m, '__len__'):
        raise TypeError('The map must be a sequence')
    # check the dtype and convert it
    fitsformat = getformat(dtype)
    #print 'format to use: "%s"'%fitsformat
    if hasattr(m[0], '__len__'):
        # we should have three maps
        if len(m) != 3 or len(m[1]) != len(m[0]) or len(m[2]) != len(m[0]):
            raise ValueError("You should give 3 maps of same size "
                             "for polarisation...")
        nside = pixelfunc.npix2nside(len(m[0]))
        if nside < 0:
            raise ValueError('Invalid healpix map : wrong number of pixel')
        m = pixelfunc.ma_to_array(m)
        cols=[]
        colnames=['I_STOKES','Q_STOKES','U_STOKES']
        for cn,mm in zip(colnames,m):
            if len(mm) > 1024 and fits_IDL:
                # I need an ndarray, for reshape:
                mm2 = np.asarray(mm)
                cols.append(pf.Column(name=cn,
                                       format='1024%s'%fitsformat,
                                       array=mm2.reshape(mm2.size/1024,1024)))
            else:
                cols.append(pf.Column(name=cn,
                                       format='%s'%fitsformat,
                                       array=mm))
    else: # we write only one map
        nside = pixelfunc.npix2nside(len(m))
        if nside < 0:
            raise ValueError('Invalid healpix map : wrong number of pixel')
        if m.size > 1024 and fits_IDL:
            cols = [pf.Column(name='I_STOKES',
                               format='1024%s'%fitsformat,
                               array=m.reshape(m.size/1024,1024))]
        else:
            cols = [pf.Column(name='I_STOKES',
                               format='%s'%fitsformat,
                               array=m)]
            
    tbhdu = pf.new_table(cols)
    # add needed keywords
    tbhdu.header.update('PIXTYPE','HEALPIX','HEALPIX pixelisation')
    if nest: ordering = 'NESTED'
    else:    ordering = 'RING'
    tbhdu.header.update('ORDERING',ordering,
                        'Pixel ordering scheme, either RING or NESTED')
    if coord:
        tbhdu.header.update('COORDSYS',coord,
                            'Ecliptic, Galactic or eQuatorial')
    tbhdu.header.update('EXTNAME','xtension',
                        'name of this binary table extension')
    tbhdu.header.update('NSIDE',nside,'Resolution parameter of HEALPIX')
    tbhdu.header.update('FIRSTPIX', 0, 'First pixel # (0 based)')
    tbhdu.header.update('LASTPIX',pixelfunc.nside2npix(nside)-1,
                        'Last pixel # (0 based)')
    tbhdu.header.update('INDXSCHM','IMPLICIT',
                        'Indexing: IMPLICIT or EXPLICIT')
    tbhdu.writeto(filename,clobber=True)


def read_map(filename,field=0,dtype=np.float64,nest=False,hdu=1,h=False,
             verbose=False):
    """Read an healpix map from a fits file.

    Parameters
    ----------
    filename : str
      the fits file name
    field : int or tuple of int, optional
      The column to read. Default: 0.
      By convention 0 is temperature, 1 is Q, 2 is U.
      Field can be a tuple to read multiple columns (0,1,2)
    dtype : data type, optional
      Force the conversion to some type. Default: np.float64
    nest : bool, optional
      If True return the map in NEST ordering, otherwise in RING ordering; 
      use fits keyword ORDERING to decide whether conversion is needed or not
      If None, no conversion is performed.
    hdu : int, optional
      the header number to look at (start at 0)
    h : bool, optional
      If True, return also the header. Default: False.
    verbose : bool, optional
      If True, print a number of diagnostic messages
    as_ma : bool, optional
      If True, return also the header. Default: False.

    Returns
    -------
    m | (m0, m1, ...) [, header] : array or a tuple of arrays, optionally with header appended
      The map(s) read from the file, and the header if *h* is True.
    """
    hdulist=pf.open(filename)
    #print hdulist[1].header
    nside = hdulist[hdu].header.get('NSIDE')
    if nside is None:
        warnings.warn("No NSIDE in the header file : will use length of array", HealpixFitsWarning)
    else:
        nside = int(nside)
    if verbose: print 'NSIDE = %d'%nside

    if not pixelfunc.isnsideok(nside):
        raise ValueError('Wrong nside parameter.')
    ordering = hdulist[hdu].header.get('ORDERING','UNDEF').strip()
    if ordering == 'UNDEF':
        ordering = (nest and 'NESTED' or 'RING')
        warnings.warn("No ORDERING keyword in header file : "
                      "assume %s"%ordering)
    if verbose: print 'ORDERING = %s in fits file'%ordering

    sz=pixelfunc.nside2npix(nside)
    if not (hasattr(field, '__len__') or isinstance(field, str)):
        field = (field,)
    ret = []

    for ff in field:
        try:
            m=hdulist[hdu].data.field(ff).astype(dtype).ravel()
        except pf.VerifyError, e:
            print(e)
            print("Trying to fix a badly formatted header")
            m=hdulist[hdu].verify("fix")
            m=hdulist[hdu].data.field(ff).astype(dtype).ravel()
            
        if (not pixelfunc.isnpixok(m.size) or (sz>0 and sz != m.size)) and verbose:
            print 'nside=%d, sz=%d, m.size=%d'%(nside,sz,m.size)
            raise ValueError('Wrong nside parameter.')
        if not nest is None: # no conversion with None
            if nest and ordering == 'RING':
                idx = pixelfunc.nest2ring(nside,np.arange(m.size,dtype=np.int32))
                m = m[idx]
                if verbose: print 'Ordering converted to NEST'
            elif (not nest) and ordering == 'NESTED':
                idx = pixelfunc.ring2nest(nside,np.arange(m.size,dtype=np.int32))
                m = m[idx]
                if verbose: print 'Ordering converted to RING'
        try:
            m[pixelfunc.mask_bad(m)] = UNSEEN
        except OverflowError:
            pass
        ret.append(m)
    
    if len(ret) == 1:
        if h:
            return ret[0],hdulist[hdu].header.items()
        else:
            return ret[0]
    else:
        if h:
            ret.append(hdulist[hdu].header.items())
            return tuple(ret)
        else:
            return tuple(ret)


def write_alm(filename,alms,out_dtype=None,lmax=-1,mmax=-1,mmax_in=-1):
    """Write alms to a fits file. 
    
    In the fits file the alms are written 
    with explicit index scheme, index = l*l + l + m +1, possibly out of order.
    By default write_alm makes a table with the same precision as the alms.
    If specified, the lmax and mmax parameters truncate the input data to
    include only alms for which l <= lmax and m <= mmax.

    Parameters
    ----------
    filename : str
      The filename of the output fits file
    alms : array, complex
      A complex ndarray holding the alms.
    lmax : int, optional
      The maximum l in the output file
    mmax : int, optional
      The maximum m in the output file
    out_dtype : data type, optional
      data type in the output file (must be a numpy dtype). Default: *alms*.real.dtype
    mmax_in : int, optional
      maximum m in the input array
    """

    l2max = Alm.getlmax(len(alms),mmax=mmax_in)
    if (lmax != -1 and lmax > l2max):
        raise ValueError("Too big lmax in parameter")
    elif lmax == -1:
        lmax = l2max

    if mmax_in == -1:
	mmax_in = l2max

    if mmax == -1:
        mmax = lmax
    if mmax > mmax_in:
	mmax = mmax_in

    if (out_dtype == None):
        out_dtype = alms.real.dtype

    l,m = Alm.getlm(lmax)
    idx = np.where((l <= lmax)*(m <= mmax))
    l = l[idx]
    m = m[idx]

    idx_in_original = Alm.getidx(l2max, l=l, m=m)
    
    index = l**2 + l + m + 1

    out_data = np.empty(len(index),
               dtype=[('index','i'),
                      ('real',out_dtype),
                      ('imag',out_dtype)])
    out_data['index'] = index
    out_data['real'] = alms.real[idx_in_original]
    out_data['imag'] = alms.imag[idx_in_original]

    cindex = pf.Column(name="index", format=getformat(np.int32), unit="l*l+l+m+1", array=out_data['index'])
    creal = pf.Column(name="real", format=getformat(out_dtype), unit="unknown", array=out_data['real'])
    cimag = pf.Column(name="imag", format=getformat(out_dtype), unit="unknown", array=out_data['imag'])

    tbhdu = pf.new_table([cindex,creal,cimag])
    tbhdu.writeto(filename,clobber=True)       
    
def read_alm(filename,hdu=1,return_mmax=False):
    """Read alm from a fits file. 
    
    In the fits file, the alm are written
    with explicit index scheme, index = l**2+l+m+1, while healpix cxx
    uses index = m*(2*lmax+1-m)/2+l. The conversion is done in this 
    function.

    Parameters
    ----------
    filename : str
      The name of the fits file to read
    hdu : int, optional
      The header to read. Start at 0. Default: hdu=1
    return_mmax : bool, optional
      If true, both the alms and mmax is returned in a tuple. Default: return_mmax=False

    Returns
    -------
    alms[, mmax] : complex array or tuple of a complex array and an int
      The alms read from the file and optionally mmax read from the file
    """
    idx, almr, almi = mrdfits(filename,hdu=hdu)
    l = np.floor(np.sqrt(idx-1)).astype(long)
    m = idx - l**2 - l - 1
    if (m<0).any():
        raise ValueError('Negative m value encountered !')
    lmax = l.max()
    mmax = m.max()
    alm = almr*(0+0j)
    i = Alm.getidx(lmax,l,m)
    alm.real[i] = almr
    alm.imag[i] = almi
    if return_mmax:
        return alm, mmax
    else:
        return alm

## Generic functions to read and write column of data in fits file

def mrdfits(filename,hdu=1):
    """Read a table in a fits file.

    Parameters
    ----------
    filename : str
      The name of the fits file to read
    hdu : int, optional
      The header to read. Start at 0. Default: hdu=1

    Returns
    -------
    cols : a list of arrays
      A list of column data in the given header
    """
    hdulist=pf.open(filename)
    if hdu>=len(hdulist):
        raise ValueError('Available hdu in [0-%d]'%len(hdulist))
    hdu=hdulist[hdu]
    val=[]
    for i in range(len(hdu.columns)):
        val.append(hdu.data.field(i))
    hdulist.close()
    del hdulist
    return val

def mwrfits(filename,data,hdu=1,colnames=None,keys=None):
    """Write columns to a fits file in a table extension.

    Parameters
    ----------
    filename : str
      The fits file name
    data : list of 1D arrays
      A list of 1D arrays to write in the table
    hdu : int, optional
      The header where to write the data. Default: 1
    colnames : list of str
      The column names
    keys : dict-like
      A dictionary with keywords to write in the header
    """
    # Check the inputs
    if colnames is not None:
        if len(colnames) != len(data):
            raise ValueError("colnames and data must the same length")
    else:
        colnames = ['']*len(data)
    cols=[]
    for line in xrange(len(data)):
        cols.append(pf.Column(name=colnames[line],
                               format=getformat(data[line]),
                               array=data[line]))
    tbhdu = pf.new_table(cols)
    if type(keys) is dict:
        for k,v in keys.items():
            tbhdu.header.update(k,v)
    # write the file
    tbhdu.writeto(filename,clobber=True)

def getformat(t):
    """Get the FITS convention format string of data type t.

    Parameters
    ----------
    t : data type
      The data type for which the FITS type is requested

    Returns
    -------
    fits_type : str or None
      The FITS string code describing the data type, or None if unknown type.
    """
    conv = {
        np.dtype(np.bool): 'L',
        np.dtype(np.uint8): 'B',
        np.dtype(np.int16): 'I',
        np.dtype(np.int32): 'J',
        np.dtype(np.int64): 'K',
        np.dtype(np.float32): 'E',
        np.dtype(np.float64): 'D',
        np.dtype(np.complex64): 'C',
        np.dtype(np.complex128): 'M'
        }
    try:
        if t in conv:
            return conv[t]
    except:
        pass
    try:
        if np.dtype(t) in conv:
            return conv[np.dtype(t)]
    except:
        pass
    try:
        if np.dtype(type(t)) in conv:
            return conv[np.dtype(type(t))]
    except:
        pass
    try:
        if np.dtype(type(t[0])) in conv:
            return conv[np.dtype(type(t[0]))]
    except:
        pass
    try:
        if t is str:
            return 'A'
    except:
        pass
    try:
        if type(t) is str:
            return 'A%d'%(len(t))
    except:
        pass
    try:
        if type(t[0]) is str:
            l=max(len(s) for s in t)
            return 'A%d'%(l)
    except:
        pass

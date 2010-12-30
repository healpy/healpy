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

import pyfits as pyf
import numpy as npy
import pixelfunc
from sphtfunc import Alm
import warnings
from _healpy_pixel_lib import UNSEEN
from exceptions import NotImplementedError

class HealpixFitsWarning(Warning):
    pass

def read_cl(filename,dtype=npy.float32,h=False):
    """Reads Cl from an healpix file, as IDL fits2cl.

    Input:
      - filename: the fits file name

    Return: 
      - cl: the cl array, currently TT only
    """
    hdulist=pyf.open(filename)
    return hdulist[1].data.field(0)

def write_cl(filename,cl,dtype=npy.float32):
    """Writes Cl into an healpix file, as IDL cl2fits.

    Input:
      - filename: the fits file name
      - cl: the cl array to write to file, currently TT only
    """
    # check the dtype and convert it
    fitsformat = getformat(dtype)
    if isinstance(cl, list):
        raise NotImplementedError('Currently it supports only temperature-only cls')
    else: # we write only one TT
        cols = [pyf.Column(name='TEMPERATURE',
                               format='%s'%fitsformat,
                               array=cl)]
            
    coldefs=pyf.ColDefs(cols)
    tbhdu = pyf.new_table(coldefs)
    # add needed keywords
    tbhdu.header.update('CREATOR','healpy')
    tbhdu.writeto(filename,clobber=True)

def write_map(filename,m,nest=False,dtype=npy.float32,fits_IDL=True):
    """Writes an healpix map into an healpix file.

    Input:
      - filename: the fits file name
      - m: the map to write. Possibly a sequence of 3 maps of same size.
           They will be considered as I, Q, U maps
      - nest=False: ordering scheme
      - fits_IDL = true reshapes columns in rows of 1024, otherwise all the data will 
        go in one column
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
        cols=[]
        colnames=['I_STOKES','Q_STOKES','U_STOKES']
        for cn,mm in zip(colnames,m):
            if len(mm) > 1024 and fits_IDL:
                # I need an ndarray, for reshape:
                mm2 = npy.asarray(mm)
                cols.append(pyf.Column(name=cn,
                                       format='1024%s'%fitsformat,
                                       array=mm2.reshape(mm2.size/1024,1024)))
            else:
                cols.append(pyf.Column(name=cn,
                                       format='%s'%fitsformat,
                                       array=mm))
    else: # we write only one map
        nside = pixelfunc.npix2nside(len(m))
        if nside < 0:
            raise ValueError('Invalid healpix map : wrong number of pixel')
        if m.size > 1024 and fits_IDL:
            cols = [pyf.Column(name='I_STOKES',
                               format='1024%s'%fitsformat,
                               array=m.reshape(m.size/1024,1024))]
        else:
            cols = [pyf.Column(name='I_STOKES',
                               format='%s'%fitsformat,
                               array=m)]
            
    coldefs=pyf.ColDefs(cols)
    tbhdu = pyf.new_table(coldefs)
    # add needed keywords
    tbhdu.header.update('PIXTYPE','HEALPIX','HEALPIX pixelisation')
    if nest: ordering = 'NESTED'
    else:    ordering = 'RING'
    tbhdu.header.update('ORDERING',ordering,
                        'Pixel ordering scheme, either RING or NESTED')
    tbhdu.header.update('EXTNAME','xtension',
                        'name of this binary table extension')
    tbhdu.header.update('NSIDE',nside,'Resolution parameter of HEALPIX')
    tbhdu.header.update('FIRSTPIX', 0, 'First pixel # (0 based)')
    tbhdu.header.update('LASTPIX',pixelfunc.nside2npix(nside)-1,
                        'Last pixel # (0 based)')
    tbhdu.header.update('INDXSCHM','IMPLICIT',
                        'Indexing: IMPLICIT or EXPLICIT')
    tbhdu.writeto(filename,clobber=True)


def read_map(filename,field=0,dtype=npy.float64,nest=False,hdu=1,h=False,
             verbose=False):
    """Read an healpix map from a fits file.

    Input:
      - filename: the fits file name
    Parameters:
      - field: the column to read Default: 0
               by convention 0 is temperature, 1 is Q, 2 is U
               field can be a tuple to read multiple columns (0,1,2)
      - dtype: force the conversion to some type. Default: npy.float64
      - nest=False: if True return the map in NEST ordering, otherwise in RING ordering; 
                    use fits keyword ORDERING to decide whether conversion is needed or not
                    if None, no conversion is performed
      - hdu=1: the header number to look at (start at 0)
      - h=False: if True, return also the header
      - verbose=False: if True, print a number of diagnostic messages
    Return:
      - an array, a tuple of array, possibly with the header at the end if h
        is True
    """
    hdulist=pyf.open(filename)
    #print hdulist[1].header
    nside = int(hdulist[hdu].header.get('NSIDE'))
    if nside is None:
        warnings.warn("No NSIDE in the header file : will use length of array",
                      HealpixFitsWarning)
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
    if not hasattr(field, '__len__'):
        field = (field,)
    ret = []

    for ff in field:
        m=hdulist[hdu].data.field(ff).astype(dtype).ravel()
        if (not pixelfunc.isnpixok(m.size) or (sz>0 and sz != m.size)) and verbose:
            print 'nside=%d, sz=%d, m.size=%d'%(nside,sz,m.size)
            raise ValueError('Wrong nside parameter.')
        if nest != None: # no conversion with None
            if nest and ordering == 'RING':
                idx = pixelfunc.nest2ring(nside,npy.arange(m.size,dtype=npy.int32))
                m = m[idx]
                if verbose: print 'Ordering converted to NEST'
            elif (not nest) and ordering == 'NESTED':
                idx = pixelfunc.ring2nest(nside,npy.arange(m.size,dtype=npy.int32))
                m = m[idx]
                if verbose: print 'Ordering converted to RING'
        m[m<-1.637e30] = UNSEEN
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
    """
    Write alms to a fits file. In the fits file the alms are written 
    with explicit index scheme, index = l*l + l + m +1, possibly out of order.
    By default write_alm makes a table with the same precision as the alms.
    If specified, the lmax and mmax parameters truncate the input data to
    include only alms for which l <= lmax and m <= mmax.

    Input:
      - filename: the filename of the output fits file
      - alms: a complex array holding the alms.

    Parameters:
     - lmax: maximum l in the output file
     - mmax: maximum m in the output file
     - out_dtype: data type in the output file (must be a numpy dtype)
     - mmax_in: maximum m in the input array
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
    idx = npy.where((l <= lmax)*(m <= mmax))
    l = l[idx]
    m = m[idx]

    idx_in_original = Alm.getidx(l2max, l=l, m=m)
    
    index = l**2 + l + m + 1

    out_data = npy.empty(len(index),\
                             dtype=[('index','i'),\
                                        ('real',out_dtype),('imag',out_dtype)])
    out_data['index'] = index
    out_data['real'] = alms.real[idx_in_original]
    out_data['imag'] = alms.imag[idx_in_original]

    cindex = pyf.Column(name="index", format=getformat(npy.int32), unit="l*l+l+m+1", array=out_data['index'])
    creal = pyf.Column(name="real", format=getformat(out_dtype), unit="unknown", array=out_data['real'])
    cimag = pyf.Column(name="imag", format=getformat(out_dtype), unit="unknown", array=out_data['imag'])

    coldefs=pyf.ColDefs([cindex,creal,cimag])
    tbhdu = pyf.new_table(coldefs)
    tbhdu.writeto(filename,clobber=True)       
    
def read_alm(filename,hdu=1,return_mmax=False):
    """Read alm from a fits file. In the fits file, the alm are written
    with explicit index scheme, index = l**2+l+m+1, while healpix cxx
    uses index = m*(2*lmax+1-m)/2+l. The conversion is done in this 
    function.
    """
    idx, almr, almi = mrdfits(filename,hdu=hdu)
    l = npy.floor(npy.sqrt(idx-1)).astype(long)
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

    Input:
      - filename: the name of the fits file to read
    Parameters:
      - hdu: the header to read. Start at 0. Default: hdu=1
    Return:
      - a list of column data in the given header
    """
    hdulist=pyf.open(filename)
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

    Input:
      - filename: the fits file name
      - data: a list of 1D arrays to write in the table
    Parameters:
      - hdu: header where to write the data. Default: 1
      - colnames: the column names
      - keys: a dictionary with keywords to write in the header
    """
    # Check the inputs
    if colnames is not None:
        if len(colnames) != len(data):
            raise ValueError("colnames and data must the same length")
    else:
        colnames = ['']*len(data)
    cols=[]
    for line in xrange(len(data)):
        cols.append(pyf.Column(name=colnames[line],
                               format=getformat(data[line]),
                               array=data[line]))
    coldefs=pyf.ColDefs(cols)
    tbhdu = pyf.new_table(coldefs)
    if type(keys) is dict:
        for k,v in keys.items():
            tbhdu.header.update(k,v)
    # write the file
    tbhdu.writeto(filename,clobber=True)

def getformat(t):
    """Get the format string of type t.
    """
    conv = {
        npy.dtype(npy.bool): 'L',
        npy.dtype(npy.uint8): 'B',
        npy.dtype(npy.int16): 'I',
        npy.dtype(npy.int32): 'J',
        npy.dtype(npy.int64): 'K',
        npy.dtype(npy.float32): 'E',
        npy.dtype(npy.float64): 'D',
        npy.dtype(npy.complex64): 'C',
        npy.dtype(npy.complex128): 'M'
        }
    try:
        if t in conv:
            return conv[t]
    except:
        pass
    try:
        if npy.dtype(t) in conv:
            return conv[npy.dtype(t)]
    except:
        pass
    try:
        if npy.dtype(type(t)) in conv:
            return conv[npy.dtype(type(t))]
    except:
        pass
    try:
        if npy.dtype(type(t[0])) in conv:
            return conv[npy.dtype(type(t[0]))]
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
        

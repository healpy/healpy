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
from __future__ import division

import six
import warnings
import astropy.io.fits as pf
import numpy as np

from . import pixelfunc
from .sphtfunc import Alm
from .pixelfunc import UNSEEN
from . import cookbook as cb

standard_column_names = {
    1: "TEMPERATURE",
    2: ["Q_POLARISATION", "U_POLARISATION"],
    3: ["TEMPERATURE", "Q_POLARISATION", "U_POLARISATION"],
    6: ["II", "IQ", "IU", "QQ", "QU", "UU"],
}


class HealpixFitsWarning(Warning):
    pass


def read_cl(filename, dtype=np.float64, h=False):
    """Reads Cl from a healpix file, as IDL fits2cl.

    Parameters
    ----------
    filename : str or HDUList or HDU
      the fits file name
    dtype : data type, optional
      the data type of the returned array

    Returns
    -------
    cl : array
      the cl array
    """
    fits_hdu = _get_hdu(filename, hdu=1)
    cl = np.array([fits_hdu.data.field(n) for n in range(len(fits_hdu.columns))])
    if len(cl) == 1:
        return cl[0]
    else:
        return cl


def write_cl(filename, cl, dtype=np.float64, overwrite=False):
    """Writes Cl into a healpix file, as IDL cl2fits.

    Parameters
    ----------
    filename : str
      the fits file name
    cl : array
      the cl array to write to file
    overwrite : bool, optional
      If True, existing file is silently overwritten. Otherwise trying to write
      an existing file raises an OSError (IOError for Python 2).
    """
    # check the dtype and convert it
    fitsformat = getformat(dtype)
    column_names = ["TEMPERATURE", "GRADIENT", "CURL", "G-T", "C-T", "C-G"]
    if len(np.shape(cl)) == 2:
        cols = [
            pf.Column(name=column_name, format="%s" % fitsformat, array=column_cl)
            for column_name, column_cl in zip(column_names[: len(cl)], cl)
        ]
    elif len(np.shape(cl)) == 1:
        # we write only TT
        cols = [pf.Column(name="TEMPERATURE", format="%s" % fitsformat, array=cl)]
    else:
        raise RuntimeError("write_cl: Expected one or more vectors of equal length")

    tbhdu = pf.BinTableHDU.from_columns(cols)
    # add needed keywords
    tbhdu.header["CREATOR"] = "healpy"
    tbhdu.writeto(filename, overwrite=overwrite)


def write_map(
    filename,
    m,
    nest=False,
    dtype=np.float32,
    fits_IDL=True,
    coord=None,
    partial=False,
    column_names=None,
    column_units=None,
    extra_header=(),
    overwrite=False,
):
    """Writes a healpix map into a healpix file.

    Parameters
    ----------
    filename : str
      the fits file name
    m : array or sequence of 3 arrays
      the map to write. Possibly a sequence of 3 maps of same size.
      They will be considered as I, Q, U maps.
      Supports masked maps, see the `ma` function.
    nest : bool, optional
      If True, ordering scheme is assumed to be NESTED, otherwise, RING. Default: RING.
      The map ordering is not modified by this function, the input map array
      should already be in the desired ordering (run `ud_grade` beforehand).
    fits_IDL : bool, optional
      If True, reshapes columns in rows of 1024, otherwise all the data will
      go in one column. Default: True
    coord : str
      The coordinate system, typically 'E' for Ecliptic, 'G' for Galactic or 'C' for
      Celestial (equatorial)
    partial : bool, optional
      If True, fits file is written as a partial-sky file with explicit indexing.
      Otherwise, implicit indexing is used.  Default: False.
    column_names : str or list
      Column name or list of column names, if None here the default column names based on
      the number of columns:
      1 : "TEMPERATURE",
      2 : ["Q_POLARISATION", "U_POLARISATION"],
      3 : ["TEMPERATURE", "Q_POLARISATION", "U_POLARISATION"],
      6 : ["II", "IQ", "IU", "QQ", "QU", "UU"]
      COLUMN_1, COLUMN_2... otherwise (FITS is 1-based)
    column_units : str or list
      Units for each column, or same units for all columns.
    extra_header : list
      Extra records to add to FITS header.
    dtype: np.dtype or list of np.dtypes, optional
      The datatype in which the columns will be stored. Will be converted
      internally from the numpy datatype to the fits convention. If a list,
      the length must correspond to the number of map arrays.
      Default: np.float32.
    overwrite : bool, optional
      If True, existing file is silently overwritten. Otherwise trying to write
      an existing file raises an OSError (IOError for Python 2).
    """
    if not hasattr(m, "__len__"):
        raise TypeError("The map must be a sequence")

    m = pixelfunc.ma_to_array(m)
    if pixelfunc.maptype(m) == 0:  # a single map is converted to a list
        m = [m]

    # check the dtype and convert it
    try:
        fitsformat = []
        for curr_dtype in dtype:
            fitsformat.append(getformat(curr_dtype))
    except TypeError:
        # dtype is not iterable
        fitsformat = [getformat(dtype)] * len(m)

    if column_names is None:
        column_names = standard_column_names.get(
            len(m), ["COLUMN_%d" % n for n in range(1, len(m) + 1)]
        )
    else:
        assert len(column_names) == len(m), "Length column_names != number of maps"

    if column_units is None or isinstance(column_units, six.string_types):
        column_units = [column_units] * len(m)

    # maps must have same length
    assert len(set(map(len, m))) == 1, "Maps must have same length"
    nside = pixelfunc.npix2nside(len(m[0]))

    if nside < 0:
        raise ValueError("Invalid healpix map : wrong number of pixel")

    cols = []
    if partial:
        fits_IDL = False
        mask = pixelfunc.mask_good(m[0])
        pix = np.where(mask)[0]
        if len(pix) == 0:
            raise ValueError("Invalid healpix map : empty partial map")
        m = [mm[mask] for mm in m]
        ff = getformat(np.min_scalar_type(-pix.max()))
        if ff is None:
            ff = "I"
        cols.append(pf.Column(name="PIXEL", format=ff, array=pix, unit=None))

    for cn, cu, mm, curr_fitsformat in zip(column_names, column_units, m, fitsformat):
        if len(mm) > 1024 and fits_IDL:
            # I need an ndarray, for reshape:
            mm2 = np.asarray(mm)
            cols.append(
                pf.Column(
                    name=cn,
                    format="1024%s" % curr_fitsformat,
                    array=mm2.reshape(mm2.size // 1024, 1024),
                    unit=cu,
                )
            )
        else:
            cols.append(
                pf.Column(name=cn, format="%s" % curr_fitsformat, array=mm, unit=cu)
            )

    tbhdu = pf.BinTableHDU.from_columns(cols)
    # add needed keywords
    tbhdu.header["PIXTYPE"] = ("HEALPIX", "HEALPIX pixelisation")
    if nest:
        ordering = "NESTED"
    else:
        ordering = "RING"
    tbhdu.header["ORDERING"] = (
        ordering,
        "Pixel ordering scheme, either RING or NESTED",
    )
    if coord:
        tbhdu.header["COORDSYS"] = (
            coord,
            "Ecliptic, Galactic or Celestial (equatorial)",
        )
    tbhdu.header["EXTNAME"] = ("xtension", "name of this binary table extension")
    tbhdu.header["NSIDE"] = (nside, "Resolution parameter of HEALPIX")
    if not partial:
        tbhdu.header["FIRSTPIX"] = (0, "First pixel # (0 based)")
        tbhdu.header["LASTPIX"] = (
            pixelfunc.nside2npix(nside) - 1,
            "Last pixel # (0 based)",
        )
    tbhdu.header["INDXSCHM"] = (
        "EXPLICIT" if partial else "IMPLICIT",
        "Indexing: IMPLICIT or EXPLICIT",
    )
    tbhdu.header["OBJECT"] = (
        "PARTIAL" if partial else "FULLSKY",
        "Sky coverage, either FULLSKY or PARTIAL",
    )

    # FIXME: In modern versions of Pyfits, header.update() understands a
    # header as an argument, and headers can be concatenated with the `+'
    # operator.
    for args in extra_header:
        tbhdu.header[args[0]] = args[1:]

    tbhdu.writeto(filename, overwrite=overwrite)


def read_map(
    filename,
    field=0,
    dtype=np.float64,
    nest=False,
    partial=False,
    hdu=1,
    h=False,
    verbose=True,
    memmap=False,
):
    """Read a healpix map from a fits file.  Partial-sky files,
    if properly identified, are expanded to full size and filled with UNSEEN.

    Parameters
    ----------
    filename : str or HDU or HDUList
      the fits file name
    field : int or tuple of int, or None, optional
      The column to read. Default: 0.
      By convention 0 is temperature, 1 is Q, 2 is U.
      Field can be a tuple to read multiple columns (0,1,2)
      If the fits file is a partial-sky file, field=0 corresponds to the
      first column after the pixel index column.
      If None, all columns are read in.
    dtype : data type or list of data types, optional
      Force the conversion to some type. Passing a list allows different
      types for each field. In that case, the length of the list must
      correspond to the length of the field parameter. Default: np.float64
    nest : bool, optional
      If True return the map in NEST ordering, otherwise in RING ordering;
      use fits keyword ORDERING to decide whether conversion is needed or not
      If None, no conversion is performed.
    partial : bool, optional
      If True, fits file is assumed to be a partial-sky file with explicit indexing,
      if the indexing scheme cannot be determined from the header.
      If False, implicit indexing is assumed.  Default: False.
      A partial sky file is one in which OBJECT=PARTIAL and INDXSCHM=EXPLICIT,
      and the first column is then assumed to contain pixel indices.
      A full sky file is one in which OBJECT=FULLSKY and INDXSCHM=IMPLICIT.
      At least one of these keywords must be set for the indexing
      scheme to be properly identified.
    hdu : int, optional
      the header number to look at (start at 0)
    h : bool, optional
      If True, return also the header. Default: False.
    verbose : bool, optional
      If True, print a number of diagnostic messages
    memmap : bool, optional
      Argument passed to astropy.io.fits.open, if True, the map is not read into memory,
      but only the required pixels are read when needed. Default: False.

    Returns
    -------
    m | (m0, m1, ...) [, header] : array or a tuple of arrays, optionally with header appended
      The map(s) read from the file, and the header if *h* is True.
    """

    fits_hdu = _get_hdu(filename, hdu=hdu, memmap=memmap)

    nside = fits_hdu.header.get("NSIDE")
    if nside is None:
        warnings.warn(
            "No NSIDE in the header file : will use length of array", HealpixFitsWarning
        )
    else:
        nside = int(nside)
    if verbose:
        print("NSIDE = {0:d}".format(nside))

    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside parameter.")
    ordering = fits_hdu.header.get("ORDERING", "UNDEF").strip()
    if ordering == "UNDEF":
        ordering = nest and "NESTED" or "RING"
        warnings.warn("No ORDERING keyword in header file : " "assume %s" % ordering)
    if verbose:
        print("ORDERING = {0:s} in fits file".format(ordering))

    sz = pixelfunc.nside2npix(nside)
    ret = []

    # partial sky: check OBJECT, then INDXSCHM
    obj = fits_hdu.header.get("OBJECT", "UNDEF").strip()
    if obj != "UNDEF":
        if obj == "PARTIAL":
            partial = True
        elif obj == "FULLSKY":
            partial = False

    schm = fits_hdu.header.get("INDXSCHM", "UNDEF").strip()
    if schm != "UNDEF":
        if schm == "EXPLICIT":
            if obj == "FULLSKY":
                raise ValueError("Incompatible INDXSCHM keyword")
            partial = True
        elif schm == "IMPLICIT":
            if obj == "PARTIAL":
                raise ValueError("Incompatible INDXSCHM keyword")
            partial = False

    if schm == "UNDEF":
        schm = partial and "EXPLICIT" or "IMPLICIT"
        warnings.warn("No INDXSCHM keyword in header file : " "assume {}".format(schm))
    if verbose:
        print("INDXSCHM = {0:s}".format(schm))

    if field is None:
        field = range(len(fits_hdu.data.columns) - 1 * partial)
    if not (hasattr(field, "__len__") or isinstance(field, str)):
        field = (field,)

    if partial:
        # increment field counters
        field = tuple(f if isinstance(f, str) else f + 1 for f in field)
        try:
            pix = fits_hdu.data.field(0).astype(int, copy=False).ravel()
        except pf.VerifyError as e:
            print(e)
            print("Trying to fix a badly formatted header")
            fits_hdu.verify("fix")
            pix = fits_hdu.data.field(0).astype(int, copy=False).ravel()

    try:
        assert len(dtype) == len(
            field
        ), "The number of dtypes are not equal to the number of fields"
    except TypeError:
        dtype = [dtype] * len(field)

    for ff, curr_dtype in zip(field, dtype):
        try:
            m = fits_hdu.data.field(ff).astype(curr_dtype, copy=False).ravel()
        except pf.VerifyError as e:
            print(e)
            print("Trying to fix a badly formatted header")
            m = fits_hdu.verify("fix")
            m = fits_hdu.data.field(ff).astype(curr_dtype, copy=False).ravel()

        if partial:
            mnew = UNSEEN * np.ones(sz, dtype=curr_dtype)
            mnew[pix] = m
            m = mnew

        if (not pixelfunc.isnpixok(m.size) or (sz > 0 and sz != m.size)) and verbose:
            print("nside={0:d}, sz={1:d}, m.size={2:d}".format(nside, sz, m.size))
            raise ValueError("Wrong nside parameter.")
        if not nest is None:  # no conversion with None
            if nest and ordering == "RING":
                idx = pixelfunc.nest2ring(nside, np.arange(m.size, dtype=np.int32))
                m = m[idx]
                if verbose:
                    print("Ordering converted to NEST")
            elif (not nest) and ordering == "NESTED":
                idx = pixelfunc.ring2nest(nside, np.arange(m.size, dtype=np.int32))
                m = m[idx]
                if verbose:
                    print("Ordering converted to RING")
        try:
            m[pixelfunc.mask_bad(m)] = UNSEEN
        except OverflowError:
            pass
        ret.append(m)

    if h:
        header = []
        for (key, value) in fits_hdu.header.items():
            header.append((key, value))

    if len(ret) == 1:
        if h:
            return ret[0], header
        else:
            return ret[0]
    else:
        if all(dt == dtype[0] for dt in dtype):
            ret = np.array(ret)
        if h:
            return ret, header
        else:
            return ret


def write_alm(
    filename, alms, out_dtype=None, lmax=-1, mmax=-1, mmax_in=-1, overwrite=False
):
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
    alms : array, complex or list of arrays
      A complex ndarray holding the alms, index = m*(2*lmax+1-m)/2+l, see Alm.getidx
    lmax : int, optional
      The maximum l in the output file
    mmax : int, optional
      The maximum m in the output file
    out_dtype : data type, optional
      data type in the output file (must be a numpy dtype). Default: *alms*.real.dtype
    mmax_in : int, optional
      maximum m in the input array
    """

    if not cb.is_seq_of_seq(alms):
        alms = [alms]

    l2max = Alm.getlmax(len(alms[0]), mmax=mmax_in)
    if lmax != -1 and lmax > l2max:
        raise ValueError("Too big lmax in parameter")
    elif lmax == -1:
        lmax = l2max

    if mmax_in == -1:
        mmax_in = l2max

    if mmax == -1:
        mmax = lmax
    if mmax > mmax_in:
        mmax = mmax_in

    if out_dtype is None:
        out_dtype = alms[0].real.dtype

    l, m = Alm.getlm(lmax)
    idx = np.where((l <= lmax) * (m <= mmax))
    l = l[idx]
    m = m[idx]

    idx_in_original = Alm.getidx(l2max, l=l, m=m)

    index = l ** 2 + l + m + 1

    hdulist = pf.HDUList()
    for alm in alms:
        out_data = np.empty(
            len(index), dtype=[("index", "i"), ("real", out_dtype), ("imag", out_dtype)]
        )
        out_data["index"] = index
        out_data["real"] = alm.real[idx_in_original]
        out_data["imag"] = alm.imag[idx_in_original]

        cindex = pf.Column(
            name="index",
            format=getformat(np.int32),
            unit="l*l+l+m+1",
            array=out_data["index"],
        )
        creal = pf.Column(
            name="real",
            format=getformat(out_dtype),
            unit="unknown",
            array=out_data["real"],
        )
        cimag = pf.Column(
            name="imag",
            format=getformat(out_dtype),
            unit="unknown",
            array=out_data["imag"],
        )

        tbhdu = pf.BinTableHDU.from_columns([cindex, creal, cimag])
        hdulist.append(tbhdu)
    hdulist.writeto(filename, overwrite=overwrite)


def read_alm(filename, hdu=1, return_mmax=False):
    """Read alm from a fits file.

    In the fits file, the alm are written
    with explicit index scheme, index = l**2+l+m+1, while healpix cxx
    uses index = m*(2*lmax+1-m)/2+l. The conversion is done in this
    function.

    Parameters
    ----------
    filename : str or HDUList or HDU
      The name of the fits file to read
    hdu : int, or tuple of int, optional
      The header to read. Start at 0. Default: hdu=1
    return_mmax : bool, optional
      If true, both the alms and mmax is returned in a tuple. Default: return_mmax=False

    Returns
    -------
    alms[, mmax] : complex array or tuple of a complex array and an int
      The alms read from the file and optionally mmax read from the file
    """
    alms = []
    lmaxtot = None
    mmaxtot = None
    for unit in np.atleast_1d(hdu):
        idx, almr, almi = mrdfits(filename, hdu=unit)
        l = np.floor(np.sqrt(idx - 1)).astype(np.long)
        m = idx - l ** 2 - l - 1
        if (m < 0).any():
            raise ValueError("Negative m value encountered !")
        lmax = l.max()
        mmax = m.max()
        if lmaxtot is None:
            lmaxtot = lmax
            mmaxtot = mmax
        else:
            if lmaxtot != lmax or mmaxtot != mmax:
                raise RuntimeError(
                    "read_alm: harmonic expansion order in {} HDUs {} does not "
                    "match".format(filename, unit, hdu)
                )
        alm = almr * (0 + 0j)
        i = Alm.getidx(lmax, l, m)
        alm.real[i] = almr
        alm.imag[i] = almi
        alms.append(alm)
    if len(alms) == 1:
        alm = alms[0]
    else:
        alm = np.array(alms)
    if return_mmax:
        return alm, mmax
    else:
        return alm


## Generic functions to read and write column of data in fits file


def _get_hdu(input_data, hdu=None, memmap=None):
    """
    Return an HDU from a FITS file

    Parameters
    ----------
    input_data : str or HDUList or HDU instance
        The input FITS file, either as a filename, HDU list, or HDU instance.

    Returns
    -------
    fits_hdu : HDU
        The extracted HDU
    """

    if isinstance(input_data, six.string_types):
        hdulist = pf.open(input_data, memmap=memmap)
        return _get_hdu(hdulist, hdu=hdu)

    if isinstance(input_data, pf.HDUList):
        if isinstance(hdu, int) and hdu >= len(input_data):
            raise ValueError("Available hdu in [0-%d]" % len(input_data))
        else:
            fits_hdu = input_data[hdu]
    elif isinstance(
        input_data,
        (pf.PrimaryHDU, pf.ImageHDU, pf.BinTableHDU, pf.TableHDU, pf.GroupsHDU),
    ):
        fits_hdu = input_data
    else:
        raise TypeError(
            "First argument should be a input_data, HDUList instance, or HDU instance"
        )

    return fits_hdu


def mrdfits(filename, hdu=1):
    """
    Read a table in a fits file.

    Parameters
    ----------
    filename : str or HDUList or HDU
      The name of the fits file to read, or an HDUList or HDU instance.
    hdu : int, optional
      The header to read. Start at 0. Default: hdu=1

    Returns
    -------
    cols : a list of arrays
      A list of column data in the given header
    """
    fits_hdu = _get_hdu(filename, hdu=hdu)
    val = []
    for i in range(len(fits_hdu.columns)):
        val.append(fits_hdu.data.field(i))
    return val


def mwrfits(filename, data, hdu=1, colnames=None, keys=None):
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
        colnames = [""] * len(data)
    cols = []
    for line in six.moves.xrange(len(data)):
        cols.append(
            pf.Column(
                name=colnames[line], format=getformat(data[line]), array=data[line]
            )
        )
    tbhdu = pf.BinTableHDU.from_columns(cols)
    if type(keys) is dict:
        for k, v in keys.items():
            tbhdu.header[k] = v
    # write the file
    tbhdu.writeto(filename)


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
        np.dtype(np.bool): "L",
        np.dtype(np.uint8): "B",
        np.dtype(np.int16): "I",
        np.dtype(np.int32): "J",
        np.dtype(np.int64): "K",
        np.dtype(np.float32): "E",
        np.dtype(np.float64): "D",
        np.dtype(np.complex64): "C",
        np.dtype(np.complex128): "M",
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
            return "A"
    except:
        pass
    try:
        if type(t) is str:
            return "A%d" % (len(t))
    except:
        pass
    try:
        if type(t[0]) is str:
            l = max(len(s) for s in t)
            return "A%d" % (l)
    except:
        pass

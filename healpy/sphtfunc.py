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
import logging

log = logging.getLogger("healpy")
import numpy as np

import astropy.io.fits as pf
from .utils.deprecation import deprecated_renamed_argument
from astropy.utils import data

DATAURL = "https://healpy.github.io/healpy-data/"
DATAURL_MIRROR = "https://github.com/healpy/healpy-data/releases/download/"

from . import _healpy_sph_transform_lib as sphtlib
from . import _sphtools as _sphtools
from . import cookbook as cb

import os.path
from . import pixelfunc

from .pixelfunc import maptype, UNSEEN, ma_to_array, accept_ma

from ._sphtools import alm2map_spin_healpy as alm2map_spin
from ._sphtools import map2alm_spin_healpy as map2alm_spin


class FutureChangeWarning(UserWarning):
    pass


DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MAX_NSIDE = (
    8192  # The maximum nside up to which most operations (e.g. map2alm) will work
)

# Spherical harmonics transformation
def anafast(
    map1,
    map2=None,
    nspec=None,
    lmax=None,
    mmax=None,
    iter=3,
    alm=False,
    pol=True,
    use_weights=False,
    datapath=None,
    gal_cut=0,
    use_pixel_weights=False,
):
    """Computes the power spectrum of a Healpix map, or the cross-spectrum
    between two maps if *map2* is given.
    No removal of monopole or dipole is performed. The input maps must be
    in ring-ordering.
    Spherical harmonics transforms in HEALPix are always on the full sky,
    if the map is masked, those pixels are set to 0. It is recommended to
    remove monopole from the map before running `anafast` to reduce
    boundary effects.

    For recommendations about how to set `lmax`, `iter`, and weights, see the
    `Anafast documentation <https://healpix.sourceforge.io/html/fac_anafast.htm>`_

    Parameters
    ----------
    map1 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps. Must be in ring ordering.
    map2 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps. Must be in ring ordering.
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
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TEB cl's and
      correlations (input must be 1 or 3 maps).
      If False, maps are assumed to be described by spin 0 spherical harmonics.
      (input can be any number of maps)
      If there is only one input map, it has no effect. Default: True.
    datapath : None or str, optional
      If given, the directory where to find the weights data.
      See the docstring of `map2alm` for details on how to set it up
    gal_cut : float [degrees]
      pixels at latitude in [-gal_cut;+gal_cut] are not taken into account
    use_pixel_weights: bool, optional
      If True, use pixel by pixel weighting, healpy will automatically download the weights, if needed
      See the map2alm docs for details about weighting

    Returns
    -------
    res : array or sequence of arrays
      If *alm* is False, returns cl or a list of cl's (TT, EE, BB, TE, EB, TB for
      polarized input map)
      Otherwise, returns a tuple (cl, alm), where cl is as above and
      alm is the spherical harmonic transform or a list of almT, almE, almB
      for polarized input
    """

    map1 = ma_to_array(map1)
    alms1 = map2alm(
        map1,
        lmax=lmax,
        mmax=mmax,
        pol=pol,
        iter=iter,
        use_weights=use_weights,
        datapath=datapath,
        gal_cut=gal_cut,
        use_pixel_weights=use_pixel_weights,
    )
    if map2 is not None:
        map2 = ma_to_array(map2)
        alms2 = map2alm(
            map2,
            lmax=lmax,
            mmax=mmax,
            pol=pol,
            iter=iter,
            use_weights=use_weights,
            datapath=datapath,
            gal_cut=gal_cut,
            use_pixel_weights=use_pixel_weights,
        )
    else:
        alms2 = None

    cls = alm2cl(alms1, alms2=alms2, lmax=lmax, mmax=mmax, lmax_out=lmax, nspec=nspec)

    if alm:
        if map2 is not None:
            return (cls, alms1, alms2)
        else:
            return (cls, alms1)
    else:
        return cls


@deprecated_renamed_argument("verbose", None, "1.15.0")
def map2alm(
    maps,
    lmax=None,
    mmax=None,
    iter=3,
    pol=True,
    use_weights=False,
    datapath=None,
    gal_cut=0,
    use_pixel_weights=False,
    verbose=True,
):
    """Computes the alm of a Healpix map. The input maps must all be
    in ring ordering.

    For recommendations about how to set `lmax`, `iter`, and weights, see the
    `Anafast documentation <https://healpix.sourceforge.io/html/fac_anafast.htm>`_

    Pixel values are weighted before applying the transform:

    * when you don't specify any weights, the uniform weight value 4*pi/n_pix is used
    * with ring weights enabled (use_weights=True), pixels in every ring
      are weighted with a uniform value similar to the one above, ring weights are
      included in healpy
    * with pixel weights (use_pixel_weights=True), every pixel gets an individual weight

    Pixel weights provide the most accurate transform, so you should always use them if
    possible. However they are not included in healpy and will be automatically downloaded
    and cached in ~/.astropy the first time you compute a trasform at a specific nside.

    If datapath is specified, healpy will first check that local folder before downloading
    the weights.
    The easiest way to setup the folder is to clone the healpy-data repository:

    git clone --depth 1 https://github.com/healpy/healpy-data
    cd healpy-data
    bash download_weights_8192.sh

    and set datapath to the root of the repository.

    Parameters
    ----------
    maps : array-like, shape (Npix,) or (n, Npix)
      The input map or a list of n input maps. Must be in ring ordering.
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    iter : int, scalar, optional
      Number of iteration (default: 3)
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TEB alm's.
      (input must be 1 or 3 maps)
      If False, apply spin 0 harmonic transform to each map.
      (input can be any number of maps)
      If there is only one input map, it has no effect. Default: True.
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    datapath : None or str, optional
      If given, the directory where to find the pixel weights.
      See in the docstring above details on how to set it up.
    gal_cut : float [degrees]
      pixels at latitude in [-gal_cut;+gal_cut] are not taken into account
    use_pixel_weights: bool, optional
      If True, use pixel by pixel weighting, healpy will automatically download the weights, if needed
    verbose : bool, optional
      Deprecated, has not effect.

    Returns
    -------
    alms : array or tuple of array
      alm or a tuple of 3 alm (almT, almE, almB) if polarized input.

    Notes
    -----
    The pixels which have the special `UNSEEN` value are replaced by zeros
    before spherical harmonic transform. They are converted back to `UNSEEN`
    value, so that the input maps are not modified. Each map have its own,
    independent mask.
    """
    maps = ma_to_array(maps)
    info = maptype(maps)
    nside = pixelfunc.get_nside(maps)
    check_max_nside(nside)

    pixel_weights_filename = None
    if use_pixel_weights:
        if use_weights:
            raise RuntimeError("Either use pixel or ring weights")
        filename = "full_weights/healpix_full_weights_nside_%04d.fits" % nside
        if datapath is not None:
            pixel_weights_filename = os.path.join(datapath, filename)
            if os.path.exists(pixel_weights_filename):
                log.info("Accessing pixel weights from %s", pixel_weights_filename)
            else:
                raise RuntimeError(
                    "You specified datapath but pixel weights file"
                    "is missing at {}".format(pixel_weights_filename)
                )
        if pixel_weights_filename is None:
            with data.conf.set_temp("dataurl", DATAURL), data.conf.set_temp(
                "dataurl_mirror", DATAURL_MIRROR
            ), data.conf.set_temp("remote_timeout", 30):
                pixel_weights_filename = data.get_pkg_data_filename(
                    filename, package="healpy"
                )

    if pol or info in (0, 1):
        alms = _sphtools.map2alm(
            maps,
            niter=iter,
            datapath=datapath,
            use_weights=use_weights,
            lmax=lmax,
            mmax=mmax,
            gal_cut=gal_cut,
            pixel_weights_filename=pixel_weights_filename,
        )
    else:
        # info >= 2 and pol is False : spin 0 spht for each map
        alms = [
            _sphtools.map2alm(
                mm,
                niter=iter,
                datapath=datapath,
                use_weights=use_weights,
                lmax=lmax,
                mmax=mmax,
                gal_cut=gal_cut,
                pixel_weights_filename=pixel_weights_filename,
            )
            for mm in maps
        ]
    return np.array(alms)


def map2alm_lsq(maps, lmax, mmax, pol=True, tol=1e-10, maxiter=20):
    """Runs an iterative map analysis up to (lmax, mmax) and returns the result
    including its quality.

    Healpix map analysis is often interpreted as "compute the `alm` for which
    `alm2map(alm) == map`". Unfortunately this inversion problem is not solvable
    in many cases, since `alm` typically has fewer elements than `map`, which
    makes the equation system overdetermined, so that a solution only exists
    for a vanishingly small subset of all possible maps. (Even if there were
    more `alm` elements than map pixels it can happen that the problem is
    unsolvable, but this situation should not be relevant for practical use.)

    This function aims to compute instead the `alm` for which the L2 norm of
    `alm2map(alm) - map` is minimal, i.e. it tries to find the best possible
    least-squares approximation. Since this is an iterative procedure, the user
    needs to specify a tolerance at which the iteration is stopped and a maximum
    iteration count to avoid excessive run times in extreme cases.

    Compared to `map2alm` this algorithm has the following advantages and
    disadvantages:
     - the exact number of iterations need not be specified beforehand, it will
       stop automatically once the desired accuracy is reached.
     - during calculation it requires more memory than `map2alm`.

    Parameters
    ----------
    maps : array-like, shape (Npix,) or (n, Npix)
      The input map or a list of n input maps. Must be in ring ordering.
    lmax, mmax : int
        The desired lmax and mmax parameters for the analysis
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TEB alm's.
      (input must be 1 or 3 maps)
      If False, apply spin 0 harmonic transform to each map.
      (input can be any number of maps)
      If there is only one input map, it has no effect. Default: True.
    tol : float
        The desired accuracy for the result. Once this is reached, the iteration
        stops.
    maxiter : int
        maximum iteration count after which the minimization is stopped

    Returns
    -------
    alm : numpy.ndarray(complex)
        The reconstructed a_lm coefficients
    rel_res : float
        The norm of the residual map (i.e. `map-alm2map(alm)`), divided by the
        norm of `maps`. This is a measure for the fraction of the map signal that
        could not be modeled by the a_lm
    n_iter : int
        the number of iterations required
    """
    from scipy.sparse.linalg import LinearOperator, lsqr, lsmr
    from scipy.linalg import norm

    maps = ma_to_array(maps)
    info = maptype(maps)
    ncomp = 1 if info == 0 else info
    maps = np.array(maps)
    npix = maps.size / ncomp
    nside = pixelfunc.get_nside(maps)
    check_max_nside(nside)

    # helper functions to convert between real- and complex-valued a_lm
    def alm2realalm(alm):
        alm = alm.reshape((ncomp, -1))
        res = np.zeros((ncomp, alm.shape[1] * 2 - lmax - 1))
        for i in range(ncomp):
            res[i, 0 : lmax + 1] = alm[i, 0 : lmax + 1].real
            res[i, lmax + 1 :] = alm[i, lmax + 1 :].view(np.float64) * np.sqrt(2.0)
        return res.reshape((-1,))

    def realalm2alm(alm):
        alm = np.array(alm).reshape((ncomp, -1))
        res = np.zeros((ncomp, (alm.shape[1] + lmax + 1) // 2), dtype=np.complex128)
        for i in range(ncomp):
            res[i, 0 : lmax + 1] = alm[i, 0 : lmax + 1]
            res[i, lmax + 1 :] = alm[i, lmax + 1 :].view(np.complex128) * (
                np.sqrt(2.0) / 2
            )
        if ncomp == 1:
            res = res[0]
        return res

    def a2m2(x):
        talm = realalm2alm(x)
        res = alm2map(talm, lmax=lmax, nside=nside, pol=pol)
        return res.reshape((-1,))

    def m2a2(x):
        talm = map2alm(x.reshape((ncomp, -1)), lmax=lmax, iter=0, pol=pol)
        talm *= (12 * nside ** 2) / (4 * np.pi)
        res = alm2realalm(talm)
        return res

    # initial guess
    alm0 = m2a2(maps) / npix * (4 * np.pi)
    op = LinearOperator(
        matvec=a2m2, rmatvec=m2a2, shape=(len(maps) * maps[0].size, alm0.size)
    )
    res = lsqr(
        A=op,
        b=maps.reshape((-1,)),
        x0=alm0,
        atol=tol,
        btol=tol,
        iter_lim=maxiter,
        show=False,
    )
    #    res = lsmr(A=op, b=maps, x0=alm0, atol=tol, btol=tol, maxiter=maxiter, show=False)
    return realalm2alm(res[0]), res[3] / norm(maps.reshape((-1,))), res[2]


@deprecated_renamed_argument("verbose", None, "1.15.0")
def alm2map(
    alms,
    nside,
    lmax=None,
    mmax=None,
    pixwin=False,
    fwhm=0.0,
    sigma=None,
    pol=True,
    inplace=False,
    verbose=True,
):
    """Computes a Healpix map given the alm.

    The alm are given as a complex array. You can specify lmax
    and mmax, or they will be computed from array size (assuming
    lmax==mmax).

    Parameters
    ----------
    alms : complex, array or sequence of arrays
      A complex array or a sequence of complex arrays.
      Each array must have a size of the form: mmax * (2 * lmax + 1 - mmax) / 2 + lmax + 1
    nside : int, scalar
      The nside of the output map.
    lmax : None or int, scalar, optional
      Explicitly define lmax (needed if mmax!=lmax)
    mmax : None or int, scalar, optional
      Explicitly define mmax (needed if mmax!=lmax)
    pixwin : bool, optional
      Smooth the alm using the pixel window functions. Default: False.
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    pol : bool, optional
      If True, assumes input alms are TEB. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, apply spin 0 harmonic transform to each alm.
      (input can be any number of alms)
      If there is only one input alm, it has no effect. Default: True.
    inplace : bool, optional
      If True, input alms may be modified by pixel window function and beam
      smoothing (if alm(s) are complex128 contiguous arrays).
      Otherwise, input alms are not modified. A copy is made if needed to
      apply beam smoothing or pixel window.

    Returns
    -------
    maps : array or list of arrays
      A Healpix map in RING scheme at nside or a list of T,Q,U maps (if
      polarized input)

    Notes
    -----
    Running map2alm then alm2map will not return exactly the same map if the discretized field you construct on the sphere is not band-limited (for example, if you have a map containing pixel-based noise rather than beam-smoothed noise). If you need a band-limited map, you have to start with random numbers in lm space and transform these via alm2map. With such an input, the accuracy of map2alm->alm2map should be quite good, depending on your choices of lmax, mmax and nside (for some typical values, see e.g., section 5.1 of https://arxiv.org/pdf/1010.2084).
    """
    if not cb.is_seq(alms):
        raise TypeError("alms must be a sequence")

    check_max_nside(nside)

    alms = smoothalm(alms, fwhm=fwhm, sigma=sigma, pol=pol, inplace=inplace)

    if not cb.is_seq_of_seq(alms):
        alms = [alms]
        lonely = True
    else:
        lonely = False

    if pixwin:
        pw = globals()["pixwin"](nside, True)
        alms_new = []
        for ialm, alm in enumerate(alms):
            pixelwindow = pw[1] if ialm >= 1 and pol else pw[0]
            alms_new.append(almxfl(alm, pixelwindow, inplace=inplace))
    else:
        alms_new = alms

    if lmax is None:
        lmax = -1
    if mmax is None:
        mmax = -1
    if pol:
        output = sphtlib._alm2map(
            alms_new[0] if lonely else tuple(alms_new), nside, lmax=lmax, mmax=mmax
        )
        if lonely:
            output = [output]
    else:
        output = [
            sphtlib._alm2map(alm, nside, lmax=lmax, mmax=mmax) for alm in alms_new
        ]
    if lonely:
        return output[0]
    else:
        return np.array(output)


@deprecated_renamed_argument("verbose", None, "1.15.0")
def synalm(cls, lmax=None, mmax=None, new=False, verbose=True):
    """Generate a set of alm given cl.
    The cl are given as a float array. Corresponding alm are generated.
    If lmax is None, it is assumed lmax=cl.size-1
    If mmax is None, it is assumed mmax=lmax.

    Parameters
    ----------
    cls : float, array or tuple of arrays
      Either one cl (1D array) or a tuple of either 4 cl
      or of n*(n+1)/2 cl.
      Some of the cl may be None, implying no
      cross-correlation. See *new* parameter.
    lmax : int, scalar, optional
      The lmax (if None or <0, the largest size-1 of cls)
    mmax : int, scalar, optional
      The mmax (if None or <0, =lmax)
    new : bool, optional
      If True, use the new ordering of cl's, ie by diagonal
      (e.g. TT, EE, BB, TE, EB, TB or TT, EE, BB, TE if 4 cl as input).
      If False, use the old ordering, ie by row
      (e.g. TT, TE, TB, EE, EB, BB or TT, TE, EE, BB if 4 cl as input).

    Returns
    -------
    alms : array or list of arrays
      the generated alm if one spectrum is given, or a list of n alms
      (with n(n+1)/2 the number of input cl, or n=3 if there are 4 input cl).

    Notes
    -----
    We don't plan to change the default order anymore, that would break old
    code in a way difficult to debug.
    """
    if not cb.is_seq(cls):
        raise TypeError("cls must be an array or a sequence of arrays")

    if not cb.is_seq_of_seq(cls, True):
        # Only one spectrum
        if lmax is None or lmax < 0:
            lmax = len(cls) - 1
        if mmax is None or mmax < 0:
            mmax = lmax
        cls_list = [np.asarray(cls, dtype=np.float64)]
        szalm = Alm.getsize(lmax, mmax)
        alm = np.zeros(szalm, "D")
        alm.real = np.random.standard_normal(szalm)
        alm.imag = np.random.standard_normal(szalm)
        alms_list = [alm]
        sphtlib._synalm(cls_list, alms_list, lmax, mmax)
        return alm

    # From here, we interpret cls as a list of spectra
    cls_list = list(cls)
    maxsize = max([len(c) for c in cls if c is not None])

    if lmax is None or lmax < 0:
        lmax = maxsize - 1
    if mmax is None or mmax < 0:
        mmax = lmax

    Nspec = sphtlib._getn(len(cls_list))

    if Nspec <= 0:
        if len(cls_list) == 4:
            if new:  ## new input order: TT EE BB TE -> TT EE BB TE 0 0
                cls_list = [cls[0], cls[1], cls[2], cls[3], None, None]
            else:  ## old input order: TT TE EE BB -> TT TE 0 EE 0 BB
                cls_list = [cls[0], cls[1], None, cls[2], None, cls[3]]
            Nspec = 3
        else:
            raise TypeError(
                "The sequence of arrays must have either 4 elements "
                "or n(n+1)/2 elements (some may be None)"
            )

    szalm = Alm.getsize(lmax, mmax)
    alms_list = []
    for i in range(Nspec):
        alm = np.zeros(szalm, "D")
        alm.real = np.random.standard_normal(szalm)
        alm.imag = np.random.standard_normal(szalm)
        alms_list.append(alm)
    if new:  # new input order: input given by diagonal, should be given by row
        cls_list = new_to_old_spectra_order(cls_list)
    # ensure cls are float64
    cls_list = [
        (np.asarray(cl, dtype=np.float64) if cl is not None else None)
        for cl in cls_list
    ]
    sphtlib._synalm(cls_list, alms_list, lmax, mmax)
    return np.array(alms_list)


@deprecated_renamed_argument("verbose", None, "1.15.0")
def synfast(
    cls,
    nside,
    lmax=None,
    mmax=None,
    alm=False,
    pol=True,
    pixwin=False,
    fwhm=0.0,
    sigma=None,
    new=False,
    verbose=True,
):
    """Create a map(s) from cl(s).

    You can choose a random seed using `numpy.random.seed(SEEDVALUE)`
    before calling `synfast`.

    Parameters
    ----------
    cls : array or tuple of array
      A cl or a list of cl (either 4 or 6, see :func:`synalm`)
    nside : int, scalar
      The nside of the output map(s)
    lmax : int, scalar, optional
      Maximum l for alm. Default: min of 3*nside-1 or length of the cls - 1
    mmax : int, scalar, optional
      Maximum m for alm.
    alm : bool, scalar, optional
      If True, return also alm(s). Default: False.
    pol : bool, optional
      If True, assumes input cls are TEB and correlation. Output will be TQU maps.
      (input must be 1, 4 or 6 cl's)
      If False, fields are assumed to be described by spin 0 spherical harmonics.
      (input can be any number of cl's)
      If there is only one input cl, it has no effect. Default: True.
    pixwin : bool, scalar, optional
      If True, convolve the alm by the pixel window function. Default: False.
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    new : bool, optional
      If True, use the new ordering of cl's, ie by diagonal
      (e.g. TT, EE, BB, TE, EB, TB or TT, EE, BB, TE if 4 cl as input).
      If False, use the old ordering, ie by row
      (e.g. TT, TE, TB, EE, EB, BB or TT, TE, EE, BB if 4 cl as input).

    Returns
    -------
    maps : array or tuple of arrays
      The output map (possibly list of maps if polarized input).
      or, if alm is True, a tuple of (map,alm)
      (alm possibly a list of alm if polarized input)

    Notes
    -----
    We don't plan to change the default order anymore, that would break old
    code in a way difficult to debug.
    """
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    check_max_nside(nside)

    cls_lmax = cb.len_array_or_arrays(cls) - 1
    if lmax is None or lmax < 0:
        lmax = min(cls_lmax, 3 * nside - 1)
    alms = synalm(cls, lmax=lmax, mmax=mmax, new=new)
    maps = alm2map(
        alms,
        nside,
        lmax=lmax,
        mmax=mmax,
        pixwin=pixwin,
        pol=pol,
        fwhm=fwhm,
        sigma=sigma,
        inplace=True,
    )
    if alm:
        return np.array(maps), np.array(alms)
    else:
        return np.array(maps)


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
    def getlm(lmax, i=None):
        """Get the l and m from index and lmax.

        Parameters
        ----------
        lmax : int
          The maximum l defining the alm layout
        i : int or None
          The index for which to compute the l and m.
          If None, the function return l and m for i=0..Alm.getsize(lmax)
        """
        szalm = Alm.getsize(lmax, lmax)
        if i is None:
            i = np.arange(szalm)
        assert (
            np.max(i) < szalm
        ), "Invalid index, it should less than the max alm array length of {}".format(
            szalm
        )

        with np.errstate(all="raise"):
            m = (
                np.ceil(
                    ((2 * lmax + 1) - np.sqrt((2 * lmax + 1) ** 2 - 8 * (i - lmax))) / 2
                )
            ).astype(int)
            l = i - m * (2 * lmax + 1 - m) // 2
        return (l, m)

    @staticmethod
    def getidx(lmax, l, m):
        r"""Returns index corresponding to (l,m) in an array describing alm up to lmax.

        In HEALPix C++ and healpy, :math:`a_{lm}` coefficients are stored ordered by
        :math:`m`. I.e. if :math:`\ell_{max}` is 16, the first 16 elements are
        :math:`m=0, \ell=0-16`, then the following 15 elements are :math:`m=1, \ell=1-16`,
        then :math:`m=2, \ell=2-16` and so on until the last element, the 153th, is
        :math:`m=16, \ell=16`.

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
        return m * (2 * lmax + 1 - m) // 2 + l

    @staticmethod
    def getsize(lmax, mmax=None):
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
        return mmax * (2 * lmax + 1 - mmax) // 2 + lmax + 1

    @staticmethod
    def getlmax(s, mmax=None):
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


def alm2cl(alms1, alms2=None, lmax=None, mmax=None, lmax_out=None, nspec=None):
    """Computes (cross-)spectra from alm(s). If alm2 is given, cross-spectra between
    alm and alm2 are computed. If alm (and alm2 if provided) contains n alm,
    then n(n+1)/2 auto and cross-spectra are returned.

    Parameters
    ----------
    alm : complex, array or sequence of arrays
      The alm from which to compute the power spectrum. If n>=2 arrays are given,
      computes both auto- and cross-spectra.
    alms2 : complex, array or sequence of 3 arrays, optional
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
    cls = _sphtools.alm2cl(alms1, alms2=alms2, lmax=lmax, mmax=mmax, lmax_out=lmax_out)
    if nspec is None:
        return np.array(cls)
    else:
        return np.array(cls[:nspec])


def almxfl(alm, fl, mmax=None, inplace=False):
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
      The modified alm, either a new array or a reference to input alm,
      if inplace is True.
    """
    # FIXME: Should handle multidimensional input
    almout = _sphtools.almxfl(alm, fl, mmax=mmax, inplace=inplace)
    return almout


@deprecated_renamed_argument("verbose", None, "1.15.0")
def smoothalm(
    alms,
    fwhm=0.0,
    sigma=None,
    beam_window=None,
    pol=True,
    mmax=None,
    verbose=True,
    inplace=True,
):
    """Smooth alm with a Gaussian symmetric beam function.

    Parameters
    ----------
    alms : array or sequence of 3 arrays
      Either an array representing one alm, or a sequence of arrays.
      See *pol* parameter.
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
      [in radians]
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.
      [in radians]
    beam_window: array, optional
      Custom beam window function. Override fwhm and sigma.
    pol : bool, optional
      If True, assumes input alms are TEB. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, apply spin 0 harmonic transform to each alm.
      (input can be any number of alms)
      If there is only one input alm, it has no effect. Default: True.
    mmax : None or int, optional
      The maximum m for alm. Default: mmax=lmax
    inplace : bool, optional
      If True, the alm's are modified inplace if they are contiguous arrays
      of type complex128. Otherwise, a copy of alm is made. Default: True.
    verbose : bool, optional
      Deprecated, has not effect.

    Returns
    -------
    alms : array or sequence of 3 arrays
      The smoothed alm. If alm[i] is a contiguous array of type complex128,
      and *inplace* is True the smoothing is applied inplace.
      Otherwise, a copy is made.
    """
    if (sigma is None) & (beam_window is None):
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    if beam_window is None:
        log.info("Sigma is %f arcmin (%f rad) ", sigma * 60 * 180 / np.pi, sigma)
        log.info(
            "-> fwhm is %f arcmin",
            sigma * 60 * 180 / np.pi * (2.0 * np.sqrt(2.0 * np.log(2.0))),
        )
    else:
        if not isinstance(beam_window, (np.ndarray, list)):
            raise ValueError("beam_window must be an array or list")

    # Check alms
    if not cb.is_seq(alms):
        raise ValueError("alm must be a sequence")

    if (sigma == 0) & (beam_window is None):
        # nothing to be done
        return alms

    lonely = False
    if not cb.is_seq_of_seq(alms):
        alms = [alms]
        lonely = True

    # we have 3 alms -> apply smoothing to each map.
    # polarization has different B_l from temperature
    # exp{-[ell(ell+1) - s**2] * sigma**2/2}
    # with s the spin of spherical harmonics
    # s = 2 for pol, s=0 for temperature
    retalm = []
    for ialm, alm in enumerate(alms):
        lmax = Alm.getlmax(len(alm), mmax)
        if lmax < 0:
            raise TypeError(
                "Wrong alm size for the given "
                "mmax (len(alms[%d]) = %d)." % (ialm, len(alm))
            )
        ell = np.arange(lmax + 1.0)
        s = 2 if ialm >= 1 and pol else 0
        if beam_window is None:
            fact = np.exp(-0.5 * (ell * (ell + 1) - s ** 2) * sigma ** 2)
        else:
            fact = np.copy(beam_window)
        res = almxfl(alm, fact, mmax=mmax, inplace=inplace)
        retalm.append(res)
    # Test what to return (inplace/not inplace...)
    # Case 1: 1d input, return 1d output
    if lonely:
        return retalm[0]
    # case 2: 2d input, check if in-place smoothing for all alm's
    for i in range(len(alms)):
        samearray = alms[i] is retalm[i]
        if not samearray:
            # Case 2a:
            # at least one of the alm could not be smoothed in place:
            # return the list of alm
            return np.array(retalm)
    # Case 2b:
    # all smoothing have been performed in place:
    # return the input alms.  If the input was a tuple, so will the output be.
    return alms


@deprecated_renamed_argument("verbose", None, "1.15.0")
@accept_ma
def smoothing(
    map_in,
    fwhm=0.0,
    sigma=None,
    beam_window=None,
    pol=True,
    iter=3,
    lmax=None,
    mmax=None,
    use_weights=False,
    use_pixel_weights=False,
    datapath=None,
    verbose=True,
    nest=False,
):
    """Smooth a map with a Gaussian symmetric beam.

    No removal of monopole or dipole is performed.

    Parameters
    ----------
    map_in : array or sequence of 3 arrays
      Either an array representing one map, or a sequence of
      3 arrays representing 3 maps, accepts masked arrays
    fwhm : float, optional
      The full width half max parameter of the Gaussian [in
      radians]. Default:0.0
    sigma : float, optional
      The sigma of the Gaussian [in radians]. Override fwhm.
    beam_window: array, optional
      Custom beam window function. Override fwhm and sigma.
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, each map is assumed to be a spin 0 map and is
      treated independently (input can be any number of alms).
      If there is only one input map, it has no effect. Default: True.
    iter : int, scalar, optional
      Number of iteration (default: 3)
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    use_pixel_weights: bool, optional
      If True, use pixel by pixel weighting, healpy will automatically download the weights, if needed
      See the map2alm docs for details about weighting
    datapath : None or str, optional
      If given, the directory where to find the pixel weights data.
      See the docstring of `map2alm` for details on how to set it up
    verbose : bool, optional
      Deprecated, has not effect.
    nest : bool, optional
      If True, the input map ordering is assumed to be NESTED. Default: False (RING)
      This function will temporary reorder the NESTED map into RING to perform the
      smoothing and order the output back to NESTED. If the map is in RING ordering
      no internal reordering will be performed.

    Returns
    -------
    maps : array or list of 3 arrays
      The smoothed map(s)
    """

    if not cb.is_seq(map_in):
        raise TypeError("map_in must be a sequence")

    # save the masks of inputs
    masks = pixelfunc.mask_bad(map_in)

    if cb.is_seq_of_seq(map_in):
        nside = pixelfunc.get_nside(map_in)
        n_maps = len(map_in)
    else:
        nside = pixelfunc.get_nside(map_in)
        n_maps = 0

    check_max_nside(nside)

    if nest:
        map_in = pixelfunc.reorder(map_in, inp=None, out=None, r2n=None, n2r=True)

    if pol or n_maps in (0, 1):
        # Treat the maps together (1 or 3 maps)
        alms = map2alm(
            map_in,
            lmax=lmax,
            mmax=mmax,
            iter=iter,
            pol=pol,
            use_weights=use_weights,
            use_pixel_weights=use_pixel_weights,
            datapath=datapath,
        )
        smoothalm(
            alms,
            fwhm=fwhm,
            sigma=sigma,
            beam_window=beam_window,
            pol=pol,
            mmax=mmax,
            inplace=True,
        )
        output_map = alm2map(alms, nside, lmax=lmax, mmax=mmax, pixwin=False, pol=pol)
    else:
        # Treat each map independently (any number)
        output_map = []
        for m in map_in:
            alm = map2alm(
                m,
                lmax=lmax,
                mmax=mmax,
                iter=iter,
                pol=pol,
                use_weights=use_weights,
                use_pixel_weights=use_pixel_weights,
                datapath=datapath,
            )
            smoothalm(
                alm, fwhm=fwhm, sigma=sigma, beam_window=beam_window, inplace=True
            )
            output_map.append(alm2map(alm, nside, pixwin=False))
        output_map = np.array(output_map)
    output_map[masks] = UNSEEN

    if nest:
        output_map = pixelfunc.reorder(
            output_map, inp=None, out=None, r2n=True, n2r=False
        )

    return output_map


def pixwin(nside, pol=False, lmax=None):
    """Return the pixel window function for the given nside.

    Parameters
    ----------
    nside : int
      The nside for which to return the pixel window function
    pol : bool, optional
      If True, return also the polar pixel window. Default: False
    lmax : int, optional
        Maximum l of the power spectrum (default: 3*nside-1)

    Returns
    -------
    pw or pwT,pwP : array or tuple of 2 arrays
      The temperature pixel window function, or a tuple with both
      temperature and polarisation pixel window functions.
    """

    if lmax is None:
        lmax = 3 * nside - 1

    datapath = DATAPATH
    if not pixelfunc.isnsideok(nside):
        raise ValueError("Wrong nside value (must be a power of two).")
    fname = os.path.join(datapath, "pixel_window_n%04d.fits" % nside)
    if not os.path.isfile(fname):
        raise ValueError("No pixel window for this nside " "or data files missing")
    # return hfitslib._pixwin(nside,datapath,pol)  ## BROKEN -> seg fault...
    pw = pf.getdata(fname)
    pw_temp, pw_pol = pw.field(0), pw.field(1)
    if pol:
        return pw_temp[: lmax + 1], pw_pol[: lmax + 1]
    else:
        return pw_temp[: lmax + 1]


def alm2map_der1(alm, nside, lmax=None, mmax=None):
    """Computes a Healpix map and its first derivatives given the alm.

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
      theta and phi. d_phi is already divided by sin(theta)
    """
    check_max_nside(nside)
    if lmax is None:
        lmax = -1
    if mmax is None:
        mmax = -1
    return np.array(sphtlib._alm2map_der1(alm, nside, lmax=lmax, mmax=mmax))


def new_to_old_spectra_order(cls_new_order):
    """Reorder the cls from new order (by diagonal) to old order (by row).
    For example : TT, EE, BB, TE, EB, BB => TT, TE, TB, EE, EB, BB
    """
    Nspec = sphtlib._getn(len(cls_new_order))
    if Nspec < 0:
        raise ValueError("Input must be a list of n(n+1)/2 arrays")
    cls_old_order = []
    for i in range(Nspec):
        for j in range(i, Nspec):
            p = j - i
            q = i
            idx_new = p * (2 * Nspec + 1 - p) // 2 + q
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
    cls = np.loadtxt(os.path.join(DATAPATH, "totcls.dat"), unpack=True)
    ell = cls[0]
    f = ell * (ell + 1) / 2 / np.pi
    cls[1:, 1:] /= f[1:]
    return ell, f, cls[1:]


def gauss_beam(fwhm, lmax=512, pol=False):
    """Gaussian beam window function

    Computes the spherical transform of an axisimmetric gaussian beam

    For a sky of underlying power spectrum C(l) observed with beam of
    given FWHM, the measured power spectrum will be
    C(l)_meas = C(l) B(l)^2
    where B(l) is given by gaussbeam(Fwhm,Lmax).
    The polarization beam is also provided (when pol = True ) assuming
    a perfectly co-polarized beam
    (e.g., Challinor et al 2000, astro-ph/0008228)

    Parameters
    ----------
    fwhm : float
        full width half max in radians
    lmax : integer
        ell max
    pol : bool
        if False, output has size (lmax+1) and is temperature beam
        if True output has size (lmax+1, 4) with components:
        * temperature beam
        * grad/electric polarization beam
        * curl/magnetic polarization beam
        * temperature * grad beam

    Returns
    -------
    beam : array
        beam window function [0, lmax] if dim not specified
        otherwise (lmax+1, 4) contains polarized beam
    """

    sigma = fwhm / np.sqrt(8.0 * np.log(2.0))
    ell = np.arange(lmax + 1)
    sigma2 = sigma ** 2
    g = np.exp(-0.5 * ell * (ell + 1) * sigma2)

    if not pol:  # temperature-only beam
        return g
    else:  # polarization beam
        # polarization factors [1, 2 sigma^2, 2 sigma^2, sigma^2]
        pol_factor = np.exp([0.0, 2 * sigma2, 2 * sigma2, sigma2])
        return g[:, np.newaxis] * pol_factor


def bl2beam(bl, theta):
    """Computes a circular beam profile b(theta) in real space from
    its transfer (or window) function b(l) in spherical harmonic space.

    Parameters
    ----------
    bl : array
        Window function b(l) of the beam.
    theta : array
        Radius at which the beam profile will be computed.
        Has to be given in radians.

    Returns
    -------
    beam : array
        (Circular) beam profile b(theta).
    """

    lmax = len(bl) - 1
    nx = len(theta)
    x = np.cos(theta)
    p0 = np.zeros(nx) + 1
    p1 = x

    beam = bl[0] * p0 + bl[1] * p1 * 3

    for l in np.arange(2, lmax + 1):
        p2 = x * p1 * (2 * l - 1) / l - p0 * (l - 1) / l
        p0 = p1
        p1 = p2
        beam += bl[l] * p2 * (2 * l + 1)

    beam /= 4 * np.pi

    return beam


def beam2bl(beam, theta, lmax):
    """Computes a transfer (or window) function b(l) in spherical
    harmonic space from its circular beam profile b(theta) in real
    space.

    Parameters
    ----------
    beam : array
        Circular beam profile b(theta).
    theta : array
        Radius at which the beam profile is given. Has to be given
        in radians with same size as beam.
    lmax : integer
        Maximum multipole moment at which to compute b(l).

    Returns
    -------
    bl : array
        Beam window function b(l).
    """
    from scipy.integrate import trapz

    nx = len(theta)
    nb = len(beam)
    if nb != nx:
        raise ValueError("Beam and theta must have same size!")

    x = np.cos(theta)
    st = np.sin(theta)
    window = np.zeros(lmax + 1)

    p0 = np.ones(nx)
    p1 = np.copy(x)

    window[0] = trapz(beam * p0 * st, theta)
    window[1] = trapz(beam * p1 * st, theta)

    for l in np.arange(2, lmax + 1):
        p2 = x * p1 * (2 * l - 1) / l - p0 * (l - 1) / l
        window[l] = trapz(beam * p2 * st, theta)
        p0 = p1
        p1 = p2

    window *= 2 * np.pi

    return window


def check_max_nside(nside):
    """Checks whether the nside used in a certain operation does not exceed the
    maximum supported nside. The maximum nside is saved in MAX_NSIDE.

    Parameters
    ----------
    nside : int
        nside of the map that is being checked
    """

    if nside > MAX_NSIDE:
        raise ValueError(
            "nside {nside} of map cannot be larger than "
            "MAX_NSIDE {max_nside}".format(nside=nside, max_nside=MAX_NSIDE)
        )

    return 0


def blm_gauss(fwhm, lmax, pol=False):
    """Computes spherical harmonic coefficients of a circular Gaussian beam
    pointing towards the North Pole

    See an example of usage
    `in the documentation <https://healpy.readthedocs.io/en/latest/blm_gauss_plot.html>`_

    Parameters
    ----------
    fwhm : float, scalar
        desired FWHM of the beam, in radians
    lmax : int, scalar
        maximum l multipole moment to compute
    pol : bool, scalar
        if True, E and B coefficients will also be computed

    Returns
    -------
    blm : array with dtype numpy.complex128
          lmax will be as specified
          mmax is 0 for pol==False, else 2
    """
    fwhm = float(fwhm)
    lmax = int(lmax)
    pol = bool(pol)
    mmax = 2 if pol else 0
    ncomp = 3 if pol else 1
    nval = Alm.getsize(lmax, mmax)

    if mmax > lmax:
        raise ValueError("lmax value too small")

    blm = np.zeros((ncomp, nval), dtype=np.complex128)
    sigmasq = fwhm * fwhm / (8 * np.log(2.0))

    for l in range(0, lmax + 1):
        blm[0, Alm.getidx(lmax, l, 0)] = np.sqrt((2 * l + 1) / (4.0 * np.pi)) * np.exp(
            -0.5 * sigmasq * l * l
        )

    if pol:
        for l in range(2, lmax + 1):
            blm[1, Alm.getidx(lmax, l, 2)] = np.sqrt(
                (2 * l + 1) / (32 * np.pi)
            ) * np.exp(-0.5 * sigmasq * l * l)
        blm[2] = 1j * blm[1]

    return blm


def resize_alm(alm, lmax, mmax, lmax_out, mmax_out):
    """Returns a resized copy of the input a_lm, either truncated or
    padded with zeros.

    Parameters
    ----------
    alm : array or sequence of arrays, complex
        A complex array of alm, or a sequence of such arrays.
        Size of each array must be of the form mmax*(lmax-mmax+1)/2+lmax
    lmax, mmax: int, scalar
        maximum l and m multipole moments of the input alm
    lmax_out, mmax_out: int, scalar
        maximum l and m multipole moments of the output alm

    Returns
    -------
    array or sequence of arrays, complex
        A complex array of alm, or a sequence of such arrays.
        Size of each array will be of the form mmax_out*(lmax_out-mmax_out+1)/2+lmax_out
    """
    alm = np.array(alm)
    if alm.ndim < 1 or alm.ndim > 2:
        raise ValueError("incorrect dimensionality of the input a_lm")
    if alm.ndim == 2:
        return [resize_alm(almi, lmax, mmax, lmax_out, mmax_out) for almi in alm]
    # alm is a 1D array
    if alm.shape[0] != Alm.getsize(lmax, mmax):
        raise ValueError("inconsistent number of input a_lm")
    res = np.zeros(Alm.getsize(lmax_out, mmax_out), dtype=alm.dtype)
    lmaxmin = min(lmax, lmax_out)
    mmaxmin = min(mmax, mmax_out)
    ofs_i, ofs_o = 0, 0
    for m in range(0, mmaxmin + 1):
        nval = lmaxmin - m + 1
        res[ofs_o : ofs_o + nval] = alm[ofs_i : ofs_i + nval]
        ofs_i += lmax - m + 1
        ofs_o += lmax_out - m + 1
    return res

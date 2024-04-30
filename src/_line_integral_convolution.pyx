import numpy as np
cimport numpy as np
cimport cython
from _common cimport RING, Healpix_Map, ndarray2map
from .pixelfunc import npix2nside
from .sphtfunc import alm2map, Alm

cdef extern from "alice3.h":
    cdef void lic_main(const Healpix_Map[double] &Q, const Healpix_Map[double] &U,
                       const Healpix_Map[double] &th, Healpix_Map[double] &hit,
                       Healpix_Map[double] &tex, Healpix_Map[double] &mag,
                       int steps, int kernel_steps, double step_radian,
                       double polmin, double polmax)


def line_integral_convolution(
    Q,
    U,
    steps=100,
    kernel_steps=50,
    step_radian=0.003,
    pol_min=-1e30,
    pol_max=1e30,
    rand_seed=42,
    ell=None,
    texture=None,
    modulate=False,
):
    """Computes line integral convolution of polarized map.

    This visualizes the vector field via convolution of streamlines with a
    random background texture. It thus serves as a technique for visualizing
    polarization orientation. Polarization intensity can additionally be
    visualized by modulating it with the convolved texture by setting the
    *modulate* argument. By default, the random texture is white noise, but it
    can also be generated with power at a given angular scale if *ell* is
    specified, or a custom texture can be passed to this function with the
    *texture* argument.

    See example: https://github.com/healpy/healpy/pull/617#issue-434041253

    Parameters
    ----------
    Q : array-like, shape (Npix,)
      Input Stokes Q map. Must be in RING ordering.
    U : array-like, shape (Npix,)
      Input Stokes U map. Must be in RING ordering.
    steps : int
      Number of steps to use for each line of the convolution (default: 100).
    kernel_steps : int
      Extent of the convolution kernel (in steps, <= steps) (default: 50).
    step_radian : float
      Size (in radians) of each step in the convolution (default: 0.003).
    pol_min : float
      Minimum value for the polarization magnitude (default: -1e30).
    pol_max : float
      Maximum value for the polarization magnitude (default: 1e30).
    rand_seed : int
      Seed for the random number generator; only used when texture is not
      provided (default: 42).
    ell : int or None
      The ell value at which the generated background texture has power
      (default: ``None``). If ``ell < 0``, it is set to ``ell=2*nside``. If
      ``None``, the generated background texture will contain white noise.
    texture : array-like, shape (Npix,) or None
      Background texture map (default: ``None``). Must be in RING ordering. If
      ``None``, a background texture will be generated based on the values of
      ``rand_seed`` and ``ell``.
    modulate : bool
      Whether or not to modulate the convolved texture with the polarization
      intensity (default: ``False``).

    Returns
    -------
    lic_texture : Line integral convolved texture.

    Example
    -------
    >>> import healpy as hp
    >>> import numpy as np
    >>> import matplotlib.colors
    >>> import matplotlib.pyplot as plt
    >>> I, Q, U = hp.read_map('iqu_map.fits', (0, 1, 2))
    >>> lic = hp.line_integral_convolution(Q, U)
    >>> hp.mollview(I)
    >>> cmap_colors = plt.get_cmap('binary', 256)(np.linspace(0, 1, 256))
    >>> cmap_colors[..., 3] = 0.5  # Make colormap partially transparent
    >>> cmap = matplotlib.colors.ListedColormap(cmap_colors)
    >>> hp.mollview(lic, cmap=cmap, cbar=False, reuse_axes=True)

    """

    if Q.shape != U.shape:
        raise ValueError("Q and U maps must be the same size!")
    if texture is not None and texture.shape != Q.shape:
        raise ValueError("Texture and Q / U maps must be the same size!")
    if Q.ndim != 1:
        raise ValueError("Maps must have only a single dimension!")
    if kernel_steps > steps:
        raise ValueError("Kernel steps must be fewer than steps!")

    npix = Q.size
    nside = npix2nside(npix)
    if texture is None:
        rng = np.random.RandomState(rand_seed)
        if ell is not None:
            if ell < 0:
                ell = 2 * nside
            alms = np.zeros(Alm.getsize(ell), dtype=complex)
            for m in range(ell + 1):
                alms[Alm.getidx(ell, ell, m)] = rng.normal() + 1j * rng.normal()
            texture = alm2map(alms, nside)
        else:
            texture = rng.uniform(size=npix) - 0.5

    hit = np.empty(npix)
    tex = np.empty(npix)
    mag = np.empty(npix)

    lic_main(
        ndarray2map(Q, RING)[0],
        ndarray2map(U, RING)[0],
        ndarray2map(texture, RING)[0],
        ndarray2map(hit, RING)[0],
        ndarray2map(tex, RING)[0],
        ndarray2map(mag, RING)[0],
        steps,
        kernel_steps,
        step_radian,
        pol_min,
        pol_max,
    )

    lic_texture = mag if modulate else tex
    return lic_texture

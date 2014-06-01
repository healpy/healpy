__all__ = ['mollview', 'projplot']

import numpy as np
from .pixelfunc import ang2pix, npix2nside
from .rotator import Rotator
from matplotlib.projections.geo import GeoAxes

###### WARNING #################
# this module is work in progress, the aim is to reimplement the healpy
# plot functions using the new features of matplotlib and remove most
# of the custom projection code

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi

    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return super(ThetaFormatterShiftPi, self).__call__(x, pos)

def lonlat(theta, phi):
    """Converts theta and phi to longitude and latitude

    From colatitude to latitude and from astro longitude to geo longitude"""

    longitude = -1*np.asarray(phi)
    latitude = np.pi/2 - np.asarray(theta)
    return longitude, latitude 

def mollview(m=None, rot=None, coord=None, unit='',
             xsize=1000, nest=False,
             min=None, max=None, flip='astro',
             format='%g',
             cbar=True, cmap=None,
             norm=None, 
             graticule=False, graticule_labels=False,
             **kwargs):
    """Plot an healpix map (given as an array) in Mollweide projection.
    
    Parameters
    ----------
    map : float, array-like or None
      An array containing the map, supports masked maps, see the `ma` function.
      If None, will display a blank map, useful for overplotting.
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, optional
      The minimum range value
    max : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    format : str, optional
      The format of the scale label. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    kwargs : keywords
      any additional keyword is passed to pcolormesh
    graticule : bool
      add graticule
    graticule_labels : bool
      longitude and latitude labels
    """

    # not implemented features
    if not (norm is None):
        raise NotImplementedError()

    # Create the figure
    import matplotlib.pyplot as plt

    width = 8.5
    fig = plt.figure(figsize=(width,width*.63))
    ax = fig.add_subplot(111, projection="mollweide")
    # FIXME: make a more general axes creation that works also with subplots
    #ax = plt.gcf().add_axes((.125, .1, .9, .9), projection="mollweide")

    # remove white space around the image
    plt.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.05)
    if graticule and graticule_labels:
        plt.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.05)

    if not m is None:
        # auto min and max
        if min is None:
            min = m.min()
        if max is None:
            max = m.max()

    # allow callers to override the hold state by passing hold=True|False
    washold = ax.ishold()
    hold = kwargs.pop('hold', None)
    if hold is not None:
        ax.hold(hold)

    try:
        ysize = xsize/2
        theta = np.linspace(np.pi, 0, ysize)
        phi   = np.linspace(-np.pi, np.pi, xsize)

        longitude = np.radians(np.linspace(-180, 180, xsize))
        if flip == "astro":
            longitude = longitude[::-1]
        latitude = np.radians(np.linspace(-90, 90, ysize))
        # project the map to a rectangular matrix xsize x ysize
        PHI, THETA = np.meshgrid(phi, theta)
        # coord or rotation
        if coord or rot:
            r = Rotator(coord=coord, rot=rot, inv=True)
            THETA, PHI = r(THETA.flatten(), PHI.flatten())
            THETA = THETA.reshape(ysize, xsize)
            PHI = PHI.reshape(ysize, xsize)
        nside = npix2nside(len(m))
        if not m is None:
            grid_pix = ang2pix(nside, THETA, PHI, nest=nest)
            grid_map = m[grid_pix]

            # plot
            ret = plt.pcolormesh(longitude, latitude, grid_map, vmin=min, vmax=max, rasterized=True, **kwargs)

        # graticule
        plt.grid(graticule)
        if graticule:
            longitude_grid_spacing = 60 # deg
            ax.set_longitude_grid(longitude_grid_spacing)
            if width < 10:
                ax.set_latitude_grid(45)
                ax.set_longitude_grid_ends(90)

        if graticule_labels:
            ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(longitude_grid_spacing))
        else:
            # remove longitude and latitude labels
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])

        # colorbar
        if cbar and not m is None:
            cb = fig.colorbar(ret, orientation='horizontal', shrink=.4, pad=0.05, ticks=[min, max])
            cb.ax.xaxis.set_label_text(unit)
            cb.ax.xaxis.labelpad = -8
            # workaround for issue with viewers, see colorbar docstring
            cb.solids.set_edgecolor("face")

        plt.draw()
    finally:
        ax.hold(washold)

    return ret
    
def projplot(theta, phi, fmt=None, **kwargs):
    """projplot is a wrapper around :func:`matplotlib.Axes.plot` to take into account the
    spherical projection.

    You can call this function as::
    
       projplot(theta, phi)        # plot a line going through points at coord (theta, phi)
       projplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)
    
    Parameters
    ----------
    theta, phi : float, array-like
      Coordinates of point to plot in radians.
    fmt : str
      A format string (see :func:`matplotlib.Axes.plot` for details)

    Notes
    -----
    Other keywords are passed to :func:`matplotlib.Axes.plot`.

    See Also
    --------
    projscatter, projtext
    """
    import matplotlib.pyplot as plt
    longitude, latitude = lonlat(theta, phi)
    if fmt is None:
        ret = plt.plot(longitude, latitude, **kwargs)
    else:
        ret = plt.plot(longitude, latitude, fmt, **kwargs)
    return ret

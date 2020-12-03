__all__ = ["projview", "projplot"]

import numpy as np
from .pixelfunc import ang2pix, npix2nside
from .rotator import Rotator
from matplotlib.projections.geo import GeoAxes
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

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
            x += 2 * np.pi
        return super(ThetaFormatterShiftPi, self).__call__(x, pos)

class ThetaFormatterTheta(GeoAxes.ThetaFormatter):
    """Change thetaticks in rads to degs"""

    def __call__(self, x, pos=None):
        return super(ThetaFormatterTheta, self).__call__(x, pos)


def lonlat(theta, phi):
    """Converts theta and phi to longitude and latitude

    From colatitude to latitude and from astro longitude to geo longitude"""

    longitude = -1 * np.asarray(phi)
    latitude = np.pi / 2 - np.asarray(theta)
    return longitude, latitude


def projview(
    m=None,
    rot=None,
    coord=None,
    unit="",
    xsize=1000,
    nest=False,
    min=None,
    max=None,
    flip="astro",
    format="%g",
    cbar=True,
    cmap='jet',
    norm=None,
    graticule=False,
    graticule_labels=False,
    return_only_data=False,
    projection_type = "mollweide",
    cb_orientation = "horizontal",
    xlabel=None,
    ylabel=None,
    longitude_grid_spacing = 60,
    latitude_grid_spacing = 30,
    override_plot_properties=None,
    title=None,
    lcolor='black',
    **kwargs
):
    """Plot a healpix map (given as an array) in the chosen projection.
    
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
    projection_type :  {'aitoff', 'hammer', 'lambert', 'mollweide', 'cart', '3d', 'polar'}
      type of the plot
    cb_orientation : {'horizontal', 'vertical'}
      color bar orientation
    xlabel : str
      set x axis label
    ylabel : str
      set y axis label
    longitude_grid_spacing : float
      set x axis grid spacing
    latitude_grid_spacing : float
      set y axis grid spacing
    override_plot_properties : list
      Override the following plot proporties: shrink (colorbar size), pad (colorbar padding), 
      lpad (colorbar label padding), width of the figure, ratio of the height and width of the plot
    title : str
      set title of the plot
    lcolor : str
      change the color of the longitude tick labels, some color maps make it hard to read black tick labels
    """
    if not m is None:
    # auto min and max
        if min is None:
            min = m.min()
        if max is None:
            max = m.max()

    # do this to find how many decimals are in the colorbar labels, so that the padding in the vertical cbar can done properly
    def find_number_of_decimals(number):
        try:
            return len(str(number).split('.')[1])
        except:
            return 0
    
    decs = np.max([find_number_of_decimals(min),find_number_of_decimals(max)])
    if decs >= 3:
        lpad=-27
    else:
        lpad = -9*decs

    # default plot settings
    r = 0.63
    if projection_type=="3d":
        if cb_orientation=='vertical':
            shrink=0.55;pad=0.02;lpad=lpad;width=11.5
        if cb_orientation=='horizontal':
            shrink=0.2;pad=0;lpad=-10;width=14
    if projection_type == "aitoff" or projection_type ==  "hammer" or  projection_type == "lambert" or projection_type == "mollweide":
        if cb_orientation=='vertical':
            shrink=0.6;pad=0.01;lpad=lpad;width=10
        if cb_orientation=='horizontal':
            shrink=0.6;pad=0.05;lpad=-8;width=8.5
    if projection_type =="cart":
        if cb_orientation=='vertical':
           shrink=1;pad=0.01;lpad=lpad;width=9.6 ;r=0.42;
        if  cb_orientation=='horizontal':
            shrink=0.4;pad=0.1;lpad=-12;width=8.8
            if xlabel == None:
                pad=0.01;r=0.63
    if projection_type=="polar":
        if  cb_orientation=='vertical':
            shrink=1;pad=0.01;lpad=lpad;width=10
        if cb_orientation=='horizontal':
            shrink=0.4;pad=0.01;lpad=0;width=12
    if override_plot_properties is not None:
        print('Overriding default plot properies: shrink= {}, pad={}, lpad={},width={},r={}, to one:  shrink= {}, pad={}, lpad={},width={},r={}'.format(\
        shrink,pad,lpad,width,r,override_plot_properties[0],override_plot_properties[1],override_plot_properties[2],override_plot_properties[3],override_plot_properties[4]))
        shrink=override_plot_properties[0]
        pad=override_plot_properties[1]
        lpad=override_plot_properties[2]
        width=override_plot_properties[3]
        r=override_plot_properties[4]

    # not implemented features
    if not (norm is None):
        raise NotImplementedError()

    # Create the figure
    if not return_only_data: # supress figure creation when only dumping the data
        import matplotlib.pyplot as plt
    
        width = width # 8.5
        fig = plt.figure(figsize=(width, width * r))
        print(str(projection_type))
        if projection_type == "cart":
            ax = fig.add_subplot(111)
        else:
            ax = fig.add_subplot(111, projection=projection_type)
        # FIXME: make a more general axes creation that works also with subplots
        # ax = plt.gcf().add_axes((.125, .1, .9, .9), projection="mollweide")
    
        # remove white space around the image
        plt.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.05)
    #end if not
    if graticule and graticule_labels:
        plt.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.05)

    # allow callers to override the hold state by passing hold=True|False
    #washold = ax.ishold() #  commented out
    hold = kwargs.pop("hold", None)
    #if hold is not None:
    #    ax.hold(hold)

    try:
        ysize = int(xsize / 2)
        theta = np.linspace(np.pi, 0, ysize)
        phi = np.linspace(-np.pi, np.pi, xsize)

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
            if return_only_data: # exit here when dumping the data
                return [longitude,  latitude, grid_map]
            if projection_type is not '3d': # test for 3d plot
                ret = plt.pcolormesh(
                    longitude,
                    latitude,
                    grid_map,
                    vmin=min,
                    vmax=max,
                    rasterized=True,
                    cmap=cmap,
                    shading='auto',
                    **kwargs
                )
            elif projection_type is '3d': # test for 3d plot
                LONGITUDE, LATITUDE = np.meshgrid(longitude, latitude)
                ret = ax.plot_surface(LONGITUDE, LATITUDE, grid_map, cmap=cmap, vmin=min,  vmax=max, rasterized=True, **kwargs)
        # graticule
        plt.grid(graticule)
        
        if graticule:
            if projection_type == "aitoff" or projection_type ==  "hammer" or  projection_type == "lambert" or projection_type == "mollweide":
                longitude_grid_spacing = longitude_grid_spacing  # deg 60
                ax.set_longitude_grid(longitude_grid_spacing)
                ax.set_latitude_grid(latitude_grid_spacing)
                ax.set_longitude_grid_ends(90)
         #       if width < 10:
          #          ax.set_latitude_grid(45)
           #         ax.set_longitude_grid_ends(90)
            else:
                longitude_grid_spacing = longitude_grid_spacing  # deg
                latitude_grid_spacing = latitude_grid_spacing #  deg
                ax.xaxis.set_major_locator(MultipleLocator(np.deg2rad(longitude_grid_spacing) )) # longitude
                ax.yaxis.set_major_locator(MultipleLocator(np.deg2rad(latitude_grid_spacing))) # lattitude

        if graticule_labels & graticule:
            if projection_type == "aitoff" or projection_type ==  "hammer" or  projection_type == "lambert" or projection_type == "mollweide":
                ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(longitude_grid_spacing))
                ax.tick_params(axis='x', colors=lcolor)
            else:
                ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(longitude_grid_spacing))
                ax.yaxis.set_major_formatter(ThetaFormatterTheta(latitude_grid_spacing))
        if graticule == False:
            # remove longitude and latitude labels
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.tick_params(axis=u'both', which=u'both',length=0)
        
        ax.set_title(title, fontsize= 14)
        # colorbar
        if projection_type is "cart":
            ax.set_aspect(1)
        extend = 'neither'
        if min > np.min(m):
            extend = "min"
        if max < np.max(m):
            extend = "max"
        if min > np.min(m) and max < np.max(m):
            extend = 'both'
        cb = fig.colorbar(
            ret, orientation=cb_orientation, shrink=shrink, pad=pad, ticks=[min, max], extend=extend
        )
        if cb_orientation == 'horizontal':
            cb.ax.xaxis.set_label_text(unit)
            cb.ax.xaxis.labelpad = lpad
        if cb_orientation == 'vertical':
            cb.ax.yaxis.set_label_text(unit)
            cb.ax.yaxis.labelpad = lpad
        # workaround for issue with viewers, see colorbar docstring
        cb.solids.set_edgecolor("face")
        ax.set_xlabel(xlabel,fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)
        plt.draw()
    finally:
        #ax.hold(washold)
        None

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

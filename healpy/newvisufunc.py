__all__ = ["projview", "newprojplot"]

import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections.geo import GeoAxes
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import Colormap, Normalize
from matplotlib import colors
from typing import Union, Sequence
from .rotator import Rotator
from .pixelfunc import ang2pix, npix2nside


class ThetaFormatterCounterclockwisePhi(GeoAxes.ThetaFormatter):
    """Convert tick labels from rads to degs and shifts labelling from -180|-90|0|90|180 to conterclockwise periodic 180|90|0|270|180 """

    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2 * np.pi
        return super(ThetaFormatterCounterclockwisePhi, self).__call__(x, pos)


class ThetaFormatterClockwisePhi(GeoAxes.ThetaFormatter):
    """Convert tick labels from rads to degs and shifts labelling from -180|-90|0|90|180 to clockwise periodic 180|270|0|90|180 """

    def __call__(self, x, pos=None):

        if x < 0:
            x += 2 * np.pi
        #   return super(ThetaFormatterShiftPhi, self).__call__(x, pos)
        return super(ThetaFormatterClockwisePhi, self).__call__(x, pos)


class ThetaFormatterSymmetricPhi(GeoAxes.ThetaFormatter):
    """Just convert phi ticks from rad to degs and keep the true -180|-90|0|90|180 """

    def __call__(self, x, pos=None):
        return super(ThetaFormatterSymmetricPhi, self).__call__(x, pos)


class ThetaFormatterTheta(GeoAxes.ThetaFormatter):
    """Convert theta ticks from rads to degs"""

    def __call__(self, x, pos=None):
        return super(ThetaFormatterTheta, self).__call__(x, pos)


def lonlat(theta, phi):
    """Converts theta and phi to longitude and latitude"""

    longitude = np.asarray(phi)
    latitude = np.pi / 2 - np.asarray(theta)
    return longitude, latitude


def update_dictionary(main_dict, update_dict):
    for key, key_val in main_dict.items():
        if key in update_dict:
            main_dict[key] = update_dict[key]
    return main_dict

class DummyNorm(Normalize):
    """ A template to generate a new Matplotlib Norm transformation

    Basically you only need to derive this class and set the
    `_transform` and `_inverse_transform` methods.
    """
    def _transform(self, value, out=None):
        delta = self.vmax - self.vmin
        return (value - self.vmin) / float(delta)

    def _inverse_transform(self, value, out=None):
        delta = self.vmax - self.vmin
        return delta * value + self.vmin

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        # Normalize initial recipe to clean inputs
        if clip:
            result, is_scalar = self.process_value(np.clip(value),
                                                   self.vmin,
                                                   self.vmax)
        else:
            result, is_scalar = self.process_value(value)

        resdat = result.data
        resdat = self._transform(resdat)
        result = np.ma.array(resdat, mask=result.mask, copy=False)
        if is_scalar:
            return result[0]
        return result

    def inverse(self, value):
        result, is_scalar = self.process_value(value)
        resdat = result.data
        resdat = self._inverse_transform(resdat)
        result = np.ma.array(resdat, mask=result.mask, copy=False)
        if is_scalar:
            return result[0]
        return result


class HistEq(DummyNorm):
    """
    Histogram equalizer normalization

    Attributes
    ----------
    artist: Artist instance or ndarray
        artist or image to normalize
    bins: int, optional
        number of bins to calculate the histogram
    """
    def __init__(self, artist, *args, **kwargs):
        bins = kwargs.pop('bins', 1024)
        super().__init__(*args, **kwargs)
        try:
            image = artist._A
        except AttributeError:
            image = artist
        values, bins = np.histogram(image.ravel(), bins=bins)
        centers = 0.5 * (bins[:-1] + bins[1:])
        csum = values.cumsum()
        self.histogram = centers, csum / csum[-1], values

    def _transform(self, value, out=None):
        return np.interp(value, self.histogram[0], self.histogram[1])

    def _inverse_transform(self, value, out=None):
        return np.interp(value, self.histogram[1], self.histogram[0])


class Arcsinh(DummyNorm):
    """ Arcsinh normalization """
    def _transform(self, value, out=None):
        delta = np.arcsinh(self.vmax) - np.arcsinh(self.vmin)
        return (np.arcsinh(value) - np.arcsinh(self.vmin)) / float(delta)

    def _inverse_transform(self, value, out=None):
        delta = np.arcsinh(self.vmax) - np.arcsinh(self.vmin)
        return np.sinh(delta * value + np.arcsinh(self.vmin))


class MidpointNormalize(Normalize):
    """ Normalize to a central point """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


class Sqrt(DummyNorm):
    """ Sqrt normalization """
    def _transform(self, value, out=None):
        delta = np.sqrt(self.vmax) - np.sqrt(self.vmin)
        return (np.sqrt(value) - np.sqrt(self.vmin)) / float(delta)

    def _inverse_transform(self, value, out=None):
        delta = np.power(self.vmax, 2) - np.power(self.vmin, 2)
        return np.power(delta * value + np.power(self.vmin, 2))


class Power(DummyNorm):
    """ Power normalization """
    def __init__(self, power_value=2, vmin=None, vmax=None, clip=False):
        self.power_ = power_value
        self.inv_power_ = 1. / power_value
        Normalize.__init__(self, vmin, vmax, clip)

    def _transform(self, value, out=None):
        delta = np.power(self.vmax, self.power_) - np.power(self.vmin, self.power_)
        return (np.power(value, self.power_) - np.power(self.vmin, self.power_)) / float(delta)

    def _inverse_transform(self, value, out=None):
        delta = np.power(self.vmax, self.inv_power_) - np.power(self.vmin, self.inv_power_)
        return np.power(delta * value + np.power(self.vmin, self.inv_power_))


def parse_norm(norm, data=None, **kwargs):
        """ Allows one to use a string shortcut to defined the norm keyword
        e.g.: parse_norm('log10')
        Parameters
        ----------
        kwargs: dict
            keywords given to the artist
        Returns
        -------
        norm: Normalize object
        """

        mapped = {'arcsinh': Arcsinh,
                  'log': colors.LogNorm,
                  'log10': colors.LogNorm,
                  'sqrt': Sqrt,
                  'pow': Power,
                  'histeq': HistEq,
                  'midpoint': MidpointNormalize,
                  }
        if norm in (None, 'None', 'none', ''):
            return Normalize(**kwargs)
        try:
            norm_ = eval(norm, mapped, colors.__dict__)
            if isinstance(norm_, type):
                if norm_ == HistEq:
                    return norm_(data, **kwargs)
                return norm_(**kwargs)
            return norm_
        except Exception:
            return norm


def projview(
    m: Union[float, np.array] = None,
    rot: Sequence[float] = None,
    coord: Sequence[str] = None,
    unit: str = "",
    xsize: int = 1000,
    nest: bool = False,
    norm: Union[str, Normalize] = None,
    vmin: float = None,
    vmax: float = None,
    flip: str = "astro",
    cbar: bool = True,
    cmap: Union[str, Colormap] = "viridis",
    graticule: bool = False,
    graticule_labels: bool = False,
    return_only_data: bool = False,
    projection_type: str = "mollweide",
    cb_orientation: str = "horizontal",
    xlabel: str = None,
    ylabel: str = None,
    longitude_grid_spacing: float = 60,
    latitude_grid_spacing: float = 30,
    override_plot_properties: dict = None,
    title: str = None,
    xtick_label_color: str = "black",
    ytick_label_color: str = "black",
    graticule_color: str = None,
    fontsize: dict = None,
    phi_convention: str = "counterclockwise",
    custom_xtick_labels: Sequence[str] = None,
    custom_ytick_labels: Sequence[str] = None,
    **kwargs
) -> Union[Sequence, None]:
    """Plot a healpix map (given as an array) in the chosen projection.

    See examples of using this function in the documentation under "Other tutorials".
    Overplot points or lines using :func:`newprojplot`.

    .. warning::
        this function is work in progress, the aim is to reimplement the healpy
        plot functions using the new features of matplotlib and remove most
        of the custom projection code.
        Please report bugs or submit feature requests via Github.
        The interface will change in future releases.

    Parameters
    ----------
    m : float, array-like or None
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
    vmin : float, optional
      The minimum range value
    vmax : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
      It creates the `healpy_flip` attribute on the Axes to save the convention in the figure.
    format : str, optional
      The format of the scale label. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    norm : str, colors.Normalize or None, optional
      Color normalization,
      available values: 'arcsinh', 'log', 'log10', 'sqrt', 'pow', 'histeq', 'midpoint'
      (log and log10 are base 10)
      any colors.Normalize object is also valid
    kwargs : keywords
      any additional keyword is passed to pcolormesh
    graticule : bool
      add graticule
    graticule_labels : bool
      longitude and latitude labels
    projection_type :  str
     valid values in
        {'aitoff', 'hammer', 'lambert', 'mollweide', 'cart', '3d', 'polar'}
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
    override_plot_properties : dict
      Override the following plot properties: "cbar_shrink", "cbar_pad",
      "cbar_label_pad", "figure_width": width, "figure_size_ratio": ratio.
    title : str
      set title of the plot
    lcolor : str
      change the color of the longitude tick labels, some color maps make it
      hard to read black tick labels
    fontsize:  dict
        Override fontsize of labels: "xlabel", "ylabel", "title", "xtick_label",
        "ytick_label", "cbar_label", "cbar_tick_label".
    phi_convention : string
        convention on x-axis (phi), 'counterclockwise' (default), 'clockwise',
        'symmetrical' (phi as it is truly given) if `flip` is "geo",
        `phi_convention` should be set to 'clockwise'.
    custom_xtick_labels : list
        override x-axis tick labels
    custom_ytick_labels : list
        override y-axis tick labels
    """
    geographic_projections = ["aitoff", "hammer", "lambert", "mollweide"]

    if not m is None:
        # auto min and max
        if vmin is None:
            vmin = m.min()
        if max is None:
            vmax = m.max()

    # do this to find how many decimals are in the colorbar labels,
    # so that the padding in the vertical cbar can done properly
    def find_number_of_decimals(number):
        try:
            return len(str(number).split(".")[1])
        except:
            return 0

    # default font sizes
    fontsize_defaults = {
        "xlabel": 12,
        "ylabel": 12,
        "title": 14,
        "xtick_label": 12,
        "ytick_label": 12,
        "cbar_label": 12,
        "cbar_tick_label": 12,
    }
    if fontsize is not None:
        fontsize_defaults.update(fontsize_defaults)

    # default plot settings
    decs = np.max([find_number_of_decimals(vmin), find_number_of_decimals(vmax)])
    if decs >= 3:
        lpad = -27
    else:
        lpad = -9 * decs
    ratio = 0.63
    if projection_type == "3d":
        if cb_orientation == "vertical":
            shrink = 0.55
            pad = 0.02
            lpad = lpad
            width = 11.5
        if cb_orientation == "horizontal":
            shrink = 0.2
            pad = 0
            lpad = -10
            width = 14
    if projection_type in geographic_projections:
        if cb_orientation == "vertical":
            shrink = 0.6
            pad = 0.01
            lpad = lpad
            width = 10
        if cb_orientation == "horizontal":
            shrink = 0.6
            pad = 0.05
            lpad = -8
            width = 8.5
    if projection_type == "cart":
        if cb_orientation == "vertical":
            shrink = 1
            pad = 0.01
            lpad = lpad
            width = 9.6
            ratio = 0.42
        if cb_orientation == "horizontal":
            shrink = 0.4
            pad = 0.1
            lpad = -12
            width = 8.8
            if xlabel == None:
                pad = 0.01
                ratio = 0.63
    if projection_type == "polar":
        if cb_orientation == "vertical":
            shrink = 1
            pad = 0.01
            lpad = lpad
            width = 10
        if cb_orientation == "horizontal":
            shrink = 0.4
            pad = 0.01
            lpad = 0
            width = 12
    # pass the default settings to the plot_properties dictionary
    plot_properties = {
        "cbar_shrink": shrink,
        "cbar_pad": pad,
        "cbar_label_pad": lpad,
        "figure_width": width,
        "figure_size_ratio": ratio,
    }

    if override_plot_properties is not None:
        warnings.warn(
            "\n *** Overriding default plot properies: " + str(plot_properties) + " ***"
        )
        plot_properties.update(override_plot_properties)
        warnings.warn("\n *** New plot properies: " + str(plot_properties) + " ***")

    ysize = xsize // 2
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

    if m is not None:
        nside = npix2nside(len(m))
        grid_pix = ang2pix(nside, THETA, PHI, nest=nest)
        grid_map = m[grid_pix]

        # plot
        if return_only_data:  # exit here when dumping the data
            return [longitude, latitude, grid_map]

        # Create the figure

        width = width  # 8.5
        fig = plt.figure(
            figsize=(
                plot_properties["figure_width"],
                plot_properties["figure_width"] * plot_properties["figure_size_ratio"],
            )
        )

        if projection_type == "cart":
            ax = fig.add_subplot(111)
        else:
            ax = fig.add_subplot(111, projection=projection_type)
        # FIXME: make a more general axes creation that works also with subplots
        # ax = plt.gcf().add_axes((.125, .1, .9, .9), projection="mollweide")

        # set property on ax so it can be used in newprojplot
        ax.healpy_flip = flip

        if graticule and graticule_labels:
            plt.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.05)
    # end if not

        norm_ = parse_norm(norm, data=grid_map, vmin=vmin, vmax=vmax)

        if projection_type != "3d":  # test for 3d plot
            ret = plt.pcolormesh(
                longitude,
                latitude,
                grid_map,
                norm=norm_,
                rasterized=True,
                cmap=cmap,
                shading="auto",
                **kwargs
            )
        elif projection_type == "3d":  # test for 3d plot
            LONGITUDE, LATITUDE = np.meshgrid(longitude, latitude)
            ret = ax.plot_surface(
                LONGITUDE,
                LATITUDE,
                grid_map,
                cmap=cmap,
                norm=norm_,
                rasterized=True,
                **kwargs
            )
    # graticule
    if graticule_color is None:
        plt.grid(graticule)
    else:
        plt.grid(graticule, color=graticule_color)

    if graticule:
        if projection_type in geographic_projections:
            longitude_grid_spacing = longitude_grid_spacing  # deg 60
            ax.set_longitude_grid(longitude_grid_spacing)
            ax.set_latitude_grid(latitude_grid_spacing)
            ax.set_longitude_grid_ends(90)
        else:
            longitude_grid_spacing = longitude_grid_spacing  # deg
            latitude_grid_spacing = latitude_grid_spacing  #  deg
            ax.xaxis.set_major_locator(
                MultipleLocator(np.deg2rad(longitude_grid_spacing))
            )  # longitude
            ax.yaxis.set_major_locator(
                MultipleLocator(np.deg2rad(latitude_grid_spacing))
            )  # lattitude

    # labelling
    if graticule_labels & graticule:
        if phi_convention == "counterclockwise":
            xtick_formatter = ThetaFormatterCounterclockwisePhi(longitude_grid_spacing)
        elif phi_convention == "clockwise":
            xtick_formatter = ThetaFormatterClockwisePhi(longitude_grid_spacing)
        elif phi_convention == "symmetrical":
            xtick_formatter = ThetaFormatterSymmetricPhi(longitude_grid_spacing)

        ax.xaxis.set_major_formatter(xtick_formatter)
        ax.yaxis.set_major_formatter(ThetaFormatterTheta(latitude_grid_spacing))

        if custom_xtick_labels is not None:
            try:
                ax.xaxis.set_ticklabels(custom_xtick_labels)
            except:
                warnings.warn(
                    "Put names for all "
                    + str(len(ax.xaxis.get_ticklabels()))
                    + " x-tick labels!. No re-labelling done."
                )
        if custom_ytick_labels is not None:
            try:
                ax.yaxis.set_ticklabels(custom_ytick_labels)
            except:
                warnings.warn(
                    "Put names for all "
                    + str(len(ax.yaxis.get_ticklabels()))
                    + " y-tick labels!. No re-labelling done."
                )
    if not graticule:
        # remove longitude and latitude labels
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.tick_params(axis=u"both", which=u"both", length=0)

    ax.set_title(title, fontsize=fontsize_defaults["title"])
    # tick font size
    ax.tick_params(
        axis="x", labelsize=fontsize_defaults["xtick_label"], colors=xtick_label_color
    )
    ax.tick_params(
        axis="y", labelsize=fontsize_defaults["ytick_label"], colors=ytick_label_color
    )

    # colorbar
    if cbar:
        if projection_type == "cart":
            ax.set_aspect(1)
        extend = "neither"

        if norm_.vmin > np.min(m):
            extend = "min"
        if norm_.vmax < np.max(m):
            extend = "max"
        if norm_.vmin > np.min(m) and norm_.vmax < np.max(m):
            extend = "both"

        cb = fig.colorbar(
            ret,
            orientation=cb_orientation,
            shrink=plot_properties["cbar_shrink"],
            pad=plot_properties["cbar_pad"],
            extend=extend,
        )

        if cb_orientation == "horizontal":
            cb.ax.xaxis.set_label_text(unit, fontsize=fontsize_defaults["cbar_label"])
            cb.ax.tick_params(axis="x", labelsize=fontsize_defaults["cbar_tick_label"])
            cb.ax.xaxis.labelpad = plot_properties["cbar_label_pad"]
        if cb_orientation == "vertical":
            cb.ax.yaxis.set_label_text(unit, fontsize=fontsize_defaults["cbar_label"])
            cb.ax.tick_params(axis="y", labelsize=fontsize_defaults["cbar_tick_label"])
            cb.ax.yaxis.labelpad = plot_properties["cbar_label_pad"]
        # workaround for issue with viewers, see colorbar docstring
        cb.solids.set_edgecolor("face")
    ax.set_xlabel(xlabel, fontsize=fontsize_defaults["xlabel"])
    ax.set_ylabel(ylabel, fontsize=fontsize_defaults["ylabel"])

    return ret


def newprojplot(theta, phi, fmt=None, lonlat=True, ax=None, **kwargs):
    """newprojplot is a wrapper around :func:`matplotlib.Axes.plot` to support
    colatitude theta and longitude phi and take into account the longitude convention
    (see the `flip` keyword of :func:`projview`)

    You can call this function as::

       newprojplot(theta, phi)        # plot a line going through points at coord (theta, phi)
       newprojplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)

    Parameters
    ----------
    theta, phi : float, array-like
      Coordinates of point to plot in radians.
    fmt : str
      A format string (see :func:`matplotlib.Axes.plot` for details)
    lonlat : bool
        Assume theta and phi to longitude and latitude in degrees.

    Notes
    -----
    Other keywords are passed to :func:`matplotlib.Axes.plot`.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    flip = getattr(ax, "healpy_flip", "astro")

    if lonlat:
        longitude, latitude = np.deg2rad(np.array(theta)), np.deg2rad(np.array(phi))
    else:
        longitude, latitude = lonlat(theta, phi)

    if flip == "astro":
        longitude = longitude * -1
    if fmt is None:
        ret = plt.plot(longitude, latitude, **kwargs)
    else:
        ret = plt.plot(longitude, latitude, fmt, **kwargs)
    return ret


def newprojtext(theta, phi, text, lonlat=True, ax=None, **kwargs):
    """newprojplot is a wrapper around :func:`matplotlib.Axes.plot` to support
    colatitude theta and longitude phi and take into account the longitude convention
    (see the `flip` keyword of :func:`projview`)

    You can call this function as::

       newprojplot(theta, phi)        # plot a line going through points at coord (theta, phi)
       newprojplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)

    Parameters
    ----------
    theta, phi : float, array-like
      Coordinates of point to plot in radians.
    text: str
        text to show
    lonlat : bool
        Assume theta and phi to longitude and latitude in degrees.

    Notes
    -----
    Other keywords are passed to :func:`matplotlib.Axes.plot`.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()
    flip = getattr(ax, "healpy_flip", "astro")

    if lonlat:
        longitude, latitude = np.deg2rad(np.array(theta)), np.deg2rad(np.array(phi))
    else:
        longitude, latitude = lonlat(theta, phi)
    if flip == "astro":
        longitude = longitude * -1
    return ax.text(longitude, latitude, text, **kwargs)


def CreateRotatedGraticule(rot, t_step=30, p_step=30):

    phi = rot[0]
    try:
        theta = rot[1]
    except:
        theta = 0
    pointDensity = 100
    conventionThetaOffset = np.pi / 2
    thetaSpacing = np.arange(-90, 90 + t_step, t_step)
    phiSpacing = np.arange(-180, 180 + p_step, p_step)
    where_zero = np.hstack(
        (
            np.where(thetaSpacing == 0)[0],
            thetaSpacing.size + np.where(phiSpacing == 0)[0],
        )
    )
    gline_phi_fixed = np.deg2rad(np.linspace(-180, 180, pointDensity))
    gline_theta_fixed = (
        np.deg2rad(np.linspace(-90, 90, pointDensity)) + conventionThetaOffset
    )
    rotated_grid_lines = []

    for thetaSpace in thetaSpacing:
        gline_theta = (
            np.deg2rad(np.zeros(pointDensity) + thetaSpace) + conventionThetaOffset
        )
        r = Rotator(rot=(0, theta), inv=True)
        gline_theta_rot, gline_phi_fixed_rot = r(gline_theta, gline_phi_fixed)
        gline_theta = gline_theta - conventionThetaOffset
        gline_theta_rot = gline_theta_rot - conventionThetaOffset
        rotated_grid_lines.append([gline_phi_fixed_rot, gline_theta_rot])

    for phiSpace in phiSpacing:
        gline_phi = np.deg2rad(np.zeros(pointDensity) + phiSpace)
        r = Rotator(rot=(phi, 0), inv=True)
        gline_theta_fixed_rot, gline_phi_rot = r(gline_theta_fixed, gline_phi)
        r = Rotator(rot=(0, theta), inv=True)
        gline_theta_fixed_rot, gline_phi_rot = r(gline_theta_fixed_rot, gline_phi_rot)
        gline_theta_fixed_rot = gline_theta_fixed_rot - conventionThetaOffset
        rotated_grid_lines.append([gline_phi_rot, gline_theta_fixed_rot])

    for g_lines in rotated_grid_lines:
        mask = np.where((np.abs(np.diff(g_lines[0]))) > np.deg2rad(45))
        g_lines[0] = np.ma.array(g_lines[0])
        g_lines[0][mask] = np.ma.masked

    return rotated_grid_lines, where_zero

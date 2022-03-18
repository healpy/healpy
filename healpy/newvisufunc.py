__all__ = ["projview", "newprojplot"]

from curses import cbreak
import numpy as np
from .pixelfunc import ang2pix, npix2nside, remove_dipole, remove_monopole
from .rotator import Rotator, coordsys2euler_zyz, vec2dir
import matplotlib.pyplot as plt
from matplotlib.projections.geo import GeoAxes
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoLocator
import warnings


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


def projview(
    m=None,
    fig=None,
    rot=None,
    coord=None,
    unit="",
    xsize=1000,
    nest=False,
    min=None,
    max=None,
    ticks=[None,None],
    flip="astro",
    format="%g",
    cbar=True,
    cmap="viridis",
    norm=None,
    graticule=False,
    graticule_labels=False,
    rot_graticule=False,
    graticule_coord = None,
    override_rot_graticule_properties=None,
    return_only_data=False,
    projection_type="mollweide",
    cb_orientation="horizontal",
    xlabel=None,
    ylabel=None,
    longitude_grid_spacing=60,
    latitude_grid_spacing=30,
    override_plot_properties=None,
    title=None,
    rlabel=None,
    llabel=None,
    xtick_label_color="black",
    ytick_label_color="black",
    graticule_color=None,
    serif=True,
    fontsize=None,
    phi_convention="counterclockwise",
    custom_xtick_labels=None,
    custom_ytick_labels=None,
    invRot=True,
    sub=111,
    reuse_axes=False,
    margins= None,
    hold=False,
    remove_dip=False,
    remove_mono=False,
    gal_cut=0,
    **kwargs
):
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
    fig : int or None, optional
      The figure number to use. Default: create a new figure
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
      default: 'G'
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
    ticks : sequence of int, optional
      Tick values for colorbar. 
      Overwrites min and max in colorbar.
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
      It creates the `healpy_flip` attribute on the Axes to save the convention in the figure.
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
    rot_graticule : bool
      rotate also the graticule when rotating the map
    graticule_coord : str
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the graticule
    override_rot_graticule_properties : dict
      Override the following rotated graticule properties: "g_linestyle", "g_linewidth", "g_color", 
      "g_alpha", "t_step", "p_step".
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
    override_plot_properties : dict
      Override the following plot properties: "cbar_shrink", "cbar_pad", "cbar_label_pad", "cbar_tick_direction"
      "figure_width": width, "figure_size_ratio": ratio.
    title : str
      set title of the plot
    rlabel : str
      set label at top right corner of axis
    llabel : str
      set label at top left corner of axis
    lcolor : str
      change the color of the longitude tick labels, some color maps make it hard to read black tick labels
    serif : bool
      Set fontstyle to serif
      Default: True
    fontsize:  dict
      Override fontsize of labels: "xlabel", "ylabel", "title", "xtick_label", "ytick_label", 
      "cbar_label", "cbar_tick_label".
    phi_convention : string
      convention on x-axis (phi), 'counterclockwise' (default), 'clockwise', 'symmetrical' (phi as it is truly given)
      if `flip` is "geo", `phi_convention` should be set to 'clockwise'.
    custom_xtick_labels : list
      override x-axis tick labels
    custom_ytick_labels : list
      override y-axis tick labels
    invRot : bool
      invert rotation
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: 111
    reuse_axes : bool, optional
      If True, reuse the current Axes (should be a MollweideAxes). This is
      useful if you want to overplot with a partially transparent colormap,
      such as for plotting a line integral convolution. Default: False
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    hold : bool, optional
      If True, replace the current Axes by new axis.
      use this if you want to have multiple maps on the same
      figure. Default: False
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    """
    geographic_projections = ["aitoff", "hammer", "lambert", "mollweide"]
    if serif:
        plt.rc(
            "font",
            family="serif",
        )
        plt.rcParams["mathtext.fontset"] = "stix"
        plt.rc(
            "text.latex",
            preamble=r"\usepackage{sfmath}",
        )
    # If no min or max, set to ticks value. If ticks is None, values will be None
    vmin = min
    vmax = max
    if ticks != [None, None]:
        if min is None:
            vmin=ticks[0]
        if max is None:
            vmax=ticks[-1]
    if not m is None:
        if remove_dip:
            m = remove_dipole(
                m, gal_cut=gal_cut, nest=nest, copy=True
            )
        elif remove_mono:
            m = remove_monopole(
                m, gal_cut=gal_cut, nest=nest, copy=True
            )
            
        # auto min and max
        percentile = 97.5
        if vmin is None:
            vmin = np.percentile(m, 100.0 - percentile)
        if vmax is None:
            vmax = np.percentile(m, percentile)


    # do this to find how many decimals are in the colorbar labels, so that the padding in the vertical cbar can done properly
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
        "cbar_tick_label": 10,
    }
    if fontsize is not None:
        fontsize_defaults = update_dictionary(fontsize_defaults, fontsize)

    # default plot settings
    decs = np.max([find_number_of_decimals(vmin), find_number_of_decimals(vmax)])
    if decs >= 3:
        lpad = -27
    else:
        lpad = -9 * decs

    ratio = 0.63
    ratio = 0.5
    if projection_type == "3d":
        if cb_orientation == "vertical":
            shrink = 0.55
            pad = 0.02
            lpad = 4 #lpad
            width = 11.5
        if cb_orientation == "horizontal":
            shrink = 0.2
            pad = 0
            lpad = -10
            width = 14
    if projection_type in geographic_projections:
        if cb_orientation == "vertical":
            shrink = 0.7
            pad = 0.01
            lpad = 4  #lpad
            width = 10
        if cb_orientation == "horizontal":
            shrink = 0.4
            pad = 0.05
            lpad = 0
            width = 8.5
    if projection_type == "cart":
        if cb_orientation == "vertical":
            shrink = 1
            pad = 0.01
            lpad = 4 #lpad
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
            lpad = 4 #lpad
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
        "cbar_tick_direction": "in",
        "figure_width": width,
        "figure_size_ratio": ratio,
    }

    if override_plot_properties is not None:
        warnings.warn(
            "\n *** Overriding default plot properies: " + str(plot_properties) + " ***"
        )
        plot_properties = update_dictionary(plot_properties, override_plot_properties)
        warnings.warn("\n *** New plot properies: " + str(plot_properties) + " ***")

    rot_graticule_properties = {
        "g_linestyle": "-",
        "g_color": "w",
        "g_alpha": 0.75,
        "g_linewidth": 0.75,
        "t_step": 30,
        "p_step": 30,
    }

    if override_rot_graticule_properties is not None:
        warnings.warn(
            "\n *** Overriding rotated graticule properies: "
            + str(rot_graticule_properties)
            + " ***"
        )
        rot_graticule_properties = update_dictionary(
            rot_graticule_properties, override_rot_graticule_properties
        )
        warnings.warn(
            "\n *** New rotated graticule properies: "
            + str(rot_graticule_properties)
            + " ***"
        )

    # Create the figure
    if not return_only_data:  # supress figure creation when only dumping the data
        if not (hold or reuse_axes) and sub==111:
            fig = plt.figure(fig, figsize=(
                plot_properties["figure_width"],
                (plot_properties["figure_width"] * plot_properties["figure_size_ratio"]),
            ))
            extent = (0.02, 0.05, 0.96, 0.9)
        elif hold:
            fig = plt.gcf()
            left, bottom, right, top = np.array(fig.gca().get_position()).ravel()
            extent = (left, bottom, right - left, top - bottom)
            fig.delaxes(fig.gca())
        elif reuse_axes:
            fig = plt.gcf()
        else:  # using subplot syntax
            if hasattr(sub, "__len__"):
                nrows, ncols, idx = sub
            else:
                nrows, ncols, idx = sub // 100, (sub % 100) // 10, (sub % 10)
            if idx < 1 or idx > ncols * nrows:
                raise ValueError("Wrong values for sub: %d, %d, %d" % (nrows, ncols, idx))

            if not plt.get_fignums():
                # Scale height depending on subplots
                fig = plt.figure(fig, figsize=(
                    plot_properties["figure_width"],
                   (plot_properties["figure_width"]*plot_properties["figure_size_ratio"])*(nrows/ncols),
                ))
            else:
                fig = plt.gcf()

            """
            # Subplot method 1, copied from mollview
            c, r = (idx - 1) % ncols, (idx - 1) // ncols
            if not margins:
                right_adjust = 0.045 if cb_orientation=="vertical" else 0.0
                margins = (0.01, 0.0, 0.01-right_adjust, 0.0)

            extent = (
                c * 1.0 / ncols + margins[0],
                1.0 - (r + 1) * 1.0 / nrows + margins[1],
                1.0 / ncols - margins[2] - margins[0],
                1.0 / nrows - margins[3] - margins[1],
            )
            extent = (
                extent[0] + margins[0],
                extent[1] + margins[1],
                extent[2] - margins[2] - margins[0],
                extent[3] - margins[3] - margins[1],
            )
            """
        # FIXME: make a more general axes creation that works also with subplots
        #ax = fig.add_axes(extent, projection=projection_type)
        if projection_type == "cart":
            ax = fig.add_subplot(sub)
        else:
            ax = fig.add_subplot(sub, projection=projection_type)
        

    # Parameters for subplots
    left=0.02
    right=0.98
    top=0.95
    bottom=0.05

    # end if not
    if graticule and graticule_labels:
        left+=0.02
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom,)

    ysize = xsize // 2
    theta = np.linspace(np.pi, 0, ysize)
    phi = np.linspace(-np.pi, np.pi, xsize)

    longitude = np.radians(np.linspace(-180, 180, xsize))
    if flip == "astro":
        longitude = longitude[::-1]
    if not return_only_data:
        # set property on ax so it can be used in newprojplot
        ax.healpy_flip = flip

    latitude = np.radians(np.linspace(-90, 90, ysize))
    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    # coord or rotation
    if coord or rot:
        r = Rotator(coord=coord, rot=rot, inv=invRot)
        THETA, PHI = r(THETA.flatten(), PHI.flatten())
        THETA = THETA.reshape(ysize, xsize)
        PHI = PHI.reshape(ysize, xsize)
    nside = npix2nside(len(m))
    if not m is None:
        grid_pix = ang2pix(nside, THETA, PHI, nest=nest)
        grid_map = m[grid_pix]

        # plot
        if return_only_data:  # exit here when dumping the data
            return [longitude, latitude, grid_map]
        if projection_type != "3d":  # test for 3d plot
            vmin_, vmax_ = (None, None) if norm is not None else (vmin, vmax)
            ret = plt.pcolormesh(
                longitude,
                latitude,
                grid_map,
                vmin=vmin_,
                vmax=vmax_,
                rasterized=True,
                cmap=cmap,
                norm=norm,
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
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                rasterized=True,
                **kwargs
            )
    # graticule
    if rot_graticule or graticule_coord is not None:
        graticule_labels=False

    if rot_graticule or graticule_coord is None:
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

    if not graticule or not graticule_labels:
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
    if projection_type == "cart":
        ax.set_aspect(1)

    if cbar:
        # Override automatic tick generation with tick variable
        if ticks == [None, None]:
            ticks = [vmin, vmax]
            if vmin<0 and vmax>0: ticks.append(0)

        extend = "neither"
        if vmin > np.min(m):
            extend = "min"
        if vmax < np.max(m):
            extend = "max"
        if vmin > np.min(m) and vmax < np.max(m):
            extend = "both"

        cb = fig.colorbar(
            ret,
            orientation=cb_orientation,
            shrink=plot_properties["cbar_shrink"],
            pad=plot_properties["cbar_pad"],
            extend=extend,
        )

        # Hide all tickslabels not in tick variable. Do not delete tick-markers
        cbar_ticks = list(set(cb.get_ticks()) | set(ticks))
        labels = [tick if tick in ticks else "" for tick in cbar_ticks]
        args = np.argsort(cbar_ticks)
        cbar_ticks = list(np.array(cbar_ticks)[args])
        labels = list(np.array(labels)[args])
        
        cb.set_ticks(cbar_ticks, labels)
        cb.set_ticklabels(labels)
        
        if cb_orientation == "horizontal":
            cb.ax.xaxis.set_label_text(unit, fontsize=fontsize_defaults["cbar_label"])
            cb.ax.tick_params(axis="x", labelsize=fontsize_defaults["cbar_tick_label"], direction=plot_properties["cbar_tick_direction"], )
            cb.ax.xaxis.labelpad = plot_properties["cbar_label_pad"]
        if cb_orientation == "vertical":
            # Weird fix to keep labels from dissapearing
            labels = cb.ax.get_yticklabels() if norm is not None else labels
            cb.ax.set_yticklabels(labels, rotation=90, va="center",)
            cb.ax.yaxis.set_label_text(unit, fontsize=fontsize_defaults["cbar_label"], rotation=90)
            cb.ax.tick_params(axis="y", labelsize=fontsize_defaults["cbar_tick_label"], direction=plot_properties["cbar_tick_direction"], )
            cb.ax.yaxis.labelpad = plot_properties["cbar_label_pad"]
            
        # workaround for issue with viewers, see colorbar docstring
        cb.solids.set_edgecolor("face")
    ax.set_xlabel(xlabel, fontsize=fontsize_defaults["xlabel"])
    ax.set_ylabel(ylabel, fontsize=fontsize_defaults["ylabel"])

    # Separate graticule coordinate rotation
    if rot_graticule or graticule_coord is not None:
        if coord is None: coord="G"
        rotated_grid_lines, where_zero = CreateRotatedGraticule(
            rot,
            coordtransform=coord+graticule_coord,
            t_step=rot_graticule_properties["t_step"],
            p_step=rot_graticule_properties["p_step"],
        )
        for i, g_line in enumerate(rotated_grid_lines):
            if i in where_zero:
                linewidth = rot_graticule_properties["g_linewidth"] * 2.5
            else:
                linewidth = rot_graticule_properties["g_linewidth"]
            plt.plot(
                *g_line,
                linewidth=linewidth,
                linestyle=rot_graticule_properties["g_linestyle"],
                color=rot_graticule_properties["g_color"],
                alpha=rot_graticule_properties["g_alpha"]
            )

    #### Right label ####
    plt.text(
        0.975,
        0.925,
        rlabel,
        ha="right",
        va="center",
        fontsize=fontsize_defaults["cbar_label"],
        transform=ax.transAxes,
    )
    #### Left label ####
    plt.text(
        0.025,
        0.925,
        llabel,
        ha="left",
        va="center",
        fontsize=fontsize_defaults["cbar_label"],
        transform=ax.transAxes,
    )

    plt.draw()
    return ret


def newprojplot(theta, phi, fmt=None, **kwargs):
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

    Notes
    -----
    Other keywords are passed to :func:`matplotlib.Axes.plot`.
    """
    import matplotlib.pyplot as plt

    ax = plt.gca()
    flip = getattr(ax, "healpy_flip", "astro")

    longitude, latitude = lonlat(theta, phi)
    if flip == "astro":
        longitude = longitude * -1
    if fmt is None:
        ret = plt.plot(longitude, latitude, **kwargs)
    else:
        ret = plt.plot(longitude, latitude, fmt, **kwargs)
    return ret


def CreateRotatedGraticule(rot, t_step=30, p_step=30, coordtransform=None):
    if rot is None: rot=(0,0)
    # Transform graticule coordinate system
    #zyz = np.array(coordsys2euler_zyz(coordtransform))
    #rot = vec2dir(zyz*(180/np.pi),lonlat=True)
    coordtransform = "CG"
    if coordtransform == "GE":
        rot=(-90,0,-60)
    elif coordtransform == "EG":
        rot=(95,-60,0)
    elif coordtransform == "GC":
        rot=(-135,-45,-45)
    elif coordtransform == "CG":
        rot=(120,-60,10)

    phi = rot[0]
    try:
        theta = rot[1]
    except:
        theta = 0
    try:
        psi = rot[2]
    except:
        psi = 0

    pointDensity = 100
    conventionThetaOffset = np.pi / 2
    phiSpacing = np.arange(-180, 180 + p_step, p_step)
    thetaSpacing = np.arange(-90, 90 + t_step, t_step)
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
        r = Rotator(rot=(0, theta, psi), inv=True)
        gline_theta_rot, gline_phi_fixed_rot = r(gline_theta, gline_phi_fixed)
        gline_theta = gline_theta - conventionThetaOffset
        gline_theta_rot = gline_theta_rot - conventionThetaOffset

        rotated_grid_lines.append([gline_phi_fixed_rot, gline_theta_rot])

    for phiSpace in phiSpacing:
        gline_phi = np.deg2rad(np.zeros(pointDensity) + phiSpace)
        r = Rotator(rot=(phi, 0,), inv=True)
        gline_theta_fixed_rot, gline_phi_rot = r(gline_theta_fixed, gline_phi)
        r = Rotator(rot=(0, theta, psi), inv=True)
        gline_theta_fixed_rot, gline_phi_rot = r(gline_theta_fixed_rot, gline_phi_rot)
        gline_theta_fixed_rot = gline_theta_fixed_rot - conventionThetaOffset

        rotated_grid_lines.append([gline_phi_rot, gline_theta_fixed_rot])

    for g_lines in rotated_grid_lines:
        mask = np.where((np.abs(np.diff(g_lines[0]))) > np.deg2rad(45))
        g_lines[0] = np.ma.array(g_lines[0])
        g_lines[0][mask] = np.ma.masked




    return rotated_grid_lines, where_zero

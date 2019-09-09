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

from . import projaxes as PA
from . import rotator as R
import numpy as np
import matplotlib
from ._healpy_pixel_lib import UNSEEN
from . import pixelfunc

pi = np.pi
dtor = pi / 180.0


def mollzoom(
    map=None,
    fig=None,
    rot=None,
    coord=None,
    unit="",
    xsize=800,
    title="Mollweide view",
    nest=False,
    min=None,
    max=None,
    flip="astro",
    remove_dip=False,
    remove_mono=False,
    gal_cut=0,
    format="%g",
    cmap=None,
    norm=None,
    hold=False,
    margins=None,
    sub=None,
):
    """Interactive mollweide plot with zoomed gnomview.
    
    Parameters:
    -----------
    map : float, array-like shape (Npix,)
      An array containing the map, 
      supports masked maps, see the `ma` function.
      if None, use map with inf value (white map), useful for
      overplotting
    fig : a figure number. 
      Default: create a new figure
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
    title : str, optional
      The title of the plot. Default: 'Mollweide view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, optional
      The minimum range value
    max : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    format : str, optional
      The format of the scale label. Default: '%g'
    """
    import pylab

    # Ensure that the nside is valid
    nside = pixelfunc.get_nside(map)
    pixelfunc.check_nside(nside, nest=nest)

    # create the figure (if interactive, it will open the window now)
    f = pylab.figure(fig, figsize=(10.5, 5.4))
    extent = (0.02, 0.25, 0.56, 0.72)
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12) + np.inf
        map = pixelfunc.ma_to_array(map)
        ax = PA.HpxMollweideAxes(
            f, extent, coord=coord, rot=rot, format=format, flipconv=flip
        )
        f.add_axes(ax)
        if remove_dip:
            map = pixelfunc.remove_dipole(
                map, gal_cut=gal_cut, nest=nest, copy=True, verbose=True
            )
        elif remove_mono:
            map = pixelfunc.remove_monopole(
                map, gal_cut=gal_cut, nest=nest, copy=True, verbose=True
            )
        ax.projmap(
            map,
            nest=nest,
            xsize=xsize,
            coord=coord,
            vmin=min,
            vmax=max,
            cmap=cmap,
            norm=norm,
        )
        im = ax.get_images()[0]
        b = im.norm.inverse(np.linspace(0, 1, im.cmap.N + 1))
        v = np.linspace(im.norm.vmin, im.norm.vmax, im.cmap.N)
        if matplotlib.__version__ >= "0.91.0":
            cb = f.colorbar(
                ax.get_images()[0],
                ax=ax,
                orientation="horizontal",
                shrink=0.5,
                aspect=25,
                ticks=PA.BoundaryLocator(),
                pad=0.05,
                fraction=0.1,
                boundaries=b,
                values=v,
            )
        else:
            # for older matplotlib versions, no ax kwarg
            cb = f.colorbar(
                ax.get_images()[0],
                orientation="horizontal",
                shrink=0.5,
                aspect=25,
                ticks=PA.BoundaryLocator(),
                pad=0.05,
                fraction=0.1,
                boundaries=b,
                values=v,
            )
        ax.set_title(title)
        ax.text(
            0.86,
            0.05,
            ax.proj.coordsysstr,
            fontsize=14,
            fontweight="bold",
            transform=ax.transAxes,
        )
        cb.ax.text(
            1.05,
            0.30,
            unit,
            fontsize=14,
            fontweight="bold",
            transform=cb.ax.transAxes,
            ha="left",
            va="center",
        )
        f.sca(ax)

        ## Gnomonic axes
        # extent = (0.02,0.25,0.56,0.72)
        g_xsize = 600
        g_reso = 1.0
        extent = (0.60, 0.04, 0.38, 0.94)
        g_ax = PA.HpxGnomonicAxes(
            f, extent, coord=coord, rot=rot, format=format, flipconv=flip
        )
        f.add_axes(g_ax)
        if remove_dip:
            map = pixelfunc.remove_dipole(map, gal_cut=gal_cut, nest=nest, copy=True)
        elif remove_mono:
            map = pixelfunc.remove_monopole(map, gal_cut=gal_cut, nest=nest, copy=True)
        g_ax.projmap(
            map,
            nest=nest,
            coord=coord,
            vmin=min,
            vmax=max,
            xsize=g_xsize,
            ysize=g_xsize,
            reso=g_reso,
            cmap=cmap,
            norm=norm,
        )
        im = g_ax.get_images()[0]
        b = im.norm.inverse(np.linspace(0, 1, im.cmap.N + 1))
        v = np.linspace(im.norm.vmin, im.norm.vmax, im.cmap.N)
        if matplotlib.__version__ >= "0.91.0":
            cb = f.colorbar(
                g_ax.get_images()[0],
                ax=g_ax,
                orientation="horizontal",
                shrink=0.5,
                aspect=25,
                ticks=PA.BoundaryLocator(),
                pad=0.08,
                fraction=0.1,
                boundaries=b,
                values=v,
            )
        else:
            cb = f.colorbar(
                g_ax.get_images()[0],
                orientation="horizontal",
                shrink=0.5,
                aspect=25,
                ticks=PA.BoundaryLocator(),
                pad=0.08,
                fraction=0.1,
                boundaries=b,
                values=v,
            )
        g_ax.set_title(title)
        g_ax.text(
            -0.07,
            0.02,
            "%g '/pix,   %dx%d pix"
            % (
                g_ax.proj.arrayinfo["reso"],
                g_ax.proj.arrayinfo["xsize"],
                g_ax.proj.arrayinfo["ysize"],
            ),
            fontsize=12,
            verticalalignment="bottom",
            transform=g_ax.transAxes,
            rotation=90,
        )
        g_ax.text(
            -0.07,
            0.8,
            g_ax.proj.coordsysstr,
            fontsize=14,
            fontweight="bold",
            rotation=90,
            transform=g_ax.transAxes,
        )
        lon, lat = np.around(g_ax.proj.get_center(lonlat=True), g_ax._coordprec)
        g_ax.text(
            0.5,
            -0.03,
            "on (%g,%g)" % (lon, lat),
            verticalalignment="center",
            horizontalalignment="center",
            transform=g_ax.transAxes,
        )
        cb.ax.text(
            1.05,
            0.30,
            unit,
            fontsize=14,
            fontweight="bold",
            transform=cb.ax.transAxes,
            ha="left",
            va="center",
        )
        # Add graticule info axes
        grat_ax = pylab.axes([0.25, 0.02, 0.22, 0.25])
        grat_ax.axis("off")
        # Add help text
        help_ax = pylab.axes([0.02, 0.02, 0.22, 0.25])
        help_ax.axis("off")
        t = help_ax.transAxes
        help_ax.text(0.1, 0.8, "r/t .... zoom out/in", transform=t, va="baseline")
        help_ax.text(0.1, 0.65, "p/v .... print coord/val", transform=t, va="baseline")
        help_ax.text(0.1, 0.5, "c ...... go to center", transform=t, va="baseline")
        help_ax.text(0.1, 0.35, "f ...... next color scale", transform=t, va="baseline")
        help_ax.text(
            0.1, 0.2, "k ...... save current scale", transform=t, va="baseline"
        )
        help_ax.text(0.1, 0.05, "g ...... toggle graticule", transform=t, va="baseline")
        f.sca(g_ax)
        # Set up the zoom capability
        zt = ZoomTool(map, fig=f.number, nest=nest, cmap=cmap, norm=norm, coord=coord)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()


def set_g_clim(vmin, vmax):
    """Set min/max value of the gnomview part of a mollzoom.
    """
    import pylab

    f = pylab.gcf()
    if not hasattr(f, "zoomtool"):
        raise TypeError("The current figure has no zoomtool")
    f.zoomtool.save_min = vmin
    f.zoomtool.save_max = vmax
    f.zoomtool._range_status = 2
    f.zoomtool.draw_gnom()


class ZoomTool(object):
    """A class providing zoom capability to a figure containing a Mollweide
    and a Gnomonic axis.
    """

    def __init__(self, m, fig=None, nest=False, cmap=None, norm=None, coord=None):
        """m: the map to be zoomed (already plotted in Mollweide view)
        fig: the figure to instrument (None->gcf())
        """
        import pylab

        self.reso_list = [
            0.05,
            0.1,
            0.2,
            0.3,
            0.5,
            0.75,
            1.0,
            1.5,
            3.0,
            5.0,
            10.0,
            15.0,
            30.0,
            45.0,
            60.0,
        ]
        self._map = m
        self._nest = nest
        self._cmap = cmap
        self._norm = norm
        self._coord = coord
        self._range_status = 0  # 0:normal, 1:global map min,max, 2: saved
        self.save_min = self.save_max = None
        self._graton = False
        # find min, max of map
        if isinstance(m, dict):
            if len(m) == 0:
                self._mapmin, self._mapmax = -1.0, 1.0
            else:
                self._mapmin, self._mapmax = min(m.values()), max(m.values())
        else:
            mgood = m[m != UNSEEN]
            if mgood.size == 0:
                self._mapmin, self._mapmax = -1.0, 1.0
            else:
                self._mapmin, self._mapmax = mgood.min(), mgood.max()
            del mgood
        if fig is None:
            f = pylab.gcf()
        else:
            f = pylab.figure(fig)
        self.f = f
        f.zoomtool = self
        (
            self._moll_ax,
            self._moll_cb_ax,
            self._gnom_ax,
            self._gnom_cb_ax,
        ) = f.get_axes()[:4]
        self._grat_ax = f.get_axes()[4]
        self._text_reso, self._text_coord, self._text_loc = self._gnom_ax.texts
        self._xsize = self._gnom_ax.proj.arrayinfo["xsize"]
        self._ysize = self._gnom_ax.proj.arrayinfo["ysize"]
        try:
            self._reso_idx = self.reso_list.index(self._gnom_ax.proj._arrayinfo["reso"])
        except ValueError as e:
            raise ValueError("Resolution not in %s" % self.reso_list)
        self.zoomcenter, = self._moll_ax.plot([0], [0], "ok", mew=1, ms=15, alpha=0.1)
        self.zoomcenter2, = self._moll_ax.plot([0], [0], "xr", ms=15, alpha=0.5, mew=3)
        self._text_range = self._gnom_ax.text(
            -0.4,
            -0.2,
            "scale mode: loc",
            transform=self._gnom_ax.transAxes,
            va="baseline",
            ha="left",
        )
        self.draw_gnom(0, 0)
        self._connected = False
        self.connect_callbacks()

    def _zoom_on_click(self, ev):
        import pylab

        try:
            ax = ev.inaxes
            lon, lat = ax.get_lonlat(ev.xdata, ev.ydata)
            if np.isnan(lon) or np.isnan(lat):
                raise ValueError("invalid position")
            val = ax.get_value(ev.xdata, ev.ydata)
            self.lastval = val
            self._move_zoom_center(lon, lat)
            self.draw_gnom(lon, lat)
        except Exception as s:
            self._move_zoom_center(0, 0, False)
            pylab.draw_if_interactive()
            # print s
        return

    def _reso_on_key(self, ev):
        if ev.key == "r":
            self._decrease_reso()
        elif ev.key == "t":
            self._increase_reso()
        elif ev.key == "p":
            print("lon,lat = %.17g,%.17g" % (self.lon, self.lat))
        elif ev.key == "c":
            self._move_zoom_center(0, 0)
            self.draw_gnom(0, 0)
        elif ev.key == "v":
            print("val = %.17g" % (self.lastval))
        elif ev.key == "f":
            self._range_status += 1
            self._range_status %= 3
            self.draw_gnom()
        elif ev.key == "k":
            self.save_min = self._gnom_ax.images[0].norm.vmin
            self.save_max = self._gnom_ax.images[0].norm.vmax
        elif ev.key == "g":
            if hasattr(self, "_graton") and self._graton == True:
                self._gnom_ax.delgraticules()
                self._moll_ax.delgraticules()
                self._graton = False
            else:
                (self._g_dpar, self._g_dmer) = self._gnom_ax.graticule(
                    local=False, verbose=False
                )
                (self._m_dpar, self._m_dmer) = self._moll_ax.graticule(verbose=False)
                self._graton = True
            self.draw_gnom()

    def _update_grat_info(self):
        self._grat_ax.cla()
        self._grat_ax.axis("off")
        if self._graton:
            a = self._grat_ax
            t = a.transAxes
            a.text(0.1, 0.8, "moll. grat.:", transform=t, weight="bold")
            vdeg = np.floor(np.around(self._m_dpar / dtor, 10))
            varcmin = (self._m_dpar / dtor - vdeg) * 60.0
            a.text(0.1, 0.65, "   -par: %d d %.2f '" % (vdeg, varcmin), transform=t)
            vdeg = np.floor(np.around(self._m_dmer / dtor, 10))
            varcmin = (self._m_dmer / dtor - vdeg) * 60.0
            a.text(0.1, 0.5, "   -mer: %d d %.2f '" % (vdeg, varcmin), transform=t)
            a.text(0.1, 0.35, "gnom. grat.:", transform=t, weight="bold")
            vdeg = np.floor(np.around(self._g_dpar / dtor, 10))
            varcmin = (self._g_dpar / dtor - vdeg) * 60.0
            a.text(0.1, 0.2, "   -par: %d d %.2f '" % (vdeg, varcmin), transform=t)
            vdeg = np.floor(np.around(self._g_dmer / dtor, 10))
            varcmin = (self._g_dmer / dtor - vdeg) * 60.0
            a.text(0.1, 0.05, "   -mer: %d d %.2f '" % (vdeg, varcmin), transform=t)

    def _increase_reso(self):
        if self._reso_idx > 0:
            self._reso_idx -= 1
            self.draw_gnom(self.lon, self.lat)

    def _decrease_reso(self):
        if self._reso_idx < len(self.reso_list) - 1:
            self._reso_idx += 1
            self.draw_gnom(self.lon, self.lat)

    def get_reso(self):
        return self.reso_list[self._reso_idx]

    def connect_callbacks(self):
        if not self._connected:
            self._callbacks_id = []
            cid = self.f.canvas.mpl_connect("button_press_event", self._zoom_on_click)
            self._callbacks_id.append(cid)
            cid = self.f.canvas.mpl_connect("key_press_event", self._reso_on_key)
            self._callbacks_id.append(cid)
            self._connected = True

    def disconnect_callbacks(self):
        if self._connected:
            for cid in self._callbacks_id:
                self.figure.canvas.mpl_disconnect(cid)

    def _move_zoom_center(self, lon, lat, visible=True):
        # Move the zoom center marker.
        if self.zoomcenter:
            x, y = self._moll_ax.proj.ang2xy(lon, lat, lonlat=True)
            self.zoomcenter.set_xdata([x])
            self.zoomcenter.set_ydata([y])
            self.zoomcenter.set_visible(visible)
        if self.zoomcenter2:
            x, y = self._moll_ax.proj.ang2xy(lon, lat, lonlat=True)
            self.zoomcenter2.set_xdata([x])
            self.zoomcenter2.set_ydata([y])
            self.zoomcenter2.set_visible(visible)

    def draw_gnom(self, lon=None, lat=None):
        import pylab

        wasinteractive = pylab.isinteractive()
        pylab.ioff()
        try:
            # modify rot of the gnom_ax
            if lon is None:
                lon = self._lon
            else:
                self._lon = lon
            if lat is None:
                lat = self._lat
            else:
                self._lat = lat
            self._gnom_ax.proj.rotator._rots.pop()
            self._gnom_ax.proj.rotator._rots.append(
                R.normalise_rot((lon, lat), deg=True)
            )
            self._gnom_ax.proj.rotator._update_matrix()
            if self._range_status == 0:
                vmin = vmax = None
            elif self._range_status == 1:
                vmin, vmax = self._mapmin, self._mapmax
            elif self._range_status == 2:
                vmin, vmax = self.save_min, self.save_max
            self._gnom_ax.images.pop()
            self._gnom_ax.projmap(
                self._map,
                nest=self._nest,
                coord=self._coord,
                vmin=vmin,
                vmax=vmax,
                xsize=self._xsize,
                ysize=self._ysize,
                reso=self.get_reso(),
                cmap=self._cmap,
                norm=self._norm,
            )
            if hasattr(self._gnom_ax, "_scatter_data"):
                l = [x for x in self._gnom_ax._scatter_data]
                # print l
                for sd in l:
                    s, input_data = sd
                    # print input_data
                    self._gnom_ax.collections.remove(s)
                    self._gnom_ax._scatter_data.remove(sd)
                    theta, phi, args, kwds = input_data
                    self._gnom_ax.projscatter(theta, phi=phi, *args, **kwds)
                del l
            if self._graton:
                self._gnom_ax.delgraticules()
                (self._g_dpar, self._g_dmer) = self._gnom_ax.graticule(
                    local=False, verbose=False
                )
            self._gnom_cb_ax.cla()
            im = self._gnom_ax.images[0]
            if matplotlib.__version__ >= "0.91.0":
                cb = self.f.colorbar(
                    im,
                    ax=self._gnom_ax,
                    cax=self._gnom_cb_ax,
                    orientation="horizontal",
                    ticks=PA.BoundaryLocator(),
                )
            else:
                cb = self.f.colorbar(
                    im,
                    cax=self._gnom_cb_ax,
                    orientation="horizontal",
                    ticks=PA.BoundaryLocator(),
                )
            lon, lat = np.around(
                self._gnom_ax.proj.get_center(lonlat=True), self._gnom_ax._coordprec
            )
            self._text_loc.set_text("on (%g,%g)" % (lon, lat))
            reso = self._gnom_ax.proj.arrayinfo["reso"]
            xsize = self._gnom_ax.proj.arrayinfo["xsize"]
            ysize = self._gnom_ax.proj.arrayinfo["ysize"]
            self._text_reso.set_text("%g '/pix,   %dx%d pix" % (reso, xsize, ysize))
            mode = ["loc", "map", "sav"][self._range_status]
            self._text_range.set_text("scale mode: %s" % mode)
            self.lon, self.lat = lon, lat
            self._update_grat_info()
        except Exception as e:
            pass  # print e
        finally:
            if wasinteractive:
                pylab.ion()
                pylab.draw()
                pylab.show()

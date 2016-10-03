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

"""
=====================================================
visufunc.py : Healpix visualization functions
=====================================================

This module provides functions to display Healpix map.

Map projections
---------------

- :func:`mollview` displays a map using Mollweide projection (full sky)
- :func:`gnomview` displays a map using Gnomonic projection (local map)
- :func:`cartview` displays a map using Cartesian projection
- :func:`orthview` displays a map using Orthographic projection (full or half sky)
- :func:`azeqview` displays a map using Azimuthal equidistant projection

Graticules
----------

- :func:`graticule` adds a graticule to the current map
- :func:`delgraticules` deletes all graticules of a map

Tracing lines or points
-----------------------

- :func:`projplot` plots data points on the current map
- :func:`projscatter` displays scatter points
- :func:`projtext` display a text on the current map
"""

__all__ = ['mollview', 'gnomview', 'cartview', 'orthview', 'azeqview',
           'graticule', 'delgraticules',
           'projplot', 'projscatter', 'projtext']

from . import projaxes as PA
import numpy as np
import matplotlib
from . import pixelfunc

pi = np.pi
dtor = pi/180.

def mollview(map=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,title='Mollweide view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None,
             return_projected_map=False):
    """Plot an healpix map (given as an array) in Mollweide projection.
    
    Parameters
    ----------
    map : float, array-like or None
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
    format2 : str, optional
      Format of the pixel value under mouse. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by a MollweideAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    gnomview, cartview, orthview, azeqview
    """
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub//100, (sub%100)//10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        ax=PA.HpxMollweideAxes(f,extent,coord=coord,rot=rot,
                               format=format2,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        img = ax.projmap(map,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                   cmap=cmap,norm=norm)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0,1,im.cmap.N+1))
            v = np.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(im,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(im,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def gnomview(map=None,fig=None,rot=None,coord=None,unit='',
             xsize=200,ysize=None,reso=1.5,
             title='Gnomonic view',nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,flip='astro',
             format='%.3g',cbar=True,
             cmap=None, norm=None,
             hold=False,sub=None,margins=None,notext=False,
             return_projected_map=False):
    """Plot an healpix map (given as an array) in Gnomonic projection.

    Parameters
    ----------
    map : array-like
      The map to project, supports masked maps, see the `ma` function.
      If None, use a blank map, useful for
      overplotting.
    fig : None or int, optional
      A figure number. Default: None= create a new figure
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
      The size of the image. Default: 200
    ysize : None or int, optional
      The size of the image. Default: None= xsize
    reso : float, optional
      Resolution (in arcmin). Default: 1.5 arcmin
    title : str, optional
      The title of the plot. Default: 'Gnomonic view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, scalar, optional
      The minimum range value
    max : float, scalar, optional
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
    hold : bool, optional
      If True, replace the current Axes by a MollweideAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    notext: bool, optional
      If True: do not add resolution info text. Default=False
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    mollview, cartview, orthview, azeqview
    """
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(5.8,6.4))
        if not margins:
                margins = (0.075,0.05,0.075,0.05)
        extent = (0.0,0.0,1.0,1.0)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(pylab.gca().get_position()).ravel()
        if not margins:
            margins = (0.0,0.0,0.0,0.0)
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(pylab.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub//100, (sub%100)//10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
    extent = (extent[0]+margins[0],
              extent[1]+margins[1],
              extent[2]-margins[2]-margins[0],
              extent[3]-margins[3]-margins[1])
    #f=pylab.figure(fig,figsize=(5.5,6))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        ax=PA.HpxGnomonicAxes(f,extent,coord=coord,rot=rot,
                              format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)
        img = ax.projmap(map,nest=nest,coord=coord,vmin=min,vmax=max,
                   xsize=xsize,ysize=ysize,reso=reso,cmap=cmap,norm=norm)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0,1,im.cmap.N+1))
            v = np.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(im,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                cb=f.colorbar(im,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                  
        ax.set_title(title)
        if not notext:
            ax.text(-0.07,0.02,
                     "%g '/pix,   %dx%d pix"%(ax.proj.arrayinfo['reso'],
                                              ax.proj.arrayinfo['xsize'],
                                              ax.proj.arrayinfo['ysize']),
                     fontsize=12,verticalalignment='bottom',
                     transform=ax.transAxes,rotation=90)
            ax.text(-0.07,0.6,ax.proj.coordsysstr,fontsize=14,
                     fontweight='bold',rotation=90,transform=ax.transAxes)
            lon,lat = np.around(ax.proj.get_center(lonlat=True),ax._coordprec)
            ax.text(0.5,-0.03,'(%g,%g)'%(lon,lat),
                    verticalalignment='center', horizontalalignment='center',
                    transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def cartview(map=None,fig=None,rot=None,zat=None,coord=None,unit='',
             xsize=800,ysize=None,lonra=None,latra=None,
             title='Cartesian view',nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,flip='astro',
             format='%.3g',cbar=True,
             cmap=None, norm=None,aspect=None,
             hold=False,sub=None,margins=None,notext=False,
             return_projected_map=False):
    """Plot an healpix map (given as an array) in Cartesian projection.

    Parameters
    ----------
    map : float, array-like or None
      An array containing the map, 
      supports masked maps, see the `ma` function.
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
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    lonra : sequence, optional
      Range in longitude. Default: [-180,180]
    latra : sequence, optional
      Range in latitude. Default: [-90,90]
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
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}, optional
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by a CartesianAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    mollview, gnomview, orthview, azeqview
    """
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        if not margins:
                margins = (0.075,0.05,0.075,0.05)
        extent = (0.0,0.0,1.0,1.0)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(pylab.gca().get_position()).ravel()
        if not margins:
            margins = (0.0,0.0,0.0,0.0)
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(pylab.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub//100, (sub%100)//10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
    extent = (extent[0]+margins[0],
              extent[1]+margins[1],
              extent[2]-margins[2]-margins[0],
              extent[3]-margins[3]-margins[1])

    #f=pylab.figure(fig,figsize=(5.5,6))
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        if zat and rot:
            raise ValueError('Only give rot or zat, not both')
        if zat:
            rot = np.array(zat,dtype=np.float64)
            rot.resize(3)
            rot[1] -= 90
        ax=PA.HpxCartesianAxes(f,extent,coord=coord,rot=rot,
                               format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)
        img = ax.projmap(map,nest=nest,coord=coord,vmin=min,vmax=max,
                   xsize=xsize,ysize=ysize,lonra=lonra,latra=latra,
                   cmap=cmap,norm=norm,aspect=aspect)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0,1,im.cmap.N+1))
            v = np.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(im,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                cb=f.colorbar(im,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                 
        ax.set_title(title)
        if not notext:
            ax.text(-0.07,0.6,ax.proj.coordsysstr,fontsize=14,
                     fontweight='bold',rotation=90,transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
    if return_projected_map:
        return img

def orthview(map=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,half_sky=False,
             title='Orthographic view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None,
             return_projected_map=False):
    """Plot an healpix map (given as an array) in Orthographic projection.
    
    Parameters
    ----------
    map : float, array-like or None
      An array containing the map.
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
    half_sky : bool, optional
      Plot only one side of the sphere. Default: False
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    title : str, optional
      The title of the plot. Default: 'Orthographic view'
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
    format2 : str, optional
      Format of the pixel value under mouse. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by an OrthographicAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array
    
    See Also
    --------
    mollview, gnomview, cartview, azeqview
    """
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub//100, (sub%100)//10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0],
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))
    
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        ax=PA.HpxOrthographicAxes(f,extent,coord=coord,rot=rot,
                                  format=format2,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        img = ax.projmap(map,nest=nest,xsize=xsize,half_sky=half_sky,
                         coord=coord,vmin=min,vmax=max,
                         cmap=cmap,norm=norm)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0,1,im.cmap.N+1))
            v = np.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(im,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(im,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                  
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def azeqview(map=None,fig=None,rot=None,zat=None,coord=None,unit='',
             xsize=800,ysize=None,reso=1.5,lamb=False,half_sky=False,
             title=None,nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,flip='astro',
             format='%.3g',cbar=True,
             cmap=None, norm=None,aspect=None,
             hold=False,sub=None,margins=None,notext=False,
             return_projected_map=False):
    """Plot a healpix map (given as an array) in Azimuthal equidistant projection
    or Lambert azimuthal equal-area projection.

    Parameters
    ----------
    map : float, array-like or None
      An array containing the map,
      supports masked maps, see the `ma` function.
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
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    ysize : None or int, optional
      The size of the image. Default: None= xsize
    reso : float, optional
      Resolution (in arcmin). Default: 1.5 arcmin
    lamb : bool, optional
      If True, plot Lambert azimuthal equal area instead of azimuthal
      equidistant. Default: False (az equidistant)
    half_sky : bool, optional
      Plot only one side of the sphere. Default: False
    title : str, optional
      The title of the plot. Default: 'Azimuthal equidistant view'
      or 'Lambert azimuthal equal-area view' (if lamb is True)
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
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by an Equidistant AzimuthalAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    mollview, gnomview, cartview, orthview
    """
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub//100, (sub%100)//10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0],
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        ax=PA.HpxAzimuthalAxes(f,extent,coord=coord,rot=rot,
                               format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        img = ax.projmap(map,nest=nest,xsize=xsize,ysize=ysize,reso=reso,lamb=lamb,
                         half_sky=half_sky,coord=coord,vmin=min,vmax=max,
                         cmap=cmap,norm=norm)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(np.linspace(0,1,im.cmap.N+1))
            v = np.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(im,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(im,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)
        if title is None:
            if lamb:
                title = 'Lambert azimuthal equal-area view'
            else:
                title = 'Azimuthal equidistant view'
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img
        

def graticule(dpar=None,dmer=None,coord=None,local=None,**kwds):
    """Draw a graticule on the current Axes.

    Parameters
    ----------
    dpar, dmer : float, scalars
      Interval in degrees between meridians and between parallels
    coord : {'E', 'G', 'C'}
      The coordinate system of the graticule (make rotation if needed,
      using coordinate system of the map if it is defined).
    local : bool
      If True, draw a local graticule (no rotation is performed, useful for
      a gnomonic view, for example)

    Notes
    -----
    Other keyword parameters will be transmitted to the projplot function.

    See Also
    --------
    delgraticules
    """
    import pylab
    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if len(f.get_axes()) == 0:
            ax=PA.HpxMollweideAxes(f,(0.02,0.05,0.96,0.9),coord=coord)
            f.add_axes(ax)
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.graticule(dpar=dpar,dmer=dmer,coord=coord,
                             local=local,**kwds)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
    
def delgraticules():
    """Delete all graticules previously created on the Axes.

    See Also
    --------
    graticule
    """
    import pylab
    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.delgraticules()
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()

def projplot(*args,**kwds):
    import pylab
    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret = None
    try:
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projplot(*args,**kwds)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    return ret
projplot.__doc__ = PA.SphericalProjAxes.projplot.__doc__
    
def projscatter(*args,**kwds):
    import pylab
    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret=None
    try:
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projscatter(*args,**kwds)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    return ret
projscatter.__doc__ = PA.SphericalProjAxes.projscatter.__doc__

def projtext(*args,**kwds):
    import pylab
    f = pylab.gcf()
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret = None
    try:
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projtext(*args,**kwds)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    return ret
projtext.__doc__ = PA.SphericalProjAxes.projtext.__doc__

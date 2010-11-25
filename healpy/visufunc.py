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
import projaxes as PA
import pylab
import numpy as npy
import matplotlib
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import pixelfunc

pi = npy.pi
dtor = pi/180.

def mollview(map=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,title='Mollweide view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None):
    """Plot an healpix map (given as an array) in Mollweide projection.
    
    Input:
      - map : an ndarray containing the map
              if None, use map with inf value (white map), useful for
              overplotting
    Parameters:
      - fig: a figure number. Default: create a new figure
      - rot: rotation, either 1,2 or 3 angles describing the rotation
             Default: None
      - coord: either one of 'G', 'E' or 'C' to describe the coordinate
               system of the map, or a sequence of 2 of these to make
               rotation from the first to the second coordinate system.
               Default: None
      - unit: a text describing the unit. Default: ''
      - xsize: the size of the image. Default: 800
      - title: the title of the plot. Default: 'Mollweide view'
      - nest: if True, ordering scheme is NEST. Default: False (RING)
      - min: the minimum range value
      - max: the maximum range value
      - flip: 'astro' (default, east towards left, west towards right) or 'geo'
      - remove_dip: if True, remove the dipole+monopole
      - remove_mono: if True, remove the monopole
      - gal_cut: galactic cut for the dipole/monopole fit
      - format: the format of the scale label. Default: '%g'
      - format2: format of the pixel value under mouse. Default: '%g'
      - cbar: display the colorbar. Default: True
      - notext: if True, no text is printed around the map
      - norm: color normalization, hist= histogram equalized color mapping, log=
              logarithmic color mapping, default: None (linear color mapping)
      - hold: if True, replace the current Axes by a MollweideAxes.
              use this if you want to have multiple maps on the same
              figure. Default: False
      - sub: use a part of the current figure (same syntax as subplot).
             Default: None
      - margins: either None, or a sequence (left,bottom,right,top)
                 giving the margins on left,bottom,right and top
                 of the axes. Values are relative to figure (0-1).
                 Default: None
    """
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = npy.zeros(12)+npy.inf
            cbar=False
        if not (hold or sub):
            f=pylab.figure(fig,figsize=(8.5,5.4))
            extent = (0.02,0.05,0.96,0.9)
        elif hold:
            f=pylab.gcf()
            left,bottom,right,top = npy.array(f.gca().get_position()).ravel()
            extent = (left,bottom,right-left,top-bottom)
            f.delaxes(f.gca())
        else: # using subplot syntax
            f=pylab.gcf()
            if hasattr(sub,'__len__'):
                nrows, ncols, idx = sub
            else:
                nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
            if idx < 1 or idx > ncols*nrows:
                raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                     ncols,
                                                                     idx))
            c,r = (idx-1)%ncols,(idx-1)/ncols
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
        ax.projmap(map,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                   cmap=cmap,norm=norm)
        if cbar:
            im = ax.get_images()[0]
            b = im.norm.inverse(npy.linspace(0,1,im.cmap.N+1))
            v = npy.linspace(im.norm.vmin,im.norm.vmax,im.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(ax.get_images()[0],ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(ax.get_images()[0],orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,0.10,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()


def gnomview(map=None,fig=None,rot=None,coord=None,unit='',
             xsize=200,ysize=None,reso=1.5,degree=False,
             title='Gnomonic view',nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,flip='astro',
             format='%g',cbar=True,
             cmap=None, norm=None,
             hold=False,sub=None,margins=None,notext=False):
    """Plot an healpix map (given as an array) in Gnomonic projection.

    Input:
      - map : an ndarray containing the map.
              if None, use map with inf value (white map), useful for
              overplotting
    Parameters:
      - fig: a figure number. Default: create a new figure
      - rot: rotation, either 1,2 or 3 angles describing the rotation
             Default: None
      - coord: either one of 'G', 'E' or 'C' to describe the coordinate
               system of the map, or a sequence of 2 of these to make
               rotation from the first to the second coordinate system.
               Default: None
      - unit: a text describing the unit. Default: ''
      - xsize: the size of the image. Default: 200
      - ysize: the size of the image. Default: xsize
      - reso: resolution in arcmin if degree is False. Default: 1.5 arcmin
      - degree: if True, reso is in degree. Default: False
      - title: the title of the plot. Default: 'Mollweide view'
      - nest: if True, ordering scheme is NEST. Default: False (RING)
      - min: the minimum range value
      - max: the maximum range value
      - flip: 'astro' (default, east towards left, west towards right) or 'geo'
      - remove_dip: if True, remove the dipole+monopole
      - remove_mono: if True, remove the monopole
      - gal_cut: galactic cut for the dipole/monopole fit
      - format: the format of the scale. Default: '%.3g'
      - hold: if True, replace the current Axes by a MollweideAxes.
              use this if you want to have multiple maps on the same
              figure. Default: False
      - sub: use a part of the current figure (same syntax as subplot).
             Default: None
      - margins: either None, or a sequence (left,bottom,right,top)
                 giving the margins on left,bottom,right and top
                 of the axes. Values are relative to figure (0-1).
                 Default: None
      - notext: True: do not add resolution info text
                Default=False
    """
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = npy.zeros(12)+npy.inf
            cbar=False
        if not (hold or sub):
            f=pylab.figure(fig,figsize=(5.5,6))
            if not margins:
                    margins = (0.075,0.05,0.075,0.05)
            extent = (0.0,0.0,1.0,1.0)
        elif hold:
            f=pylab.gcf()
            left,bottom,right,top = npy.array(pylab.gca().get_position()).ravel()
            if not margins:
                margins = (0.0,0.0,0.0,0.0)
            extent = (left,bottom,right-left,top-bottom)
            f.delaxes(pylab.gca())
        else: # using subplot syntax
            f=pylab.gcf()
            if hasattr(sub,'__len__'):
                nrows, ncols, idx = sub
            else:
                nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
            if idx < 1 or idx > ncols*nrows:
                raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                     ncols,
                                                                     idx))
            c,r = (idx-1)%ncols,(idx-1)/ncols
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
        ax=PA.HpxGnomonicAxes(f,extent,coord=coord,rot=rot,
                              format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)
        ax.projmap(map,nest=nest,coord=coord,vmin=min,vmax=max,
                   xsize=xsize,ysize=ysize,reso=reso,cmap=cmap,norm=norm)
        if cbar:
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(ax.get_images()[0],ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1)
            else:
                cb=f.colorbar(ax.get_images()[0],orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1)
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
            lon,lat = npy.around(ax.proj.get_center(lonlat=True),ax._coordprec)
            ax.text(0.5,-0.03,'on (%g,%g)'%(lon,lat),
                    verticalalignment='center', horizontalalignment='center',
                    transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()


def cartview(map=None,fig=None,rot=None,zat=None,coord=None,unit='',
             xsize=800,ysize=None,lonra=None,latra=None,
             title='Cartesian view',nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,flip='astro',
             format='%g',cbar=True,
             cmap=None, norm=None,aspect=None,
             hold=False,sub=None,margins=None,notext=False):
    """Plot an healpix map (given as an array) in Cartesian projection.

    Input:
      - map : an ndarray containing the map.
              if None, use map with inf value (white map), useful for
              overplotting
    Parameters:
      - fig: a figure number. Default: create a new figure
      - rot: rotation, either 1,2 or 3 angles describing the rotation
             Default: None
      - coord: either one of 'G', 'E' or 'C' to describe the coordinate
               system of the map, or a sequence of 2 of these to make
               rotation from the first to the second coordinate system.
               Default: None
      - unit: a text describing the unit. Default: ''
      - xsize: the size of the image. Default: 200
      - lonra: range in longitude. Default: [-180,180]
      - latra: range in latitude. Default: [-90,90]
      - title: the title of the plot. Default: 'Mollweide view'
      - nest: if True, ordering scheme is NEST. Default: False (RING)
      - min: the minimum range value
      - max: the maximum range value
      - flip: 'astro' (default, east towards left, west towards right) or 'geo'
      - remove_dip: if True, remove the dipole+monopole
      - remove_mono: if True, remove the monopole
      - gal_cut: galactic cut for the dipole/monopole fit
      - format: the format of the scale. Default: '%.3g'
      - hold: if True, replace the current Axes by a MollweideAxes.
              use this if you want to have multiple maps on the same
              figure. Default: False
      - sub: use a part of the current figure (same syntax as subplot).
             Default: None
      - margins: either None, or a sequence (left,bottom,right,top)
                 giving the margins on left,bottom,right and top
                 of the axes. Values are relative to figure (0-1).
                 Default: None
      - notext: True: do not add resolution info text
                Default=False
    """
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = npy.zeros(12)+npy.inf
            cbar=False
        if not (hold or sub):
            f=pylab.figure(fig,figsize=(8.5,5.4))
            if not margins:
                    margins = (0.075,0.05,0.075,0.05)
            extent = (0.0,0.0,1.0,1.0)
        elif hold:
            f=pylab.gcf()
            left,bottom,right,top = npy.array(pylab.gca().get_position()).ravel()
            if not margins:
                margins = (0.0,0.0,0.0,0.0)
            extent = (left,bottom,right-left,top-bottom)
            f.delaxes(pylab.gca())
        else: # using subplot syntax
            f=pylab.gcf()
            if hasattr(sub,'__len__'):
                nrows, ncols, idx = sub
            else:
                nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
            if idx < 1 or idx > ncols*nrows:
                raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                     ncols,
                                                                     idx))
            c,r = (idx-1)%ncols,(idx-1)/ncols
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
        if zat and rot:
            raise ValueError('Only give rot or zat, not both')
        if zat:
            rot = npy.array(zat,dtype=npy.float64)
            rot.resize(3)
            rot[1] -= 90
        ax=PA.HpxCartesianAxes(f,extent,coord=coord,rot=rot,
                               format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)
        ax.projmap(map,nest=nest,coord=coord,vmin=min,vmax=max,
                   xsize=xsize,ysize=ysize,lonra=lonra,latra=latra,
                   cmap=cmap,norm=norm,aspect=aspect)
        if cbar:
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(ax.get_images()[0],ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1)
            else:
                cb=f.colorbar(ax.get_images()[0],orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1)
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

def graticule(dpar=None,dmer=None,coord=None,local=None,**kwds):
    """Create a graticule, either on an existing mollweide map or not.

    Parameters:
      - dpar, dmer: interval in degrees between meridians and between parallels
      - coord: the coordinate system of the graticule (make rotation if needed,
               using coordinate system of the map if it is defined)
      - local: True if local graticule (no rotation is performed)
    Return:
      None
    """
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f = pylab.gcf()
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
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
graticule.__doc__ = PA.SphericalProjAxes.graticule.__doc__
    
def delgraticules():
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.delgraticules()
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
delgraticules.__doc__ = PA.SphericalProjAxes.delgraticules.__doc__

def projplot(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret = None
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projplot(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
    return ret
projplot.__doc__ = PA.SphericalProjAxes.projplot.__doc__
    
def projscatter(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret=None
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projscatter(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
    return ret
projscatter.__doc__ = PA.SphericalProjAxes.projscatter.__doc__

def projtext(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    ret = None
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ret = ax.projtext(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
    return ret
projtext.__doc__ = PA.SphericalProjAxes.projtext.__doc__



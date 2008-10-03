import projaxes as PA
import pylab
import numpy as npy
import matplotlib
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import pixelfunc

def mollview(map,fig=None,rot=None,coord=None,unit='',
             xsize=800,title='Mollweide view',nest=False,
             min=None,max=None,remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',cbar=True,cmap=None,
             norm=None):
    """Plot an healpix map (given as an array) in Mollweide projection.
    
    Input:
      - map : an ndarray containing the map
    Parameters:
      - fig: a figure number. Default: create a new figure
      - rot: rotation, either 1,2 or 3 angles describing the rotation
             Default: None
      - coord: either one of 'G', 'E' or 'C' to describe the coordinate
               systm of the map, or a sequence of 2 of these to make
               rotation from the first to the second coordinate system.
               Default: None
      - unit: a text describing the unit. Default: ''
      - xsize: the size of the image. Default: 800
      - title: the title of the plot. Default: 'Mollweide view'
      - nest: if True, ordering scheme is NEST. Default: False (RING)
      - min: the minimum range value
      - max: the maximum range value
      - remove_dip: if True, remove the dipole+monopole
      - remove_mono: if True, remove the monopole
      - gal_cut: galactic cut for the dipole/monopole fit
      - format: the format of the scale. Default: '%g'
    """
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f=pylab.figure(fig,figsize=(8.5,5.4))
        ax=PA.HpxMollweideAxes(f,(0.02,0.05,0.96,0.9),coord=coord,rot=rot,
                               format=format)
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
                              pad=0.05,fraction=0.1,boundaries=b,values=v)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(ax.get_images()[0],orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v)
        ax.set_title(title)
        ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            pylab.show()


def gnomview(map,fig=None,rot=None,coord=None,unit='',
             xsize=200,ysize=None,reso=1.5,degree=False,
             title='Gnomonic view',nest=False,remove_dip=False,
             remove_mono=False,gal_cut=0,
             min=None,max=None,format='%g',cbar=True,
             cmap=None, norm=None):
    """Plot an healpix map (given as an array) in Gnomonic projection.

    Input:
      - map : an ndarray containing the map
    Parameters:
      - fig: a figure number. Default: create a new figure
      - rot: rotation, either 1,2 or 3 angles describing the rotation
             Default: None
      - coord: either one of 'G', 'E' or 'C' to describe the coordinate
               systm of the map, or a sequence of 2 of these to make
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
      - remove_dip: if True, remove the dipole+monopole
      - remove_mono: if True, remove the monopole
      - gal_cut: galactic cut for the dipole/monopole fit
      - format: the format of the scale. Default: '%.3g'
    """
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f=pylab.figure(fig,figsize=(5.5,6))
        ax=PA.HpxGnomonicAxes(f,(0.0,0.05,1.0,0.85),coord=coord,rot=rot,
                              format=format)
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
            pylab.show()


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
            pylab.show()
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
            pylab.show()
delgraticules.__doc__ = PA.SphericalProjAxes.delgraticules.__doc__

def projplot(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.projplot(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            pylab.show()
projplot.__doc__ = PA.SphericalProjAxes.projplot.__doc__
    
def projscatter(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.projscatter(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            pylab.show()
projscatter.__doc__ = PA.SphericalProjAxes.projscatter.__doc__

def projtext(*args,**kwds):
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        f = pylab.gcf()
        for ax in f.get_axes():
            if isinstance(ax,PA.SphericalProjAxes):
                ax.projtext(*args,**kwds)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            pylab.show()
projtext.__doc__ = PA.SphericalProjAxes.projtext.__doc__



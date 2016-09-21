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
from . import projector as P
from . import rotator as R
from . import pixelfunc
import matplotlib
import matplotlib.axes
import numpy as np
import six

from ._healpy_pixel_lib import UNSEEN

pi = np.pi
dtor = pi/180.

class SphericalProjAxes(matplotlib.axes.Axes):
    """Define a special Axes to take care of spherical projection.

    Parameters
    ----------
    projection : a SphericalProj class or a class derived from it.
        type of projection
    rot : list or string
        define rotation. See rotator.
    coord : list or string
        define coordinate system. See rotator.
    coordprec : number of digit after floating point for coordinates display.
    format : format string for value display.
      
    Notes
    -----
    Other keywords from Axes (see Axes).
    """
    def __init__(self, ProjClass, *args, **kwds):
        if not issubclass(ProjClass, P.SphericalProj):
            raise TypeError("First argument must be a SphericalProj class "
                            "(or derived from)")
        self.proj = ProjClass(rot   = kwds.pop('rot',None),
                              coord = kwds.pop('coord',None),
                              flipconv = kwds.pop('flipconv',None),
                              **kwds.pop('arrayinfo', {}))
        kwds.setdefault('format','%g')
        kwds.setdefault('coordprec',2)
        kwds['aspect'] = 'equal'
        super(SphericalProjAxes,self).__init__(*args, **kwds)
        self.axis('off')
        self.set_autoscale_on(False)
        xmin,xmax,ymin,ymax = self.proj.get_extent()
        self.set_xlim(xmin,xmax)
        self.set_ylim(ymin,ymax)
        dx,dy = self.proj.ang2xy(pi/2.,1.*dtor,direct=True)
        self._segment_threshold = 16.*np.sqrt(dx**2+dy**2)
        self._segment_step_rad = 0.1*pi/180
        self._do_border = True
        self._gratdef = {}
        self._gratdef['local'] = False
        self._gratdef['dpar'] = 30.
        
    def set_format(self,f):
        """Set the format string for value display
        """
        self._format=f
        return f
    
    def set_coordprec(self,n):
        """Set the number of digits after floating point for coord display.
        """
        self._coordprec = n

    def format_coord(self,x,y):
        """Format the coordinate for display in status bar. Take projection
        into account.
        """
        format=self._format+' at '
        pos=self.get_lonlat(x,y)
        if pos is None or np.isnan(pos).any(): return ''
        lon,lat = np.around(pos,decimals=self._coordprec)
        val = self.get_value(x,y)
        if val is None:
            format = '%s'
            val = ''
        elif type(val) is str: format='%s @ '
        coordsys = self.proj.coordsysstr
        if coordsys != '':
            res=(format+'(%g, %g) in %s')%(val,lon,lat,
                                           coordsys[0:3])
        else:
            res=(format+'lon=%g, lat=%g')%(val,lon,lat)
        return res

    def get_lonlat(self,x,y):
        """Get the coordinate in the coord system of the image, in lon/lat in deg.
        """
        lon,lat = self.proj.xy2ang(x,y,lonlat=True)
        return lon,lat

    def get_value(self,x,y):
        """Get the value of the map at position x,y
        """
        if len(self.get_images()) < 1:
            return None
        im = self.get_images()[-1]
        arr=im.get_array()
        i,j = self.proj.xy2ij(x,y)
        if i is None or j is None:
            return None
        elif arr.mask is not np.ma.nomask and arr.mask[i,j]:
            return 'UNSEEN'
        else:
            return arr[i,j]

    def projmap(self,map,vec2pix_func,vmin=None,vmax=None,badval=UNSEEN,
                cmap=None,norm=None,rot=None,coord=None,**kwds):
        """Project a map on the SphericalProjAxes.

        Parameters
        ----------
        map : array-like
          The map to project.
        vec2pix_func : function
          The function describing the pixelisation.
        vmin, vmax : float, scalars
          min and max value to use instead of min max of the map
        badval : float
          The value of the bad pixels
        cmap : a color map
          The colormap to use (see matplotlib.cm)
        rot : sequence
          In the form (lon, lat, psi) (unit: degree):the center of the map is
          at (lon, lat) and rotated by angle psi around that direction.
        coord : {'G', 'E', 'C', None}
          The coordinate system of the map ('G','E' or 'C'), rotate
          the map if different from the axes coord syst.
        
        Notes
        -----
        Other keywords are transmitted to :func:`matplotlib.Axes.imshow`
        """
        img = self.proj.projmap(map,vec2pix_func,rot=rot,coord=coord)
        w = ~( np.isnan(img) | 
               np.isinf(img) | 
               pixelfunc.mask_bad(img, badval = badval) )
        try:
            if vmin is None: vmin = img[w].min()
        except ValueError:
            vmin = 0.
        try:
            if vmax is None: vmax = img[w].max()
        except ValueError:
            vmax = 0.
        if vmin > vmax:
            vmin = vmax
        if vmin == vmax:
            vmin -= 1.
            vmax += 1.
        cm,nn = get_color_table(vmin,vmax,img[w],cmap=cmap,norm=norm)
        ext = self.proj.get_extent()
        img = np.ma.masked_values(img, badval)
        aximg = self.imshow(img,extent = ext,cmap=cm,norm=nn,
                            interpolation='nearest',origin='lower',
                            vmin=vmin,vmax=vmax,**kwds)
        xmin,xmax,ymin,ymax = self.proj.get_extent()
        self.set_xlim(xmin,xmax)
        self.set_ylim(ymin,ymax)
        return img

    def projplot(self,*args,**kwds):
        """projplot is a wrapper around :func:`matplotlib.Axes.plot` to take into account the
        spherical projection.

        You can call this function as::
        
           projplot(theta, phi)        # plot a line going through points at coord (theta, phi)
           projplot(theta, phi, 'bo')  # plot 'o' in blue at coord (theta, phi)
           projplot(thetaphi)          # plot a line going through points at coord (thetaphi[0], thetaphi[1])
           projplot(thetaphi, 'bx')    # idem but with blue 'x'
        
        Parameters
        ----------
        theta, phi : float, array-like
          Coordinates of point to plot. Can be put into one 2-d array, first line is
          then *theta* and second line is *phi*. See *lonlat* parameter for unit.
        fmt : str
          A format string (see :func:`matplotlib.Axes.plot` for details)
        lonlat : bool, optional
          If True, theta and phi are interpreted as longitude and latitude
          in degree, otherwise, as colatitude and longitude in radian
        coord : {'E', 'G', 'C', None}
          The coordinate system of the points, only used if the coordinate
          coordinate system of the Axes has been defined and in this
          case, a rotation is performed
        rot : None or sequence
          rotation to be applied =(lon, lat, psi) : lon, lat will be position of the
          new Z axis, and psi is rotation around this axis, all in degree.
          if None, no rotation is performed
        direct : bool
          if True, the rotation to center the projection is not
          taken into account

        Notes
        -----
        Other keywords are passed to :func:`matplotlib.Axes.plot`.

        See Also
        --------
        projscatter, projtext
        """
        fmt = None
        if len(args) < 1:
            raise ValueError("No argument given")
        if len(args) == 1:
            theta,phi = np.asarray(args[0])
        elif len(args) == 2:
            if type(args[1]) is str:
                fmt=args[1]
                theta,phi = np.asarray(args[0])
            else:
                theta,phi = np.asarray(args[0]),np.asarray(args[1])
        elif len(args) == 3:
            if type(args[2]) is not str:
                raise TypeError("Third argument must be a string")
            else:
                theta,phi = np.asarray(args[0]),np.asarray(args[1])
                fmt = args[2]
        else:
            raise TypeError("Three args maximum")
        rot=kwds.pop('rot',None)
        if rot is not None:
            rot = np.array(np.atleast_1d(rot),copy=1)
            rot.resize(3)
            rot[1] = rot[1]-90.
        coord=self.proj.mkcoord(kwds.pop('coord',None))[::-1]
        lonlat=kwds.pop('lonlat',False)
        vec = R.dir2vec(theta,phi,lonlat=lonlat)
        vec = (R.Rotator(rot=rot,coord=coord,eulertype='Y')).I(vec)
        x,y = self.proj.vec2xy(vec,direct=kwds.pop('direct',False))
        x,y = self._make_segment(x,y,threshold=kwds.pop('threshold',
                                                        self._segment_threshold))
        thelines = []
        for xx,yy in zip(x,y):
            if fmt is not None:
                try: # works in matplotlib 1.3 and earlier
                    linestyle, marker, color = matplotlib.axes._process_plot_format(fmt)
                except: # matplotlib 1.4 and later
                    linestyle, marker, color = matplotlib.axes._axes._process_plot_format(fmt)
                kwds.setdefault('linestyle',linestyle)
                kwds.setdefault('marker',marker)
                if color is not None: kwds.setdefault('color',color)
            l = matplotlib.lines.Line2D(xx,yy,**kwds)
            self.add_line(l)
            thelines.append(l)
        return thelines

    def projscatter(self,theta, phi=None,*args,**kwds):
        """Projscatter is a wrapper around :func:`matplotlib.Axes.scatter` to take into account the
        spherical projection.
        
        You can call this function as::
        
           projscatter(theta, phi)     # plot points at coord (theta, phi)
           projplot(thetaphi)          # plot points at coord (thetaphi[0], thetaphi[1])

        Parameters
        ----------
        theta, phi : float, array-like
          Coordinates of point to plot. Can be put into one 2-d array, first line is
          then *theta* and second line is *phi*. See *lonlat* parameter for unit.
        lonlat : bool, optional
          If True, theta and phi are interpreted as longitude and latitude
          in degree, otherwise, as colatitude and longitude in radian
        coord : {'E', 'G', 'C', None}, optional
          The coordinate system of the points, only used if the coordinate
          coordinate system of the axes has been defined and in this
          case, a rotation is performed
        rot : None or sequence, optional
          rotation to be applied =(lon, lat, psi) : lon, lat will be position of the
          new Z axis, and psi is rotation around this axis, all in degree.
          if None, no rotation is performed
        direct : bool, optional
          if True, the rotation to center the projection is not
          taken into account

        Notes
        -----
        Other keywords are passed to :func:`matplotlib.Axes.plot`.

        See Also
        --------
        projplot, projtext
        """
        save_input_data = hasattr(self.figure, 'zoomtool')
        if save_input_data:
            input_data = (theta, phi, args, kwds.copy())
        if phi is None:
            theta,phi = np.asarray(theta)
        else:
            theta, phi = np.asarray(theta), np.asarray(phi)
        rot=kwds.pop('rot',None)
        if rot is not None:
            rot = np.array(np.atleast_1d(rot),copy=1)
            rot.resize(3)
            rot[1] = rot[1]-90.
        coord=self.proj.mkcoord(kwds.pop('coord',None))[::-1]
        lonlat=kwds.pop('lonlat',False)
        vec = R.dir2vec(theta,phi,lonlat=lonlat)
        vec = (R.Rotator(rot=rot,coord=coord,eulertype='Y')).I(vec)
        x,y = self.proj.vec2xy(vec,direct=kwds.pop('direct',False))
        s = self.scatter(x, y, *args, **kwds)
        if save_input_data:
            if not hasattr(self, '_scatter_data'):
                self._scatter_data = []
            self._scatter_data.append((s, input_data))
        return s

    def projtext(self, theta, phi, s, **kwds):
        """Projtext is a wrapper around :func:`matplotlib.Axes.text` to take into account the
        spherical projection.

        Parameters
        ----------
        theta, phi : float, array-like
          Coordinates of point to plot. Can be put into one 2-d array, first line is
          then *theta* and second line is *phi*. See *lonlat* parameter for unit.
        text : str
          The text to be displayed.
        lonlat : bool, optional
          If True, theta and phi are interpreted as longitude and latitude
          in degree, otherwise, as colatitude and longitude in radian
        coord : {'E', 'G', 'C', None}, optional
          The coordinate system of the points, only used if the coordinate
          coordinate system of the axes has been defined and in this
          case, a rotation is performed
        rot : None or sequence, optional
          rotation to be applied =(lon, lat, psi) : lon, lat will be position of the
          new Z axis, and psi is rotation around this axis, all in degree.
          if None, no rotation is performed
        direct : bool, optional
          if True, the rotation to center the projection is not
          taken into account

        Notes
        -----
        Other keywords are passed to :func:`matplotlib.Axes.text`.

        See Also
        --------
        projplot, projscatter
        """
        if phi is None:
            theta,phi = np.asarray(theta)
        else:
            theta, phi = np.asarray(theta), np.asarray(phi)
        rot=kwds.pop('rot',None)
        if rot is not None:
            rot = np.array(np.atleast_1d(rot),copy=1)
            rot.resize(3)
            rot[1] = rot[1]-90.
        coord=self.proj.mkcoord(kwds.pop('coord',None))[::-1]
        lonlat=kwds.pop('lonlat',False)
        vec = R.dir2vec(theta,phi,lonlat=lonlat)
        vec = (R.Rotator(rot=rot,coord=coord,eulertype='Y')).I(vec)
        x,y = self.proj.vec2xy(vec,direct=kwds.pop('direct',False))
        return self.text(x,y,s,**kwds)

    def _make_segment(self,x,y,threshold=None):
        if threshold is None:
            threshold = self._segment_threshold
        x,y=np.atleast_1d(x),np.atleast_1d(y)
        d2 = np.sqrt((np.roll(x,1)-x)**2+(np.roll(y,1)-y)**2)
        w=np.where(d2 > threshold)[0]
        #w=w[w!=0]
        xx=[]
        yy=[]
        if len(w) == 1:
            x=np.roll(x,-w[0])
            y=np.roll(y,-w[0])
            xx.append(x)
            yy.append(y)
        elif len(w) >= 2:
            xx.append(x[0:w[0]])
            yy.append(y[0:w[0]])
            for i in six.moves.xrange(len(w)-1):
                xx.append(x[w[i]:w[i+1]])
                yy.append(y[w[i]:w[i+1]])
            xx.append(x[w[-1]:])
            yy.append(y[w[-1]:])
        else:
            xx.append(x)
            yy.append(y)
        return xx,yy

    def get_parallel_interval(self,vx,vy=None,vz=None):
        """Get the min and max value of theta of the parallel to cover the
        field of view.

        Input:
          - the normalized vector of the direction of the center of the
            projection, in the reference frame of the graticule.
        Return:
          - vmin,vmax : between 0 and pi, vmin<vmax, the interval of theta
                        for the parallels crossing the field of view
        """
        if vy is None and vz is None:
            vx,vy,vz = vx
        elif vy is None or vz is None:
            raise ValueError("Both vy and vz must be given or both not given")
        a = np.arccos(vz)
        fov = self.proj.get_fov()
        vmin = max(0., a-fov/2.)
        vmax = min(pi, a+fov/2.)
        return vmin,vmax

    def get_meridian_interval(self, vx, vy=None, vz=None):
        """Get the min and max value of phi of the meridians to cover the field
        of view.

        Input:
          - the normalized vector of the direction of the center of the
            projection, in the reference frame of the graticule.
        Return:
          - vmin,vmax : the interval of phi for the
                        meridians crossing the field of view.
        """
        if vy is None and vz is None:
            vx,vy,vz = vx
        elif vy is None or vz is None:
            raise ValueError("Both vy and vz must be given or both not given")
        fov = self.proj.get_fov()
        th = np.arccos(vz)
        if th <= fov/2.: # test whether north pole is visible
            return -np.pi,np.pi
        if abs(th-pi) <= fov/2.: # test whether south pole is visible
            return -np.pi,np.pi
        sth = np.sin(th)
        phi0 = np.arctan2(vy,vx)
        return phi0 - fov/sth/2., phi0 + fov/sth/2.

    def graticule(self,dpar=None,dmer=None,coord=None,local=None,verbose=True,**kwds):
        """Draw a graticule.
        
        Input:
         - dpar: angular separation between parallels in degree
         - dmer: angular separation between meridians in degree
         - coord: coordinate system of the graticule ('G', 'E' or 'C')
         - local: if True, no rotation performed at all
        """
        gratargs = (dpar,dmer,coord,local)
        gratkwds = kwds
        if dpar is None: dpar=self._gratdef['dpar']
        if local is None: local=self._gratdef['local']
        if dmer is None: dmer = dpar
        dpar = abs(dpar)*dtor
        dmer = abs(dmer)*dtor
        if not local:
            vec = R.dir2vec(self.proj.get_center())
            vec0 = R.Rotator(coord=self.proj.mkcoord(coord=coord)).I(vec)
        else:
            vec = (1,0,0)
            vec0 = (1,0,0)
        u_pmin,u_pmax = kwds.pop('pmax',None),kwds.pop('pmin',None)
        u_mmin,u_mmax = kwds.pop('mmin',None),kwds.pop('mmax',None)
        if u_pmin: u_pmin = (pi/2.-u_pmin*dtor)%pi
        if u_pmax: u_pmax = (pi/2.-u_pmax*dtor)%pi
        if u_mmin: u_mmin = ( ((u_mmin+180.)%360)-180)*dtor
        if u_mmax: u_mmax = ( ((u_mmax+180.)%360)-180)*dtor
        pmin,pmax = self.get_parallel_interval(vec0)
        mmin,mmax = self.get_meridian_interval(vec0)
        if u_pmin: pmin = u_pmin
        if u_pmax: pmax = u_pmax
        if u_mmin: mmin = u_mmin
        if u_mmax: mmax = u_pmax
        if verbose:
            print('{0} {1} {2} {3}'.format(
                pmin/dtor, pmax/dtor, mmin/dtor, mmax/dtor))
        if not kwds.pop('force',False):
            dpar,dmer = self._get_interv_graticule(pmin,pmax,dpar,
                                                   mmin,mmax,dmer,
                                                   verbose=verbose)
        theta_list = np.around(np.arange(pmin,pmax+0.5*dpar,dpar)/dpar)*dpar
        phi_list = np.around(np.arange(mmin,mmax+0.5*dmer,dmer)/dmer)*dmer
        theta = np.arange(pmin,pmax,min((pmax-pmin)/100.,
                                         self._segment_step_rad))
        phi = np.arange(mmin,mmax,min((mmax-mmin)/100.,
                                       self._segment_step_rad))
        equator = False
        gratlines = []
        kwds.setdefault('lw',1)
        kwds.setdefault('color','k')
        for t in theta_list:
            if abs(t-pi/2.)<1.e-10:
                fmt = '-'
                equator=True
            elif abs(t) < 1.e-10: # special case: north pole
                t = 1.e-10 
                fmt = '-'
            elif abs(t-pi) < 1.e-10: # special case: south pole
                t = pi-1.e-10 
                fmt = '-'
            else:
                fmt =':'
            gratlines.append(self.projplot(phi*0.+t, phi,fmt,
                                           coord=coord,
                                           direct=local,**kwds))
        if not equator and pmin <= pi/2. and pi/2 <= pmax:
            gratlines.append(self.projplot(phi*0.+pi/2., phi,'-',
                                           coord=coord,
                                           direct=local,**kwds))
        for p in phi_list:
            if abs(p)<1.e-10: fmt = '-'
            else: fmt =':'
            gratlines.append(self.projplot(theta, theta*0.+p,fmt,
                                           coord=coord,
                                           direct=local,**kwds))
        # Now the borders (only useful for full sky projection)
        if hasattr(self,'_do_border') and self._do_border:
            theta = np.arange(0,181)*dtor
            gratlines.append(self.projplot(theta, theta*0-pi,'-k',
                                           lw=1,direct=True))
            gratlines.append(self.projplot(theta, theta*0+0.9999*pi,'-k',
                                           lw=1,direct=True))
            phi = np.arange(-180,180)*dtor
            gratlines.append(self.projplot(phi*0+1.e-10, phi,'-k',
                                           lw=1,direct=True))
            gratlines.append(self.projplot(phi*0+pi-1.e-10, phi,'-k',
                                           lw=1,direct=True))
        if hasattr(self,'_graticules'):
            self._graticules.append((gratargs,gratkwds,gratlines))
        else:
            self._graticules = [(gratargs,gratkwds,gratlines)]
        return dpar,dmer
    
    def delgraticules(self):
        """Delete all graticules previously created on the Axes.
        """
        if hasattr(self,'_graticules'):
            for dum1,dum2,g in self._graticules:
                for gl in g:
                    for l in gl: 
                        if l in self.lines:
                            self.lines.remove(l)
                        else:
                            print('line not in lines')
            del self._graticules

    def _get_interv_graticule(self,pmin,pmax,dpar,mmin,mmax,dmer,verbose=True):
        def set_prec(d,n,nn=2):
            arcmin=False
            if d/n < 1.:
                d *= 60
                arcmin = True
                nn = 1
            x = d/n
            y = nn*x
            ex = np.floor(np.log10(y))
            z = np.around(y/10**ex)*10**ex/nn
            if arcmin:
                z = 1./np.around(60./z)
            return z
        max_n_par = 18
        max_n_mer = 36
        n_par = (pmax-pmin)/dpar
        n_mer = (mmax-mmin)/dmer
        if n_par > max_n_par:
            dpar = set_prec((pmax-pmin)/dtor,max_n_par/2)*dtor
        if n_mer > max_n_mer:
            dmer = set_prec((mmax-mmin)/dtor,max_n_mer/2,nn=1)*dtor
        if dmer/dpar < 0.2 or dmer/dpar > 5.:
            dmer = dpar = max(dmer,dpar)
        vdeg = int(np.floor(np.around(dpar/dtor,10)))
        varcmin = (dpar/dtor-vdeg)*60.
        if verbose: print("The interval between parallels is {0:d} deg {1:.2f}'.".format(vdeg,varcmin))
        vdeg = int(np.floor(np.around(dmer/dtor,10)))
        varcmin = (dmer/dtor-vdeg)*60.
        if verbose: print("The interval between meridians is {0:d} deg {1:.2f}'.".format(vdeg,varcmin))
        return dpar,dmer
        
class GnomonicAxes(SphericalProjAxes):
    """Define a gnomonic Axes to handle gnomonic projection.

    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',3)
        super(GnomonicAxes,self).__init__(P.GnomonicProj, *args,**kwds)
        self._do_border = False
        self._gratdef['local'] = True
        self._gratdef['dpar'] = 1.

    def projmap(self,map,vec2pix_func,xsize=200,ysize=None,reso=1.5,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,reso=reso)
        return super(GnomonicAxes,self).projmap(map,vec2pix_func,**kwds)
        
class HpxGnomonicAxes(GnomonicAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        xsize = kwds.pop('xsize',200)
        ysize = kwds.pop('ysize',None)
        reso = kwds.pop('reso',1.5)
        return super(HpxGnomonicAxes,self).projmap(map,f,xsize=xsize,
                                            ysize=ysize,reso=reso,**kwds)


class MollweideAxes(SphericalProjAxes):
    """Define a mollweide Axes to handle mollweide projection.

    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(MollweideAxes,self).__init__(P.MollweideProj, *args,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)

    def projmap(self,map,vec2pix_func,xsize=800,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize)
        img = super(MollweideAxes,self).projmap(map,vec2pix_func,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)
        return img
        
class HpxMollweideAxes(MollweideAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxMollweideAxes,self).projmap(map,f,**kwds)

class CartesianAxes(SphericalProjAxes):
    """Define a cylindrical Axes to handle cylindrical projection.
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(CartesianAxes,self).__init__(P.CartesianProj, *args, **kwds)
        self._segment_threshold = 180
        self._segment_step_rad = 0.1*pi/180
        self._do_border = True
        
    def projmap(self,map,vec2pix_func,xsize=800,ysize=None,lonra=None,latra=None,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,lonra=lonra,latra=latra)
        return super(CartesianAxes,self).projmap(map,vec2pix_func,**kwds)
        
class HpxCartesianAxes(CartesianAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxCartesianAxes,self).projmap(map,f,**kwds)

        
class OrthographicAxes(SphericalProjAxes):
    """Define an orthographic Axes to handle orthographic projection.
    
    Input:
    - rot=, coord= : define rotation and coordinate system. See rotator.
    - coordprec= : num of digits after floating point for coordinates display.
    - format= : format string for value display.
    
    Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(OrthographicAxes,self).__init__(P.OrthographicProj, *args,**kwds)
        self._segment_threshold = 0.01
        self._do_border = False
        
    def projmap(self,map,vec2pix_func,xsize=800,half_sky=False,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,half_sky=half_sky)
        img = super(OrthographicAxes,self).projmap(map,vec2pix_func,**kwds)
        if half_sky: ratio = 1.01
        else: ratio = 2.01
        self.set_xlim(-ratio,ratio)
        self.set_ylim(-1.01,1.01)
        return img
        
class HpxOrthographicAxes(OrthographicAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(len(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxOrthographicAxes,self).projmap(map,f,**kwds)

class AzimuthalAxes(SphericalProjAxes):
    """Define an Azimuthal Axes to handle azimuthal equidistant or
       Lambert azimuthal equal-area projections.

    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.

      Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',3)
        super(AzimuthalAxes,self).__init__(P.AzimuthalProj, *args,**kwds)
        self._do_border = False

    def projmap(self,map,vec2pix_func,xsize=200,ysize=None,reso=1.5,lamb=True,half_sky=False,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,reso=reso,lamb=lamb,half_sky=half_sky)
        return super(AzimuthalAxes,self).projmap(map,vec2pix_func,**kwds)

class HpxAzimuthalAxes(AzimuthalAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        xsize = kwds.pop('xsize',800)
        ysize = kwds.pop('ysize',None)
        reso = kwds.pop('reso',1.5)
        lamb = kwds.pop('lamb',True)
        return super(HpxAzimuthalAxes,self).projmap(map,f,xsize=xsize,
                                            ysize=ysize,reso=reso,lamb=lamb,**kwds)

###################################################################
#
#   Table color for mollview, gnomview, and orthview.
#   Currently defined for so that the default colormap, found in 
#   matplotlib.rcParams['image.cmap'], the data is displayed with
#   values greater than vmax as the final element of the colormap,
#   masked indices gray, and the background set to white.
#
#   With matplotlib.rcParams['image.cmap'] assigned to a string
#   corresponding to a standard matplotlib colormap, one can call
#   hp.mollview(m) and have the map projected in the standard way,
#   whereas using just, e.g., hp.mollview(m, cmap='jet') will display
#   the data with a non-white background.
#
#   One can set the default colormap in the matplotlibrc file, or set
#   it in situ:
#   >>> matplotlib.rcParam['image.cmap'] = 'coolwarm'
#   >>> hp.mollview(m)
#   Note that custom colormaps can also be used, but they need to be 
#   registered ahead fo time, as shown in
#   http://matplotlib.org/examples/pylab_examples/custom_cmap.html

def get_color_table(vmin,vmax,val,cmap=None,norm=None):
    # Create color table
    newcmap = create_colormap(cmap)
    if type(norm) is str:
        if norm.lower().startswith('log'):
            norm = LogNorm2(clip=False)
        elif norm.lower().startswith('hist'):
            norm = HistEqNorm(clip=False)
        else:
            norm = None
    if norm is None:
        norm = LinNorm2(clip=False)

    norm.vmin = vmin
    norm.vmax = vmax
    norm.autoscale_None(val)
    
    return newcmap,norm

def create_colormap(cmap):
    if cmap is not None:
        return cmap 
    cmap0 = matplotlib.cm.get_cmap(matplotlib.rcParams['image.cmap'])
    if hasattr(cmap0, '_segmentdata'):
        newcm = matplotlib.colors.LinearSegmentedColormap('newcm',cmap0._segmentdata,
                                                    cmap0.N)
    else:
        newcm = cmap0
    newcm.set_over(newcm(1.0))
    newcm.set_under('w')
    newcm.set_bad('gray')
    return newcm

##################################################################
#
#   A Locator that gives the bounds of the interval
#
class BoundaryLocator(matplotlib.ticker.Locator):
    def __init__(self,N=2):
        if N < 2:
            raise ValueError("Number of locs must be greater than 1")
        self.Nlocs=N

    def __call__(self):
        if matplotlib.__version__ < '0.98':
            vmin,vmax = self.viewInterval.get_bounds()
        else:
            vmin, vmax = self.axis.get_view_interval()
        locs = vmin + np.arange(self.Nlocs)*(vmax-vmin)/(self.Nlocs-1.)
        return locs

    def autoscale(self):
        self.verify_intervals()
        vmin,vmax = self.dataInterval.get_bounds()
        if vmax<vmin:
            vmin,vmax = vmax,vmin
        if vmin==vmax:
            vmin -= 1
            vmax += 1
        return vmin,vmax


##################################################################
#
#   A normalization class to get color table equalised by
#   the histogram of data
#

class HistEqNorm(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, clip=False):
        matplotlib.colors.Normalize.__init__(self,vmin,vmax,clip)
        self.xval = None
        self.yval = None

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = 'array'
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = np.ma.array([value]).astype(np.float)

        self.autoscale_None(val)

        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin==vmax:
            return 0.0 * val
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax),
                                   mask=mask)
            result = np.ma.array(np.interp(val, self.xval, self.yval),
                                  mask=np.ma.getmask(val))
            result[np.isinf(val.data)] = -np.inf
        if vtype == 'scalar':
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        if matplotlib.cbook.iterable(value):
            vtype='array'
            val = np.ma.array(value)
        else:
            vtype='scalar'
            val = np.ma.array([value])
        result = np.ma.array(self._lininterp(val, self.yval, self.xval),
                              mask=np.ma.getmask(val))
        result[np.isinf(val.data)] = -np.inf
        if vtype == 'scalar':
            result = result[0]
        return result
        
    def autoscale_None(self, val):
        changed = False
        if self.vmin is None:
            self.vmin = val.min()
            changed = True
        if self.vmax is None:
            self.vmax = val.max()
            changed = True
        if changed or self.xval is None or self.yval is None:
            self._set_xyvals(val)
            
    def autoscale(self, val):
        self.vmin = val.min()
        self.vmax = val.max()
        self._set_xyvals(val)

    def _set_xyvals(self,val):
        data = np.ma.asarray(val).ravel()
        w=np.isinf(data.data)
        if data.mask is not np.ma.nomask:
            w = w|data.mask
        data2 = data.data[~w]
        bins = min(data2.size//20, 5000)
        if bins < 3: bins=data2.size
        try:
            # for numpy 1.1, use new bins format (left and right edges)
            hist, bins = np.histogram(data2,bins=bins,
                                       range=(self.vmin,self.vmax),
                                       new=True)
        except TypeError:
            # for numpy <= 1.0 or numpy >= 1.2, no new keyword
            hist, bins = np.histogram(data2,bins=bins,
                                       range=(self.vmin,self.vmax))
        if bins.size == hist.size+1:
            # new bins format, remove last point
            bins = bins[:-1]
        hist = hist.astype(np.float)/np.float(hist.sum())
        self.yval = np.concatenate([[0.], hist.cumsum(), [1.]])
        self.xval = np.concatenate([[self.vmin], bins + 0.5*(bins[1]-bins[0]),
                                    [self.vmax]])

    def _lininterp(self,x,X,Y):
        if hasattr(x,'__len__'):
            xtype = 'array'
            xx=np.asarray(x).astype(np.float)
        else:
            xtype = 'scalar'
            xx=np.asarray([x]).astype(np.float)
        idx = X.searchsorted(xx)
        yy = xx*0
        yy[idx>len(X)-1] = Y[-1]  # over
        yy[idx<=0] = Y[0]          # under
        wok = np.where((idx>0) & (idx<len(X)))  # the good ones
        iok=idx[wok]
        yywok = Y[iok-1] + ( (Y[iok]-Y[iok-1])/(X[iok]-X[iok-1])
                             * (xx[wok]-X[iok-1]) )
        w = np.where( ((X[iok]-X[iok-1]) == 0) )   # where are the nan ?
        yywok[w] = Y[iok[w]-1]    # replace by previous value
        wl = np.where(xx[wok] == X[0])
        yywok[wl] = Y[0]
        wh = np.where(xx[wok] == X[-1])
        yywok[wh] = Y[-1]
        yy[wok] = yywok
        if xtype == 'scalar':
            yy = yy[0]
        return yy


##################################################################
#
#   A normalization class to get logarithmic color table
#

class LogNorm2(matplotlib.colors.Normalize):
    """
    Normalize a given value to the 0-1 range on a log scale
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = 'array'
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = np.ma.array([value]).astype(np.float)

        val = np.ma.masked_where(np.isinf(val.data),val)

        self.autoscale_None(val)
        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin<=0:
            raise ValueError("values must all be positive")
        elif vmin==vmax:
            return type(value)(0.0 * np.asarray(value))
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax),
                                   mask=mask)
            result = (np.ma.log(val)-np.log(vmin))/(np.log(vmax)-np.log(vmin))
            result.data[result.data<0]=0.0
            result.data[result.data>1]=1.0
            result[np.isinf(val.data)] = -np.inf
            if result.mask is not np.ma.nomask:
                result.mask[np.isinf(val.data)] = False
        if vtype == 'scalar':
            result = result[0]
        return result

    def autoscale_None(self, A):
        ' autoscale only None-valued vmin or vmax'
        if self.vmin is None or self.vmax is None:
            val = np.ma.masked_where(np.isinf(A.data),A)
            matplotlib.colors.Normalize.autoscale_None(self,val)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = float(self.vmin), float(self.vmax)

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin * np.ma.power((vmax/vmin), val)
        else:
            return vmin * np.pow((vmax/vmin), value)



        
##################################################################
#
#   A normalization class to get linear color table
#

class LinNorm2(matplotlib.colors.Normalize):
    """
    Normalize a given value to the 0-1 range on a lin scale
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        if matplotlib.cbook.iterable(value):
            vtype = 'array'
            val = np.ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = np.ma.array([value]).astype(np.float)

        winf = np.isinf(val.data)
        val = np.ma.masked_where(winf,val)

        self.autoscale_None(val)
        vmin, vmax = float(self.vmin), float(self.vmax)
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin==vmax:
            return type(value)(0.0 * np.asarray(value))
        else:
            if clip:
                mask = np.ma.getmask(val)
                val = np.ma.array(np.clip(val.filled(vmax), vmin, vmax),
                                   mask=mask)
            result = (val-vmin) * (1./(vmax-vmin))
            result.data[result.data<0]=0.0
            result.data[result.data>1]=1.0
            result[winf] = -np.inf
            if result.mask is not np.ma.nomask:
                result.mask[winf] = False
        if vtype == 'scalar':
            result = result[0]
        return result

    def autoscale_None(self, A):
        ' autoscale only None-valued vmin or vmax'
        if self.vmin is None or self.vmax is None:
            val = np.ma.masked_where(np.isinf(A.data),A)
            matplotlib.colors.Normalize.autoscale_None(self,val)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = float(self.vmin), float(self.vmax)

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin + (vmax-vmin) * val
        else:
            return vmin + (vmax-vmin) * value

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
"""This module provides classes for some spherical projection.
To be used when calling SphereProjAxes class.

SphericalProj : a virtual class (do nothing). Just a template for derived 
                (useful) classes

GnomonicProj : Gnomonic projection

AzimuthalProj : Azimuthal equidistant or Lambert azimuthal equal-area projection
"""

from . import rotator as R
import numpy as np
from . import pixelfunc
from .pixelfunc import UNSEEN

pi = np.pi
dtor = np.pi/180.

class SphericalProj(object):
    """
    This class defines functions for spherical projection.
    
    This class contains class method for spherical projection computation. It 
    should not be instantiated. It should be inherited from and methods should
    be overloaded for desired projection.
    """

    name = "None"
    
    def __init__(self, rot=None, coord=None, flipconv=None, **kwds):
        self.rotator  = R.Rotator(rot=rot,  coord=None, eulertype='ZYX')
        self.coordsys = R.Rotator(coord=coord).coordout
        self.coordsysstr = R.Rotator(coord=coord).coordoutstr
        self.set_flip(flipconv)
        self.set_proj_plane_info(**kwds)

    def set_proj_plane_info(self, **kwds):
        allNone = True
        for v in kwds.values():
            if v is not None: allNone = False
        if not allNone:
            self._arrayinfo = dict(kwds)
        else:
            self._arrayinfo = None

    def get_proj_plane_info(self):
        return self._arrayinfo
    arrayinfo = property(get_proj_plane_info,
                         doc="Dictionary with information on the projection array")

    def __eq__(self, a):
        if type(a) is not type(self): return False
        return ( (self.rotator == a.rotator) and
                 (self.coordsys == a.coordsys ) )
    
    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        """From angular direction to position in the projection plane (%s).

        Input:
          - theta: if phi is None, theta[0] contains theta, theta[1] contains phi
          - phi  : if phi is not None, theta,phi are direction
          - lonlat: if True, angle are assumed in degree, and longitude, latitude
          - flipconv is either 'astro' or 'geo'. None will be default.
        Return:
          - x, y: position in %s plane.
        """
        pass
    
    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        """From unit vector direction to position in the projection plane (%s).

        Input:
          - vx: if vy and vz are None, vx[0],vx[1],vx[2] defines the unit vector.
          - vy,vz: if defined, vx,vy,vz define the unit vector
          - lonlat: if True, angle are assumed in degree, and longitude, latitude
          - flipconv is either 'astro' or 'geo'. None will be default.

        Return:
          - x, y: position in %s plane.
        """
        pass
    
    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        """From position in the projection plane to angular direction (%s).

        Input:
          - x : if y is None, x[0], x[1] define the position in %s plane.
          - y : if defined, x,y define the position in projection plane.
          - lonlat: if True, angle are assumed in degree, and longitude, latitude
          - flipconv is either 'astro' or 'geo'. None will be default.

        Return:
          - theta, phi : angular direction.
        """
        pass

    def xy2vec(self, x, y=None, direct=False):
        """From position in the projection plane to unit vector direction (%s).

        Input:
          - x : if y is None, x[0], x[1] define the position in %s plane.
          - y : if defined, x,y define the position in projection plane.
          - lonlat: if True, angle are assumed in degree, and longitude, latitude
          - flipconv is either 'astro' or 'geo'. None will be default.

        Return:
          - theta, phi : angular direction.
        """
        pass
    
    def xy2ij(self, x, y=None):
        """From position in the projection plane to image array index (%s).

        Input:
          - x : if y is None, x[0], x[1] define the position in %s plane.
          - y : if defined, x,y define the position in projection plane.
          - projinfo : additional projection information.

        Return:
          - i,j : image array indices.
        """
        pass
    
    def ij2xy(self, i=None, j=None):
        """From image array indices to position in projection plane (%s).

        Input:
          - if i and j are None, generate arrays of i and j as input
          - i : if j is None, i[0], j[1] define array indices in %s image.
          - j : if defined, i,j define array indices in image.
          - projinfo : additional projection information.

        Return:
          - x,y : position in projection plane.
        """
        pass

    def projmap(self, map, vec2pix_func,rot=None,coord=None):
        """Create an array containing the projection of the map.

        Input:
          - vec2pix_func: a function taking theta,phi and returning pixel number
          - map: an array containing the spherical map to project,
                 the pixelisation is described by vec2pix_func
        Return:
          - a 2D array with the projection of the map.

        Note: the Projector must contain information on the array.
        """
        x,y = self.ij2xy()
        if np.__version__ >= '1.1':
            matype = np.ma.core.MaskedArray
        else:
            matype = np.ma.array
        if type(x) is matype and x.mask is not np.ma.nomask:
            w = (x.mask == False)
        else:
            w = slice(None)
        img=np.zeros(x.shape,np.float64)-np.inf
        vec = self.xy2vec(np.asarray(x[w]),np.asarray(y[w]))
        vec = (R.Rotator(rot=rot,coord=self.mkcoord(coord))).I(vec)
        pix=vec2pix_func(vec[0],vec[1],vec[2])
        # support masked array for map, or a dictionnary (for explicit pixelisation)
        if isinstance(map, matype) and map.mask is not np.ma.nomask:
            mpix = map[pix]
            mpix[map.mask[pix]] = UNSEEN
        elif isinstance(map, dict):
            is_pix_seen = np.in1d(pix, map.keys()).reshape(pix.shape)
            is_pix_unseen = ~is_pix_seen
            mpix = np.zeros_like(img[w])
            mpix[is_pix_unseen] = UNSEEN
            pix_seen = pix[is_pix_seen]
            iterable = (map[p] for p in pix_seen)
            mpix[is_pix_seen] = np.fromiter(iterable, mpix.dtype,
                                             count = pix_seen.size)
        else:
            mpix = map[pix]
        img[w] = mpix
        return img
        
    def set_flip(self, flipconv):
        """flipconv is either 'astro' or 'geo'. None will be default.
        
        With 'astro', east is toward left and west toward right. 
        It is the opposite for 'geo'
        """
        if flipconv is None:
            flipconv = 'astro'  # default
        if    flipconv == 'astro': self._flip = -1
        elif  flipconv == 'geo':   self._flip = 1
        else: raise ValueError("flipconv must be 'astro', 'geo' or None for default.")
        
    def get_extent(self):
        """Get the extension of the projection plane.

        Return:
          extent = (left,right,bottom,top)
        """
        pass

    def get_fov(self):
        """Get the field of view in degree of the plane of projection

        Return:
          fov: the diameter in radian of the field of view
        """
        return 2.*pi

    def get_center(self,lonlat=False):
        """Get the center of the projection.

        Input:
          - lonlat : if True, will return longitude and latitude in degree,
                     otherwise, theta and phi in radian
        Return:
          - theta,phi or lonlat depending on lonlat keyword
        """
        lon, lat = np.asarray(self.rotator.rots[0][0:2])*180/pi
        if lonlat: return lon,lat
        else: return pi/2.-lat*dtor, lon*dtor

    def mkcoord(self,coord):
        if self.coordsys is None:
            return (coord,coord)
        elif coord is None:
            return (self.coordsys,self.coordsys)
        elif type(coord) is str:
            return (coord,self.coordsys)
        else:
            return (tuple(coord)[0],self.coordsys)
        
            
class GnomonicProj(SphericalProj):
    """This class provides class methods for Gnomonic projection.
    """
    
    name = "Gnomonic"

    def __init__(self, rot=None, coord=None, xsize=None, ysize=None, reso=None,
                 **kwds):
        super(GnomonicProj,self).__init__(rot=rot, coord=coord,
                                          xsize=xsize, ysize=ysize,reso=reso,
                                          **kwds)

    def set_proj_plane_info(self, xsize=200,ysize=None,reso=1.5):
        if xsize is None: xsize=200
        if ysize is None: ysize=xsize
        if reso is None: reso=1.5
        super(GnomonicProj,self).set_proj_plane_info(xsize=xsize,
                                                     ysize=ysize,reso=reso)
    
    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct: vec = self.rotator(vx,vy,vz)
        elif vy is None and vz is None: vec=vx
        elif vy is not None and vz is not None: vec=vx,vy,vz
        else: raise ValueError("vy and vz must be both defined or both not defined")
        flip = self._flip
        mask = (np.asarray(vec[0])<=0.)
        w = np.where(mask == False)
        if not mask.any(): mask=np.ma.nomask
        if not hasattr(vec[0],'__len__'):
            if mask is not np.ma.nomask:
                x = np.nan
                y = np.nan
            else:
                x = flip*vec[1]/vec[0]
                y = vec[2]/vec[0]
        else:
            x = np.zeros(vec[0].shape)+np.nan
            y = np.zeros(vec[0].shape)+np.nan
            x[w] = flip*vec[1][w]/vec[0][w]
            y[w] = vec[2][w]/vec[0][w]
        return x,y
    vec2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)

    def xy2vec(self, x, y=None, direct=False):
        flip = self._flip
        if y is None:
            x,y = x
        x,y=np.asarray(x),np.asarray(y)
        rm1=1./np.sqrt(1.+x**2+y**2)
        vec = (rm1,flip*rm1*x,rm1*y)
        if not direct:
            return self.rotator.I(vec)
        else:
            return vec
    xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name,name)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        vec=R.dir2vec(theta,phi,lonlat=lonlat)
        return self.vec2xy(vec,direct=direct)
    ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)
    
    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        return R.vec2dir(self.xy2vec(x,y,direct=direct),lonlat=lonlat)
    xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name,name)


    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        reso = self.arrayinfo['reso']
        if y is None: x,y = x
        dx = reso/60. * dtor
        xc,yc = 0.5*(xsize-1), 0.5*(ysize-1)
        j = np.around(xc+x/dx).astype(np.long)
        i = np.around(yc+y/dx).astype(np.long)
        return i,j
    xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name,name)

    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        reso = self.arrayinfo['reso']
        dx = reso/60. * dtor
        xc,yc = 0.5*(xsize-1), 0.5*(ysize-1)
        if i is None and j is None:
            idx=np.outer(np.ones(ysize),np.arange(xsize))
            x=(idx-xc) * dx   # astro= '-' sign, geo '+' sign
            idx=np.outer(np.arange(ysize),np.ones(xsize))
            y=(idx-yc)*dx #(idx-yc) * dx
        elif i is not None and j is not None:
            x=(np.asarray(j)-xc) * dx
            y=(np.asarray(i)-yc) * dx #(asarray(i)-yc) * dx
        elif i is not None and j is None:
            i, j = i
            x=(np.asarray(j)-xc) * dx
            y=(np.asarray(i)-yc) * dx #(i-yc) * dx
        else:
            raise TypeError("Wrong parameters")
        return x,y
    ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name,name)

    def get_extent(self):
        xsize,ysize = self.arrayinfo['xsize'],self.arrayinfo['ysize']
        left,bottom = self.ij2xy(0,0)
        right,top = self.ij2xy(ysize-1,xsize-1)
        return (left,right,bottom,top)

    def get_fov(self):
        vx,vy,vz = self.xy2vec(self.ij2xy(0,0), direct=True)
        a = np.arccos(vx)
        return 2.*a

class MollweideProj(SphericalProj):
    """This class provides class methods for Mollweide projection.
    """
    
    name = "Mollweide"
    __molldata = []

    def __init__(self, rot=None, coord=None, xsize=800, **kwds):
        self.__initialise_data()
        super(MollweideProj,self).__init__(rot=rot, coord=coord,
                                           xsize=xsize, **kwds)
        
    def set_proj_plane_info(self,xsize):
        super(MollweideProj,self).set_proj_plane_info(xsize=xsize)

    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct:
            theta,phi=R.vec2dir(self.rotator(vx,vy,vz))
        else:
            theta,phi=R.vec2dir(vx,vy,vz)
        flip = self._flip
        X,Y = MollweideProj.__molldata
        # set phi in [-pi,pi]
        phi = (phi+pi)%(2*pi)-pi
        lat = pi/2. - theta
        A = MollweideProj.__lininterp(X,Y,lat)
        x = flip*2./pi * phi * np.cos(A)
        y = np.sin(A)
        return x,y
    vec2xy.__doc__ = SphericalProj.vec2xy.__doc__ % (name,name)

    def xy2vec(self, x, y=None, direct=False):
        flip = self._flip
        if y is None: x,y = x
        mask = (np.asarray(x)**2/4.+np.asarray(y)**2 > 1.)
        w=np.where(mask == False)
        if not mask.any(): mask = np.ma.nomask
        if not hasattr(x,'__len__'):
            if mask is not np.ma.nomask:
                return np.nan,np.nan,np.nan
            else:
                s = np.sqrt((1-y)*(1+y))
                a = np.arcsin(y)
                z = 2./pi * (a + y*s)
                phi = flip * pi/2. * x/np.maximum(s,1.e-6)
                sz = np.sqrt((1-z)*(1+z))
                vec = sz*np.cos(phi),sz*np.sin(phi),z
                if not direct:
                    return self.rotator.I(vec)
                else:
                    return vec
        else:
            vec = (np.zeros(x.shape)+np.nan,
                   np.zeros(x.shape)+np.nan,
                   np.zeros(x.shape)+np.nan)
            s = np.sqrt((1-y[w])*(1+y[w]))
            a = np.arcsin(y[w])
            vec[2][w] = 2./pi * (a + y[w]*s)
            phi = flip * pi/2. * x[w]/np.maximum(s,1.e-6)
            sz = np.sqrt((1-vec[2][w])*(1+vec[2][w]))
            vec[0][w] = sz*np.cos(phi)
            vec[1][w] = sz*np.sin(phi)
            if not direct:
                return self.rotator.I(vec)
            else:
                return vec
    xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name,name)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        return self.vec2xy(R.dir2vec(theta,phi,lonlat=lonlat),direct=direct)
    ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)
    
    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        vec = self.xy2vec(x,y,direct=direct)
        return R.vec2dir(vec,lonlat=lonlat)
    xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name,name)


    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = xsize // 2
        if y is None: x,y = x
        xc,yc = (xsize-1.)/2., (ysize-1.)/2.
        if hasattr(x,'__len__'):
            j = np.around(x*xc/2.+xc).astype(np.long)
            i = np.around(yc+y*yc).astype(np.long)
            mask = (x**2/4.+y**2>1.)
            if not mask.any(): mask=np.ma.nomask
            j=np.ma.array(j,mask=mask)
            i=np.ma.array(i,mask=mask)
        else:
            if x**2/4.+y**2 > 1.:
                i,j=np.nan,np.nan
            else:
                j = np.around(x*xc/2.+xc).astype(np.long)
                i = np.around(yc+y*yc).astype(np.long)
        return i,j
    xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name,name)

    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = xsize // 2
        xc,yc=(xsize-1.)/2.,(ysize-1.)/2.
        if i is None and j is None:
            idx = np.outer(np.arange(ysize),np.ones(xsize))
            y = (idx-yc)/yc
            idx = np.outer(np.ones(ysize),np.arange(xsize))
            x = 2.*(idx-xc)/xc
            mask = x**2/4.+y**2 > 1.
            if not mask.any(): mask=np.ma.nomask
            x = np.ma.array(x,mask=mask)
            y = np.ma.array(y,mask=mask)
        elif i is not None and j is not None:
            y = (np.asarray(i)-yc)/yc
            x=2.*(np.asarray(j)-xc)/xc
            if x**2/4.+y**2 > 1.: x,y=np.nan,np.nan
        elif i is not None and j is None:
            i,j = i
            y=(np.asarray(i)-yc)/yc
            x=2.*(np.asarray(j)-xc)/xc
            if x**2/4.+y**2 > 1.: x,y=np.nan,np.nan
        else:
            raise TypeError("i and j must be both given or both not given")
        return x,y
    ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name,name)

    def get_extent(self):
        return (-2.0,2.0,-1.0,1.0)

    @staticmethod
    def __initialise_data():
        if len(MollweideProj.__molldata) == 0:
            X = (np.arange(1.,180.,1.)-90.)*dtor
            Y = MollweideProj.__findRoot(MollweideProj.__fmoll,
                                         MollweideProj.__dfmoll,
                                         X.copy(),X,niter=10)
            X = np.concatenate([[-pi/2],X,[pi/2]])
            Y = np.concatenate([[-pi/2],Y,[pi/2]])
            MollweideProj.__molldata.append( X )
            MollweideProj.__molldata.append( Y )
        return

    @staticmethod
    def __findRoot(f, df, x0, argsf=None, argsdf=None, niter=100):
        x = x0
        niter = min(abs(niter),1000)
        i = 0
        while i < niter:
            dx = -f(x,argsf)/df(x,argsdf)
            x += dx
            i += 1
        return x

    @staticmethod
    def __fmoll(x,args):
        return 2.*x+np.sin(2.*x)-pi*np.sin(args)

    @staticmethod
    def __dfmoll(x,args):
        return 2.*(1.+np.cos(2.*x))

    @staticmethod
    def __lininterp(X,Y,x):
        idx = X.searchsorted(x)
        y = Y[idx-1] + (Y[idx]-Y[idx-1])/(X[idx]-X[idx-1]) * (x-X[idx-1])
        return y


class CartesianProj(SphericalProj):
    """This class provides class methods for Cartesian projection.
    """
    
    name = "Cartesian"

    def __init__(self, rot=None, coord=None, xsize=800, ysize=None, lonra=None, 
                 latra=None, **kwds):
        super(CartesianProj,self).__init__(rot=rot, coord=coord,
                                           xsize=xsize, ysize=ysize, lonra=lonra, latra=latra, **kwds)
        
    def set_proj_plane_info(self,xsize,ysize,lonra,latra):
        if lonra is None: lonra = [-180.,180.]
        if latra is None: latra = [-90.,90.]
        if (len(lonra)!=2 or len(latra)!=2 or lonra[0]<-180. or lonra[1]>180.
            or latra[0]<-90 or latra[1]>90 or lonra[0]>=lonra[1] or latra[0]>=latra[1]):
            raise TypeError("Wrong argument lonra or latra. Must be lonra=[a,b],latra=[c,d] "
                            "with a<b, c<d, a>=-180, b<=180, c>=-90, d<=+90")
        lonra = self._flip*np.float64(lonra)[::self._flip]
        latra = np.float64(latra)
        xsize = np.long(xsize)
        if ysize is None:
            ratio = (latra[1]-latra[0])/(lonra[1]-lonra[0])
            ysize = np.long(round(ratio*xsize))
        else:
            ysize = np.long(ysize)
            ratio = float(ysize)/float(xsize)
        super(CartesianProj,self).set_proj_plane_info(xsize=xsize, lonra=lonra, latra=latra, 
                                                        ysize=ysize, ratio=ratio)

    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct:
            theta,phi=R.vec2dir(self.rotator(vx,vy,vz))
        else:
            theta,phi=R.vec2dir(vx,vy,vz)
        flip = self._flip
        # set phi in [-pi,pi]
        x = flip*((phi+pi)%(2*pi)-pi)
        x /= dtor # convert in degree
        y = pi/2. - theta
        y /= dtor # convert in degree
        return x,y
    vec2xy.__doc__ = SphericalProj.vec2xy.__doc__ % (name,name)

    def xy2vec(self, x, y=None, direct=False):
        if y is None:
            x,y = np.asarray(x)
        else:
            x,y = np.asarray(x),np.asarray(y)
        flip = self._flip
        theta = pi/2.-y*dtor # convert in radian
        phi = flip*x*dtor # convert in radian
        # dir2vec does not support 2d arrays, so first use flatten and then
        # reshape back to previous shape
        if not direct: 
            vec = self.rotator.I(R.dir2vec(theta.flatten(),phi.flatten()))
        else:
            vec = R.dir2vec(theta.flatten(),phi.flatten())
        vec = [v.reshape(theta.shape) for v in vec]
        return vec
    xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name,name)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        return self.vec2xy(R.dir2vec(theta,phi,lonlat=lonlat),direct=direct)
    ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)
    
    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        vec = self.xy2vec(x,y,direct=direct)
        return R.vec2dir(vec,lonlat=lonlat)
    xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name,name)


    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        lonra = self.arrayinfo['lonra']
        latra = self.arrayinfo['latra']
        if y is None: x,y = np.asarray(x)
        else: x,y = np.asarray(x), np.asarray(y)
        j = np.around((x-lonra[0])/(lonra[1]-lonra[0])*(xsize-1)).astype(np.int64)
        i = np.around((y-latra[0])/(latra[1]-latra[0])*(ysize-1)).astype(np.int64)
        if len(x.shape) > 0:
            mask = ((i<0)|(i>=ysize)|(j<0)|(j>=xsize))
            if not mask.any(): mask=np.ma.nomask
            j=np.ma.array(j,mask=mask)
            i=np.ma.array(i,mask=mask)
        else:
            if j<0 or j>=xsize or i<0 or i>=ysize: i=j=None
        return i,j
    xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name,name)

    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        lonra = self.arrayinfo['lonra']
        latra = self.arrayinfo['latra']
        if i is not None and j is None: i,j = np.asarray(i)
        elif i is not None and j is not None: i,j = np.asarray(i),np.asarray(j)
        if i is None and j is None:
            idx = np.outer(np.arange(ysize),np.ones(xsize))
            y = (float(latra[1]-latra[0])/(ysize-1.)) * idx
            y += latra[0]
            idx = np.outer(np.ones(ysize),np.arange(xsize))
            x = (float(lonra[1]-lonra[0])/(xsize-1.) * idx)
            x +=  lonra[0]
            x = np.ma.array(x)
            y = np.ma.array(y)
        elif i is not None and j is not None:
            y = (float(latra[1]-latra[0])/(ysize-1) ) * i 
            y += latra[0]
            x = (float(lonra[1]-lonra[0])/(xsize-1)) * j 
            x += lonra[0]
            if len(i.shape) > 0:
                mask = ((x<-180)|(x>180)|(y<-90)|(y>90))
                if not mask.any():
                    mask = np.ma.nomask
                x = np.ma.array(x,mask=mask)
                y = np.ma.array(y,mask=mask)
            else:
                if x<-180 or x>180 or y<-90 or y>90:
                    x = y = np.nan
        else:
            raise TypeError("i and j must be both given or both not given")
        return x,y
    ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name,name)

    def get_extent(self):
        lonra = self.arrayinfo['lonra']
        latra = self.arrayinfo['latra']
        return (lonra[0],lonra[1],latra[0],latra[1])
    get_extent.__doc__ = SphericalProj.get_extent.__doc__

    def get_fov(self):
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        v1 = np.asarray(self.xy2vec(self.ij2xy(0,0), direct=True))
        v2 = np.asarray(self.xy2vec(self.ij2xy(ysize-1,xsize-1), direct=True))
        a = np.arccos((v1*v2).sum())
        return 2*a

#    def get_fov(self):
#        lonra = self.arrayinfo['lonra']
#        latra = self.arrayinfo['latra']
#        return np.sqrt((lonra[1]-lonra[0])**2+(latra[1]-latra[0])**2)
        
    def get_center(self,lonlat=False):
        lonra = self.arrayinfo['lonra']
        latra = self.arrayinfo['latra']
        xc = 0.5*(lonra[1]+lonra[0])
        yc = 0.5*(latra[1]+latra[0])
        return self.xy2ang(xc,yc,lonlat=lonlat)
    get_center.__doc__ = SphericalProj.get_center.__doc__


class OrthographicProj(SphericalProj):
    """This class provides methods for orthographic projection
    """
    
    name = "Orthographic"
    
    def __init__(self, rot=None, coord=None, xsize=800, half_sky=False,**kwds):
        super(OrthographicProj,self).__init__(rot=rot, coord=coord,xsize=xsize,
                                              half_sky=half_sky,**kwds)
    
    def set_proj_plane_info(self,xsize,half_sky):
        super(OrthographicProj,self).set_proj_plane_info(xsize=xsize,
                                                         half_sky=half_sky)
    
    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct:
            theta,phi=R.vec2dir(self.rotator(vx,vy,vz))
        else:
            theta,phi=R.vec2dir(vx,vy,vz)
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        half_sky = self.arrayinfo['half_sky']
        flip = self._flip
        # set phi in [-pi,pi]
        phi = flip*(phi+pi)%(2*pi)-pi
        lat = pi/2. - theta
        x = np.cos(lat)*np.sin(phi)
        if not half_sky: x -= 1.0
        y = np.sin(lat)
        # unfold back of sphere
        cosc = np.cos(lat)*np.cos(phi)
        if np.any(cosc<0):
            hmask = (cosc<0)
            if hasattr(x,'__len__'):
                if half_sky:
                    x[hmask] = np.nan
                else:
                    x[hmask] *= -1
            elif hmask:
                if half_sky:
                    x = np.nan
                else:
                    x *= -1
        if half_sky:
            mask = (np.asarray(x)**2+np.asarray(y)**2>1.0)
        else:
            mask = ((np.mod(np.asarray(x)+2.0,2.0)-1.0)**2 + \
                    np.asarray(y)**2>1.0)
        if mask.any():
            if not hasattr(x,'__len__'):
                x = np.nan
                y = np.nan
            else:
                x[mask] = np.nan
                y[mask] = np.nan
        return x,y
    vec2xy.__doc__ = SphericalProj.vec2xy.__doc__ % (name,name)
    
    def xy2vec(self, x, y=None, direct=False):
        if y is None:
            x,y = x
        if hasattr(x,'__len__'):
            x,y = np.asarray(x),np.asarray(y)
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        half_sky = self.arrayinfo['half_sky']
        flip = self._flip
        # re-fold back of sphere
        mask = None
        if not half_sky:
            if hasattr(x,'__len__'):
                if np.any(x>0.0):
                    mask = (x>0.0)
                    x[mask] *= -1
            elif x>0:
                mask = 0
                x = -x
            x+=1.0
        r = np.sqrt(x**2+y**2)
        if hasattr(r,'__len__'):
            r[r>1] = np.nan
        elif r>1: r = np.nan
        c = np.arcsin(r)
        if hasattr(y,'__len__'):
            y[np.abs(y)>1] = np.nan
        elif np.abs(y)>1: y = np.nan
        lat = np.arcsin(y)
        phi = np.arctan2(x,np.cos(c))
        phi *= flip
        if not mask is None:
            if hasattr(phi,'__len__'):
                phi[mask] = pi-phi[mask]
            else: phi = pi-phi
        theta = pi/2. - lat
        vec = R.dir2vec(theta,phi)
        if not direct:
            return self.rotator.I(vec)
        else:
            return vec
    xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name,name)
    
    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        return self.vec2xy(R.dir2vec(theta,phi,lonlat=lonlat),direct=direct)
    ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)
    
    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        return R.vec2dir(self.xy2vec(x,y,direct=direct),lonlat=lonlat)
    xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name,name)

    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        xsize = int(self.arrayinfo['xsize'])
        half_sky = self.arrayinfo['half_sky']
        if half_sky: ratio = 1
        else: ratio = 2
        ysize = xsize // ratio
        if y is None: x,y = np.asarray(x)
        else: x,y = np.asarray(x), np.asarray(y)
        xc,yc = (xsize-1.)/2., (ysize-1.)/2.
        if hasattr(x,'__len__'):
            if half_sky:
                mask = (x**2+y**2>1.0)
            else:
                mask = ((np.mod(x+2.0,2.0)-1.0)**2+y**2>1.0)
            if not mask.any(): mask = np.ma.nomask
            j=np.ma.array(np.around(x*xc/ratio+xc).astype(np.long),mask=mask)
            i=np.ma.array(np.around(yc+y*yc).astype(np.long),mask=mask)
        else:
            if ( half_sky and x**2+y**2>1.0 ) or \
                   ( not half_sky and (np.mod(x+2.0,2.0)-1.0)**2+y**2>1.0 ):
                i,j,=np.nan,np.nan
            else:
                j = np.around(x*xc/ratio+xc).astype(np.long)
                i = np.around(yc+y*yc).astype(np.long)
        return i,j
    xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name,name)
    
    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        xsize = int(self.arrayinfo['xsize'])
        half_sky = self.arrayinfo['half_sky']
        if half_sky: ratio = 1
        else: ratio = 2
        ysize = xsize // ratio
        xc,yc=(xsize-1.)/2.,(ysize-1.)/2.
        if i is None and j is None:
            idx = np.outer(np.arange(ysize),np.ones(xsize))
            y = (idx-yc)/yc
            idx = np.outer(np.ones(ysize),np.arange(xsize))
            x = ratio*(idx-xc)/xc
        elif i is not None and j is not None:
            y = (np.asarray(i)-yc)/yc
            x = ratio*(np.asarray(j)-xc)/xc
            # if np.mod(x,1.0)**2+y**2 > 1.0: x,y=np.nan,np.nan
        elif i is not None and j is None:
            i,j = i
            y=(np.asarray(i)-yc)/yc
            x=ratio*(np.asarray(j)-xc)/xc
            # if np.mod(x,1.0)**2.+y**2 > 1.: x,y=np.nan,np.nan
        else:
            raise TypeError("i and j must be both given or both not given")
        if half_sky:
            mask = (x**2+y**2>1.)
        else:
            mask = ((np.mod(x+2.0,2.0)-1.0)**2+y**2 > 1.)
        if not mask.any(): mask=np.ma.nomask
        x = np.ma.array(x,mask=mask)
        y = np.ma.array(y,mask=mask)
        if len(x)==0: x = x[0]
        if len(y)==0: y = y[0]
        return x,y
    ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name,name)
    
    def get_extent(self):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        half_sky = self.arrayinfo['half_sky']
        if half_sky: ratio = 1.0
        else: ratio = 2.0
        return (-ratio,ratio,-1.0,1.0)
    get_extent.__doc__ = SphericalProj.get_extent.__doc__

class AzimuthalProj(SphericalProj):
    """This class provides methods for Lambert azimuthal equal-area projection and
    azimuthal equidistant projection
    """

    name = "Azimuthal"

    def __init__(self, rot=None, coord=None, xsize=None, ysize=None, reso=None, lamb=None, half_sky=None, **kwds):
        super(AzimuthalProj,self).__init__(rot=rot, coord=coord,xsize=xsize,ysize=ysize,reso=reso,lamb=lamb,half_sky=half_sky,**kwds)

    def set_proj_plane_info(self, xsize=800,ysize=None,reso=1.5,lamb=True,half_sky=False):
        if xsize is None: xsize=800
        if ysize is None: ysize=xsize
        if reso is None: reso=1.5
        if lamb is None: lamb=True
        if half_sky is None: half_sky=False
        super(AzimuthalProj,self).set_proj_plane_info(xsize=xsize,ysize=ysize,
                                                      reso=reso,lamb=lamb,half_sky=half_sky)

    def vec2xy(self, vx, vy=None, vz=None, direct=False):
        if not direct:
            theta,phi=R.vec2dir(self.rotator(vx,vy,vz))
        else:
            theta,phi=R.vec2dir(vx,vy,vz)
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        flip = self._flip
        lamb = self.arrayinfo['lamb']
        half_sky = self.arrayinfo['half_sky']
        # set phi in [-pi,pi]
        phi = flip*((phi+pi)%(2*pi)-pi)
        lat = pi/2. - theta
        if lamb:
            kprime = np.sqrt (2. / (1. + np.cos(lat) * np.cos(phi)))
        else:
            c = np.arccos(np.cos(lat) * np.cos(phi))
            kprime = c / np.sin(c)
        x = kprime * np.cos(lat) * np.sin(phi)
        y = kprime * np.sin(lat)
        if lamb: r2max = 4.
        else: r2max = pi**2
        if half_sky:
            if lamb: r2max /= 2.
            else: r2max /= 4.
        mask = (np.asarray(x)**2+np.asarray(y)**2 > r2max)
        if not hasattr(x,'__len__'):
            if mask is not np.ma.nomask:
                return np.nan,np.nan
        else:
            w = np.where(mask)
            x[w] = np.nan
            y[w] = np.nan 
        return x,y
    vec2xy.__doc__ = SphericalProj.vec2xy.__doc__ % (name,name)

    def xy2vec(self, x, y=None, direct=False):
        if y is None:
            x,y = x
        if hasattr(x,'__len__'):
            x,y = np.asarray(x),np.asarray(y)
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for"
                            " this projector")
        flip = self._flip
        lamb = self.arrayinfo['lamb']
        half_sky = self.arrayinfo['half_sky']
        if lamb: r2max = 4.
        else: r2max = pi**2
        if half_sky:
            if lamb: r2max /= 2.
            else: r2max /= 4.
        mask = (np.asarray(x)**2+np.asarray(y)**2 > r2max)
        w=np.where(mask == False)
        if not mask.any(): mask = np.ma.nomask
        if not hasattr(x,'__len__'):
            if mask is not np.ma.nomask:
                return np.nan,np.nan,np.nan
            else:
                rho = np.sqrt(x**2 + y**2)
                if lamb:
                  c = 2. * np.arcsin(rho/2.)
                else:
                  c = rho
                lat = np.arcsin(y * np.sin(c)/rho)
                phi = np.arctan2(x * np.sin(c), (rho * np.cos(c)))
                phi *= flip
                vec = R.dir2vec(pi/2.-lat,phi)
                if not direct:
                    return self.rotator.I(vec)
                else:
                    return vec
        else:
            vec = (np.zeros(x.shape)+np.nan,
                   np.zeros(x.shape)+np.nan,
                   np.zeros(x.shape)+np.nan)
            rho = np.sqrt(x[w]**2 + y[w]**2)
            if lamb:
              c = 2. * np.arcsin(rho/2.)
            else:
              c = rho
            lat = np.arcsin(y[w] * np.sin(c)/rho)
            phi = np.arctan2(x[w] * np.sin(c), (rho * np.cos(c)))
            phi *= flip
            vec[0][w] = np.cos(phi)*np.cos(lat)
            vec[1][w] = np.sin(phi)*np.cos(lat)
            vec[2][w] = np.sin(lat)
            if not direct:
                return self.rotator.I(vec)
            else:
                return vec
    xy2vec.__doc__ = SphericalProj.xy2vec.__doc__ % (name,name)

    def ang2xy(self, theta, phi=None, lonlat=False, direct=False):
        return self.vec2xy(R.dir2vec(theta,phi,lonlat=lonlat),direct=direct)
    ang2xy.__doc__ = SphericalProj.ang2xy.__doc__ % (name,name)

    def xy2ang(self, x, y=None, lonlat=False, direct=False):
        return R.vec2dir(self.xy2vec(x,y,direct=direct),lonlat=lonlat)
    xy2ang.__doc__ = SphericalProj.xy2ang.__doc__ % (name,name)

    def xy2ij(self, x, y=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        reso = self.arrayinfo['reso']
        lamb = self.arrayinfo['lamb']
        half_sky = self.arrayinfo['half_sky']
        if lamb: r2max = 4.
        else: r2max = pi**2
        if half_sky:
            if lamb: r2max /= 2.
            else: r2max /= 4.
        if y is None: x,y = x
        dx = reso/60. * dtor
        xc,yc = 0.5*(xsize-1), 0.5*(ysize-1)
        if hasattr(x,'__len__'):
            mask = (x**2+y**2>r2max)
            if not mask.any(): mask = np.ma.nomask
            j=np.ma.array(np.around(xc+x/dx).astype(np.long),mask=mask)
            i=np.ma.array(np.around(yc+y/dx).astype(np.long),mask=mask)
        else:
            if (x**2+y**2>r2max):
                i,j,=np.nan,np.nan
            else:
                j = np.around(xc+x/dx).astype(np.long)
                i = np.around(yc+y/dx).astype(np.long)
        return i,j
    xy2ij.__doc__ = SphericalProj.xy2ij.__doc__ % (name,name)

    def ij2xy(self, i=None, j=None):
        if self.arrayinfo is None:
            raise TypeError("No projection plane array information defined for "
                            "this projector")
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        reso = self.arrayinfo['reso']
        lamb = self.arrayinfo['lamb']
        half_sky = self.arrayinfo['half_sky']
        dx = reso/60. * dtor
        xc,yc = 0.5*(xsize-1), 0.5*(ysize-1)
        if lamb: r2max = 4.
        else: r2max = pi**2
        if half_sky:
            if lamb: r2max /= 2.
            else: r2max /= 4.
        if i is None and j is None:
            idx = np.outer(np.arange(ysize),np.ones(xsize))
            y = (idx-yc) * dx
            idx = np.outer(np.ones(ysize),np.arange(xsize))
            x = (idx-xc) * dx
        elif i is not None and j is not None:
            y = (np.asarray(i)-yc) * dx
            x = (np.asarray(j)-xc) * dx
        elif i is not None and j is None:
            i,j = i
            y=(np.asarray(i)-yc) * dx
            x=(np.asarray(j)-xc) * dx
        else:
            raise TypeError("i and j must be both given or both not given")
        if hasattr(x,'__len__'):
            mask = (x**2+y**2 > r2max)
            if not mask.any(): mask=np.ma.nomask
            x = np.ma.array(x,mask=mask)
            y = np.ma.array(y,mask=mask)
        else:
            if (x**2+y**2>r2max): x,y=np.nan,np.nan
        return x,y
    ij2xy.__doc__ = SphericalProj.ij2xy.__doc__ % (name,name)

    def get_extent(self):
        xsize = int(self.arrayinfo['xsize'])
        ysize = int(self.arrayinfo['ysize'])
        reso = self.arrayinfo['reso']
        dx = reso/60.0 * dtor
        xc,yc = 0.5*(xsize-1), 0.5*(ysize-1)
        left = -xc * dx
        bottom = -yc * dx
        right = (xsize-1-xc) * dx
        top = (ysize-1-yc) * dx
        return (left,right,bottom,top)
    get_extent.__doc__ = SphericalProj.get_extent.__doc__

    def get_fov(self):
        half_sky = self.arrayinfo['half_sky']
        vx,vy,vz = self.xy2vec(self.ij2xy(0,0), direct=True)
        a = np.arccos(vx)
        if np.isfinite(a):
          return 2.*a
        else:
          if half_sky: return pi
          else: return 2.*pi

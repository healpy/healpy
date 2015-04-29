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
import numpy as np
import six
import warnings

coordname = {'G': 'Galactic', 'E': 'Ecliptic', 'C': 'Equatorial'}

class ConsistencyWarning(Warning):
    """Warns for a problem in the consistency of data
    """
    pass

if __name__ != '__main__':
    warnings.filterwarnings("always", category=ConsistencyWarning, module=__name__)

class Rotator(object):
    """Rotation operator, including astronomical coordinate systems.

    This class provides tools for spherical rotations. It is meant to be used 
    in the healpy library for plotting, and for this reason reflects the
    convention used in the Healpix IDL library.

    Parameters
    ----------
    rot : None or sequence
      Describe the rotation by its euler angle. See :func:`euler_matrix_new`.
    coord : None or sequence of str
      Describe the coordinate system transform. If *rot* is also given, the
      coordinate transform is applied first, and then the rotation.
    inv : bool
      If True, the inverse rotation is defined. (Default: False)
    deg : bool
      If True, angles are assumed to be in degree. (Default: True)
    eulertype : str
      The Euler angle convention used. See :func:`euler_matrix_new`.

    Attributes
    ----------
    mat
    coordin
    coordout
    coordinstr
    coordoutstr
    rots
    coords
      
    Examples
    --------
    >>> r = Rotator(coord=['G','E'])  # Transforms galactic to ecliptic coordinates
    >>> theta_gal, phi_gal = np.pi/2., 0.
    >>> theta_ecl, phi_ecl = r(theta_gal, phi_gal)  # Apply the conversion
    >>> print(theta_ecl)
    1.66742286715
    >>> print(phi_ecl)
    -1.62596400306
    >>> theta_ecl, phi_ecl = Rotator(coord='ge')(theta_gal, phi_gal) # In one line
    >>> print(theta_ecl)
    1.66742286715
    >>> print(phi_ecl)
    -1.62596400306
    >>> vec_gal = np.array([1, 0, 0]) #Using vectors
    >>> vec_ecl = r(vec_gal)
    >>> print(vec_ecl)
    [-0.05488249 -0.99382103 -0.09647625]
    """
    ErrMessWrongPar = ("rot and coord must be single elements or "
                       "sequence of same size.")

    def __init__(self,rot=None,coord=None,inv=None,deg=True,
                 eulertype='ZYX'):
        """Create a rotator with given parameters.
        - rot: a float, a tuple of 1,2 or 3 floats or a sequence of tuples.
               If it is a sequence of tuple, it must have the same length as coord.
        - coord: a string or a tuple of 1 or 2 strings or a sequence of tuple
                 If it is a sequence of tuple, it must have same length as rot.
        - inv: whether to use inverse rotation or not
        - deg: if True, angles in rot are assumed in degree (default: True)
        - eulertype: the convention for Euler angles in rot.
        Note: the coord system conversion is applied first, then the rotation.
        """
        rot_is_seq = (hasattr(rot,'__len__') 
                      and hasattr(rot[0], '__len__'))
        coord_is_seq = (hasattr(coord,'__len__') 
                        and hasattr(coord[0],'__len__')
                        and type(coord[0]) is not str)
        if rot_is_seq and coord_is_seq:
            if len(rot) != len(coord):
                raise ValueError(Rotator.ErrMessWrongPar)
            else:
                rots = rot
                coords = coord
        elif (rot_is_seq or coord_is_seq) and (rot   is not None and
                                               coord is not None):
            raise ValueError(Rotator.ErrMessWrongPar)
        else:
            rots = [rot]
            coords = [coord]
        inv_is_seq = hasattr(inv,'__len__')
        if inv_is_seq:
            if len(inv) != len(rots):
                raise ValueError("inv must have same length as rot and/or coord")
            invs = inv
        else:
            invs = [inv]*len(rots)
        # check the argument and normalize them
        if eulertype in ['ZYX','X','Y']:
            self._eultype = eulertype
        else:
            self._eultype = 'ZYX'
        self._rots = []
        self._coords = []
        self._invs = []
        for r,c,i in zip(rots,coords,invs):
            rn = normalise_rot(r,deg=deg)
#            if self._eultype in ['X','Y']:
#                rn[1] = -rn[1]
            cn = normalise_coord(c)
            self._rots.append(rn)   # append(rn) or insert(0, rn) ?
            self._coords.append(cn) # append(cn) or insert(0, cn) ? 
            self._invs.append(bool(i))
        if not self.consistent:
            warnings.warn("The chain of coord system rotations is not consistent",
                          category=ConsistencyWarning)
        self._update_matrix()

    def _update_matrix(self):
        self._matrix = np.identity(3)
        self._do_rotation = False
        for r,c,i in zip(self._rots, self._coords,self._invs):
            rotmat,do_rot,rotnorm = get_rotation_matrix(r,
                                                        eulertype=self._eultype)
            convmat,do_conv,coordnorm = get_coordconv_matrix(c)
            r = np.dot(rotmat,convmat)
            if i: r = r.T
            self._matrix = np.dot(self._matrix, r)
            self._do_rotation = self._do_rotation or (do_rot or do_conv)

    def _is_coords_consistent(self):
        for c,i in zip(self._coords,self._invs):
            break
        for cnext,inext in zip(self._coords[1:],self._invs[1:]):
            if c[i] != cnext[not inext]:
                return False
            c,i = cnext,inext
        return True
    consistent = property(_is_coords_consistent,
                          doc="consistency of the coords transform chain")

    def __eq__(self,a):
        if type(a) is not type(self): return False
        # compare the _rots
        v = [np.allclose(x,y,rtol=0,atol=1e-15) for x,y in zip(self._rots,a._rots)]
        return ( np.array(v).all() and
                 (self._coords == a._coords) and
                 (self._invs == a._invs) )
    
    def __call__(self,*args,**kwds):
        """Use the rotator to rotate either spherical coordinates (theta, phi)
        or a vector (x,y,z). You can use lonla keyword to use longitude, latitude
        (in degree) instead of theta, phi (in radian). In this case, returns 
        longitude, latitude in degree.

        Accepted forms:

        r(x,y,z)  # x,y,z either scalars or arrays
        r(theta,phi) # theta, phi scalars or arrays 
        r(lon,lat,lonlat=True)  # lon, lat scalars or arrays
        r(vec) # vec 1-D array with 3 elements, or 2-D array 3xN
        r(direction) # direction 1-D array with 2 elements, or 2xN array

        Parameters
        ----------
        vec_or_dir : array or multiple arrays
          The direction to rotate. See above for accepted formats.
        lonlat : bool, optional
          If True, assumes the input direction is longitude/latitude in degrees.
          Otherwise, assumes co-latitude/longitude in radians. Default: False
        inv : bool, optional
          If True, applies the inverse rotation. Default: False.
        """
        if kwds.pop('inv',False): m=self._matrix.T
        else:                     m=self._matrix
        lonlat = kwds.pop('lonlat',False)
        if len(args) == 1:
            arg=args[0]
            if not hasattr(arg,'__len__') or len(arg) < 2 or len(arg) > 3:
                raise TypeError('Argument must be a sequence of 2 or 3 '
                                'elements')
            if len(arg) == 2:
                return rotateDirection(m,arg[0],arg[1],
                                       self._do_rotation,lonlat=lonlat)
            else:
                return rotateVector(m,arg[0],arg[1],arg[2],
                                    self._do_rotation)
        elif len(args) == 2:
            return rotateDirection(m,args[0],args[1],
                                   self._do_rotation,lonlat=lonlat)
        elif len(args) == 3:
            return rotateVector(m,args[0],args[1],args[2],
                                self._do_rotation)
        else:
            raise TypeError('Either 1, 2 or 3 arguments accepted')

    def __mul__(self,a):
        """Composition of rotation.
        """
        if not isinstance(a,Rotator):
            raise TypeError("A Rotator can only multiply another Rotator "
                            "(composition of rotations)")
        rots = self._rots + a._rots
        coords = self._coords + a._coords
        invs = self._invs + a._invs
        return Rotator(rot=rots,coord=coords,inv=invs,deg=False)
    
    def __rmul__(self,b):
        if not isinstance(b,Rotator):
            raise TypeError("A Rotator can only be multiplied by another Rotator "
                            "(composition of rotations)")
        rots = b._rots + self._rots
        coords = b._coords + self._coords
        invs = self._invs + a._invs
        return Rotator(rot=rots,coord=coords,inv=invs,deg=False)

    def __nonzero__(self):
        return self._do_rotation

    def get_inverse(self):
        rots = self._rots[::-1]
        coords = self._coords[::-1]
        invs = [ not i for i in self._invs[::-1]]
        return Rotator(rot=rots,coord=coords,inv=invs,deg=False)
    #I = property(get_inverse,doc='Return a new rotator representing the '
    #             'inverse rotation')

    def I(self,*args,**kwds):
        """Rotate the given vector or direction using the inverse matrix.
        rot.I(vec) <==> rot(vec,inv=True)
        """
        kwds['inv'] = True
        return self.__call__(*args,**kwds)
    
    @property
    def mat(self):
        """The matrix representing the rotation.
        """
        return np.matrix(self._matrix)

    @property
    def coordin(self):
        """The input coordinate system.
        """
        if not self.consistent: return None
        for c,i in zip(self._coords,self._invs):
            pass
        return c[i]

    @property
    def coordout(self):
        """The output coordinate system.
        """
        if not self.consistent: return None
        for c,i in zip(self._coords,self._invs):
            pass
        return c[not i]

    @property
    def coordinstr(self):
        """The input coordinate system in str.
        """
        return coordname.get(self.coordin,'')

    @property
    def coordoutstr(self):
        """The output coordinate system in str.
        """
        return coordname.get(self.coordout,'')

    @property
    def rots(self):
        """The sequence of rots defining the rotation.
        """
        return self._rots
    
    @property
    def coords(self):
        """The sequence of coords defining the rotation.
        """
        return self._coords
    
    def do_rot(self,i):
        """Returns True if rotation is not (close to) identity.
        """
        return not np.allclose(self.rots[i],np.zeros(3),rtol=0.,atol=1.e-15)

    def angle_ref(self,*args,**kwds):
        """Compute the angle between transverse reference direction of initial and final frames

        For example, if angle of polarisation is psi in initial frame, it will be psi+angle_ref in final
        frame.

        Parameters
        ----------
        dir_or_vec : array
          Direction or vector (see Rotator.__call__)
        lonlat: bool, optional
          If True, assume input is longitude,latitude in degrees. Otherwise,
          theta,phi in radian. Default: False
        inv : bool, optional
          If True, use the inverse transforms. Default: False

        Returns
        -------
        angle : float, scalar or array
          Angle in radian (a scalar or an array if input is a sequence of direction/vector)
        """
        R = self
        lonlat = kwds.get('lonlat',False)
        inv = kwds.get('inv',False)
        if len(args) == 1:
            arg=args[0]
            if not hasattr(arg,'__len__') or len(arg) < 2 or len(arg) > 3:
                raise TypeError('Argument must be a sequence of 2 or 3 '
                                'elements')
            if len(arg) == 2:
                v = dir2vec(arg[0],arg[1],lonlat=lonlat)
            else:
                v = arg
        elif len(args) == 2:
            v = dir2vec(args[0],args[1],lonlat=lonlat)
        elif len(args) == 3:
            v = args
        else:
            raise TypeError('Either 1, 2 or 3 arguments accepted')
        vp = R(v,inv=inv)
        north_pole = R([0.,0.,1.],inv=inv)
        sinalpha = north_pole[0]*vp[1]-north_pole[1]*vp[0]
        cosalpha = north_pole[2] - vp[2]*np.dot(north_pole,vp)
        return np.arctan2(sinalpha,cosalpha)

    def __repr__(self):
        return '[ '+', '.join([str(self._coords),
                               str(self._rots),
                               str(self._invs)]) +' ]'
    __str__ = __repr__



################################################################
#
#     Helpers function for rotation
#     used in the Rotator class.

def rotateVector(rotmat,vec,vy=None,vz=None, do_rot=True):
    """Rotate a vector (or a list of vectors) using the rotation matrix
    given as first argument.
    
    Parameters
    ----------
    rotmat : float, array-like shape (3,3)
      The rotation matrix
    vec : float, scalar or array-like
      The vector to transform (shape (3,) or (3,N)),
      or x component (scalar or shape (N,)) if vy and vz are given
    vy : float, scalar or array-like, optional
      The y component of the vector (scalar or shape (N,))
    vz : float, scalar or array-like, optional
      The z component of the vector (scalar or shape (N,))
    do_rot : bool, optional
      if True, really perform the operation, if False do nothing.

    Returns
    -------
    vec : float, array
      The component of the rotated vector(s).

    See Also
    --------
    Rotator
    """
    if vy is None and vz is None:
       if do_rot: return np.tensordot(rotmat,vec,axes=(1,0))
       else: return vec
    elif vy is not None and vz is not None:
       if do_rot: return np.tensordot(rotmat,np.array([vec,vy,vz]),axes=(1,0))
       else: return vec,vy,vz
    else:
       raise TypeError("You must give either vec only or vec, vy "
                       "and vz parameters")

def rotateDirection(rotmat,theta,phi=None,do_rot=True,lonlat=False):
    """Rotate the vector described by angles theta,phi using the rotation matrix
    given as first argument.
   
    Parameters
    ----------
    rotmat : float, array-like shape (3,3)
      The rotation matrix
    theta : float, scalar or array-like
      The angle theta (scalar or shape (N,)) 
      or both angles (scalar or shape (2, N)) if phi is not given.
    phi : float, scalar or array-like, optionnal
      The angle phi (scalar or shape (N,)).
    do_rot : bool, optional
      if True, really perform the operation, if False do nothing.
    lonlat : bool
      If True, input angles are assumed to be longitude and latitude in degree,
      otherwise, they are co-latitude and longitude in radians.
  
    Returns
    -------
    angles : float, array
      The angles of describing the rotated vector(s).

    See Also
    --------
    Rotator
    """
    vx,vy,vz=rotateVector(rotmat,dir2vec(theta,phi,lonlat=lonlat),do_rot=do_rot)
    return vec2dir(vx,vy,vz,lonlat=lonlat)

def vec2dir(vec,vy=None,vz=None,lonlat=False):
    """Transform a vector to angle given by theta,phi.
    
    Parameters
    ----------
    vec : float, scalar or array-like
      The vector to transform (shape (3,) or (3,N)),
      or x component (scalar or shape (N,)) if vy and vz are given
    vy : float, scalar or array-like, optional
      The y component of the vector (scalar or shape (N,))
    vz : float, scalar or array-like, optional
      The z component of the vector (scalar or shape (N,))
    lonlat : bool, optional
      If True, return angles will be longitude and latitude in degree,
      otherwise, angles will be longitude and co-latitude in radians (default)

    Returns
    -------
    angles : float, array
      The angles (unit depending on *lonlat*) in an array of 
      shape (2,) (if scalar input) or (2, N)

    See Also
    --------
    :func:`dir2vec`, :func:`pixelfunc.ang2vec`, :func:`pixelfunc.vec2ang`
    """
    if np.any(np.isnan(vec)):
        return np.nan, np.nan
    if vy is None and vz is None:
        vx,vy,vz = vec
    elif vy is not None and vz is not None:
        vx=vec
    else:
        raise TypeError("You must either give both vy and vz or none of them")
    r = np.sqrt(vx**2+vy**2+vz**2)
    ang = np.empty((2, r.size))
    ang[0, :] = np.arccos(vz / r)
    ang[1, :] = np.arctan2(vy, vx)
    if lonlat:
        ang = np.degrees(ang)
        np.negative(ang[0, :], ang[0, :])
        ang[0, :] += 90.
        return ang[::-1,:].squeeze()
    else:
        return ang.squeeze()

def dir2vec(theta,phi=None,lonlat=False):
    """Transform a direction theta,phi to a unit vector.
    
    Parameters
    ----------
    theta : float, scalar or array-like
      The angle theta (scalar or shape (N,)) 
      or both angles (scalar or shape (2, N)) if phi is not given.
    phi : float, scalar or array-like, optionnal
      The angle phi (scalar or shape (N,)).
    lonlat : bool
      If True, input angles are assumed to be longitude and latitude in degree,
      otherwise, they are co-latitude and longitude in radians.
    
    Returns
    -------
    vec : array
      The vector(s) corresponding to given angles, shape is (3,) or (3, N).
 
    See Also
    --------
    :func:`vec2dir`, :func:`pixelfunc.ang2vec`, :func:`pixelfunc.vec2ang`
    """
    if phi is None:
        theta,phi=theta
    if lonlat:
        lon,lat=theta,phi
        theta,phi = np.pi/2.-np.radians(lat),np.radians(lon)
    ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
    vec = np.empty((3, ct.size), np.float64)
    vec[0, :] = st * cp
    vec[1, :] = st * sp
    vec[2, :] = ct
    return vec.squeeze()

def angdist(dir1,dir2,lonlat=False):
    """Returns the angular distance between dir1 and dir2.

    Parameters
    ----------
    dir1, dir2 : float, array-like
      The directions between which computing the angular distance.
      Angular if len(dir) == 2 or vector if len(dir) == 3.
      See *lonlat* for unit
    lonlat : bool, scalar or sequence
      If True, angles are assumed to be longitude and latitude in degree,
      otherwise they are interpreted as colatitude and longitude in radian.
      If a sequence, lonlat[0] applies to dir1 and lonlat[1] applies to dir2.

    Returns
    -------
    angles : float, scalar or array-like
      The angle(s) between dir1 and dir2 in radian.

    Examples
    --------
    >>> import healpy as hp
    >>> hp.rotator.angdist([.2,0], [.2, 1e-6])
    array([  1.98669331e-07])
    """
    if hasattr(lonlat,'__len__') and len(lonlat) == 2:
        lonlat1,lonlat2 = lonlat
    else:
        lonlat1=lonlat2=lonlat
    dir1 = np.asarray(dir1)
    dir2 = np.asarray(dir2)
    if dir1.ndim == 2:
        if dir1.shape[0] == 2: # theta, phi -> vec
            vec1 = dir2vec(dir1, lonlat = lonlat1)
        else:
            vec1 = np.reshape(dir1, (3, -1))
            vec1 = normalize_vec(vec1)
    elif dir1.ndim == 1:
        if dir1.shape[0] == 2: # theta, phi -> vec
            vec1 = np.reshape(dir2vec(dir1, lonlat = lonlat1), (3, 1))
        else:
            vec1 = np.reshape(dir1, (3, 1))
            vec1 = normalize_vec(vec1)
    if dir2.ndim == 2:
        if dir2.shape[0] == 2: # theta, phi -> vec
            vec2 = dir2vec(dir2, lonlat = lonlat2)
        else:
            vec2 = np.reshape(dir2, (3, -1))
            vec2 = normalize_vec(vec2)
    elif dir2.ndim == 1:
        if dir2.shape[0] == 2: # theta, phi -> vec
            vec2 = np.reshape(dir2vec(dir2, lonlat = lonlat2), (3, 1))
        else:
            vec2 = np.reshape(dir2, (3, 1))
            vec2 = normalize_vec(vec2)
    # compute vec product
    vec_prod = np.sqrt((np.cross(vec1.T, vec2.T)**2).sum(axis=1))
    # compute scalar product
    scal_prod = (vec1*vec2).sum(axis=0)

    return np.arctan2(vec_prod, scal_prod)

def normalize_vec(vec):
    """Normalize the vector(s) *vec* (in-place if it is a ndarray).

    Parameters
    ----------
    vec : float, array-like of shape (D,) or (D, N)
      The D-vector(s) to normalize.

    Returns
    -------
    vec_normed : float, array
      Normalized vec, shape (D,) or (D, N)
    """
    vec = np.array(vec, np.float64)
    r = np.sqrt(np.sum(vec ** 2, axis = 0))
    vec /= r
    return vec

#######################################################
#
#   Manage the coord system conventions
#

def check_coord(c):
    """Check if parameter is a valid coord system.
    Raise a TypeError exception if it is not, otherwise returns the normalized
    coordinate system name.
    """
    if c is None:
        return c
    if not isinstance(c, six.string_types):
        raise TypeError('Coordinate must be a string (G[alactic],'
                        ' E[cliptic], C[elestial]'
                        ' or Equatorial=Celestial)')
    if c[0].upper() == 'G':
        x='G'
    elif c[0].upper() == 'E' and c != 'Equatorial':
        x='E'
    elif c[0].upper() == 'C' or c == 'Equatorial':
        x='C'
    else:
        raise ValueError('Wrong coordinate (either G[alactic],'
                         ' E[cliptic], C[elestial]'
                         ' or Equatorial=Celestial)')
    return x

def normalise_coord(coord):
    """Normalise the coord argument.
    Coord sys are either 'E','G', 'C' or 'X' if undefined.

    Input: either a string or a sequence of string.
           
    Output: a tuple of two strings, each being one of the norm coord sys name
            above.

    eg, 'E' -> ['E','E'], ['Ecliptic','G'] -> ['E','G']
    None -> ['X','X'] etc.
    """
    coord_norm = []
    if coord is None:
        coord = (None,None)
    coord=tuple(coord)
    if len(coord) > 2:
        raise TypeError('Coordinate must be a string (G[alactic],'
                        ' E[cliptic] or C[elestial])'
                        ' or a sequence of 2 strings')
    for x in coord:
        coord_norm.append(check_coord(x))
    if len(coord_norm) < 2:
        coord_norm.append(coord_norm[0])
    return tuple(coord_norm)

def normalise_rot(rot,deg=False):
   """Return rot possibly completed with zeroes to reach size 3.
   If rot is None, return a vector of 0.
   If deg is True, convert from degree to radian, otherwise assume input
   is in radian.
   """
   if deg: convert=np.pi/180.
   else: convert=1.
   if rot is None:
      rot=np.zeros(3)
   else:
      rot=np.array(rot,np.float64).flatten()*convert
      rot.resize(3)
   return rot

def get_rotation_matrix(rot, deg=False, eulertype='ZYX'):
   """Return the rotation matrix corresponding to angles given in rot.
   
   Usage: matrot,do_rot,normrot = get_rotation_matrix(rot)
   
   Input:
      - rot: either None, an angle or a tuple of 1,2 or 3 angles
             corresponding to Euler angles.
   Output:
      - matrot: 3x3 rotation matrix
      - do_rot: True if rotation is not identity, False otherwise
      - normrot: the normalized version of the input rot.
   """
   rot = normalise_rot(rot, deg=deg)
   if not np.allclose(rot,np.zeros(3),rtol=0.,atol=1.e-15):
      do_rot = True
   else:
      do_rot = False
   if eulertype == 'X':
       matrot=euler_matrix_new(rot[0],-rot[1],rot[2],X=True)
   elif eulertype == 'Y':
       matrot=euler_matrix_new(rot[0],-rot[1],rot[2],Y=True)
   else:
       matrot=euler_matrix_new(rot[0],-rot[1],rot[2],ZYX=True)
       
   return matrot,do_rot,rot
    
def get_coordconv_matrix(coord):
   """Return the rotation matrix corresponding to coord systems given
   in coord.

   Usage: matconv,do_conv,normcoord = get_coordconv_matrix(coord)

   Input:
      - coord: a tuple with initial and final coord systems.
               See normalise_coord.
   Output:
      - matconv: the euler matrix for coord sys conversion
      - do_conv: True if matconv is not identity, False otherwise
      - normcoord: the tuple of initial and final coord sys.

   History: adapted from CGIS IDL library.
   """
   
   coord_norm = normalise_coord(coord)
   
   if coord_norm[0] == coord_norm[1]:
      matconv = np.identity(3)
      do_conv = False        
   else:
      eps = 23.452294 - 0.0130125 - 1.63889E-6 + 5.02778E-7
      eps = eps * np.pi / 180.
      
      # ecliptic to galactic
      e2g = np.array([[-0.054882486, -0.993821033, -0.096476249],
                   [ 0.494116468, -0.110993846,  0.862281440],
                   [-0.867661702, -0.000346354,  0.497154957]])
      
      # ecliptic to equatorial
      e2q = np.array([[1.,     0.    ,      0.         ],
                   [0., np.cos( eps ), -1. * np.sin( eps )],
                   [0., np.sin( eps ),    np.cos( eps )   ]])
      
      # galactic to ecliptic
      g2e = np.linalg.inv(e2g)
      
      # galactic to equatorial                   
      g2q = np.dot(e2q , g2e)
      
      # equatorial to ecliptic
      q2e = np.linalg.inv(e2q)
      
      # equatorial to galactic
      q2g = np.dot(e2g , q2e)
   
      if coord_norm == ('E','G'):
         matconv = e2g
      elif coord_norm == ('G','E'):
         matconv = g2e
      elif coord_norm == ('E','C'):
         matconv = e2q
      elif coord_norm == ('C','E'):
         matconv = q2e
      elif coord_norm == ('C','G'):
         matconv = q2g
      elif coord_norm == ('G','C'):
         matconv = g2q
      else:
         raise ValueError('Wrong coord transform :',coord_norm)
      do_conv = True
      
   return matconv,do_conv,coord_norm


###################################################
##                                               ##
##                euler functions                ## 
##                                               ##
######                                      #######

def euler(ai, bi, select, FK4 = 0):
   """
   NAME:
       euler
   PURPOSE:
       Transform between Galactic, celestial, and ecliptic coordinates.
   EXPLANATION:
       Use the procedure ASTRO to use this routine interactively
   
   CALLING SEQUENCE:
        EULER, AI, BI, AO, BO, [ SELECT, /FK4, SELECT = ] 
   
   INPUTS:
         AI - Input Longitude in DEGREES, scalar or vector.  If only two 
                 parameters are supplied, then  AI and BI will be modified
                 to contain the output longitude and latitude.
         BI - Input Latitude in DEGREES
   
   OPTIONAL INPUT:
         SELECT - Integer (1-6) specifying type of coordinate
                  transformation.
   
        SELECT   From          To        |   SELECT      From         To
         1     RA-Dec (2000)  Galactic   |     4       Ecliptic     RA-Dec
         2     Galactic       RA-DEC     |     5       Ecliptic    Galactic
         3     RA-Dec         Ecliptic   |     6       Galactic    Ecliptic
   
        If not supplied as a parameter or keyword, then EULER will prompt
        for the value of SELECT
        Celestial coordinates (RA, Dec) should be given in equinox J2000 
        unless the /FK4 keyword is set.
   OUTPUTS:
         AO - Output Longitude in DEGREES
         BO - Output Latitude in DEGREES
   
   INPUT KEYWORD:
         /FK4 - If this keyword is set and non-zero, then input and output 
               celestial and ecliptic coordinates should be given in
               equinox B1950.
         /SELECT  - The coordinate conversion integer (1-6) may
                    alternatively be specified as a keyword
   NOTES:
         EULER was changed in December 1998 to use J2000 coordinates as the
         default, ** and may be incompatible with earlier versions***.
   REVISION HISTORY:
         Written W. Landsman,  February 1987
         Adapted from Fortran by Daryl Yentis NRL
         Converted to IDL V5.0   W. Landsman   September 1997
         Made J2000 the default, added /FK4 keyword
          W. Landsman December 1998
         Add option to specify SELECT as a keyword W. Landsman March 2003
         Converted to python by K. Ganga December 2007
   """

   # npar = N_params()
   #  if npar LT 2 then begin
   #     print,'Syntax - EULER, AI, BI, A0, B0, [ SELECT, /FK4, SELECT= ]'
   #     print,'    AI,BI - Input longitude,latitude in degrees'
   #     print,'    AO,BO - Output longitude, latitude in degrees'
   #     print,'    SELECT - Scalar (1-6) specifying transformation type'
   #     return
   #  endif

   PI = np.pi
   twopi   =   2.0*PI
   fourpi  =   4.0*PI
   deg_to_rad = 180.0/PI
   # 
   # ;   J2000 coordinate conversions are based on the following constants
   # ;   (see the Hipparcos explanatory supplement).
   # ;  eps = 23.4392911111 # Obliquity of the ecliptic
   # ;  alphaG = 192.85948d           Right Ascension of Galactic North Pole
   # ;  deltaG = 27.12825d            Declination of Galactic North Pole
   # ;  lomega = 32.93192d            Galactic longitude of celestial equator  
   # ;  alphaE = 180.02322d           Ecliptic longitude of Galactic North Pole
   # ;  deltaE = 29.811438523d        Ecliptic latitude of Galactic North Pole
   # ;  Eomega  = 6.3839743d          Galactic longitude of ecliptic equator
   # 
   if FK4 == 1:
      
      equinox = '(B1950)' 
      psi   = [ 0.57595865315, 4.9261918136,
                0.00000000000, 0.0000000000,
                0.11129056012, 4.7005372834]     
      stheta =[ 0.88781538514,-0.88781538514,
                0.39788119938,-0.39788119938,
                0.86766174755,-0.86766174755]    
      ctheta =[ 0.46019978478, 0.46019978478,
                0.91743694670, 0.91743694670,
                0.49715499774, 0.49715499774]    
      phi  = [ 4.9261918136,  0.57595865315,
               0.0000000000, 0.00000000000,
               4.7005372834, 0.11129056012]
   else:
      
      equinox = '(J2000)'
      psi   = [ 0.57477043300, 4.9368292465,  
                0.00000000000, 0.0000000000,  
                0.11142137093, 4.71279419371]     
      stheta =[ 0.88998808748,-0.88998808748, 
                0.39777715593,-0.39777715593, 
                0.86766622025,-0.86766622025]    
      ctheta =[ 0.45598377618, 0.45598377618, 
                0.91748206207, 0.91748206207, 
                0.49714719172, 0.49714719172]    
      phi  = [ 4.9368292465,  0.57477043300, 
               0.0000000000, 0.00000000000, 
               4.71279419371, 0.11142137093]
   # 
   i  = select - 1                         # IDL offset
   a  = ai/deg_to_rad - phi[i]
   b = bi/deg_to_rad
   sb = np.sin(b)
   cb = np.cos(b)
   cbsa = cb * np.sin(a)
   b  = -stheta[i] * cbsa + ctheta[i] * sb
   #bo    = math.asin(where(b<1.0, b, 1.0)*deg_to_rad)
   bo    = np.arcsin(b)*deg_to_rad
   #
   a = np.arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * np.cos(a) )
   ao = np.fmod( (a+psi[i]+fourpi), twopi) * deg_to_rad
   return ao, bo


def euler_matrix_new(a1,a2,a3,X=True,Y=False,ZYX=False,deg=False):
   """
   NAME:
     euler_matrix_new

   PURPOSE:
     computes the Euler matrix of an arbitrary rotation described
     by 3 Euler angles
     correct bugs present in Euler_Matrix
       
   CALLING SEQUENCE:
     result = euler_matrix_new (a1, a2, a3 [,X, Y, ZYX, DEG ])

   INPUTS:
     a1, a2, a3 = Euler angles, scalar
                  (in radian by default, in degree if DEG is set)
                  all the angles are measured counterclockwise
                  correspond to x, y, zyx-conventions (see Goldstein)
                  the default is x

   KEYWORD PARAMETERS:
      DEG : if set the angle are measured in degree

      X : 	rotation a1 around original Z 
    		rotation a2 around interm   X 
    		rotation a3 around final    Z
                DEFAULT,  classical mechanics convention

      Y : 	rotation a1 around original Z
    		rotation a2 around interm   Y
    		rotation a3 around final    Z
                quantum mechanics convention (override X)

      ZYX : 	rotation a1 around original Z
    		rotation a2 around interm   Y
    		rotation a3 around final    X
                aeronautics convention (override X)
                
      * these last three keywords are obviously mutually exclusive *

   OUTPUTS:
      result is a 3x3 matrix

   USAGE:
     if vec is an Nx3 array containing N 3D vectors,
     vec # euler_matrix_new(a1,a2,a3,/Y) will be the rotated vectors


   MODIFICATION HISTORY:
      March 2002, EH, Caltech, rewritting of euler_matrix

      convention   euler_matrix_new           euler_matrix
     X:       M_new(a,b,c,/X)  =  M_old(-a,-b,-c,/X) = Transpose( M_old(c, b, a,/X))
     Y:       M_new(a,b,c,/Y)  =  M_old(-a, b,-c,/Y) = Transpose( M_old(c,-b, a,/Y))
   ZYX:       M_new(a,b,c,/Z)  =  M_old(-a, b,-c,/Z)
   """
   
   t_k = 0
   if ZYX: t_k = t_k + 1
   #if X:   t_k = t_k + 1
   if Y:   t_k = t_k + 1
   if t_k > 1:
      raise ValueError('Choose either X, Y or ZYX convention')
   
   convert = 1.0
   if deg:
      convert = np.pi/180.
      
   c1 = np.cos(a1*convert)
   s1 = np.sin(a1*convert)
   c2 = np.cos(a2*convert)
   s2 = np.sin(a2*convert)
   c3 = np.cos(a3*convert)
   s3 = np.sin(a3*convert)
        
   if ZYX:
      m1 = np.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = np.array([[ c2,  0, s2],
                  [  0,  1,  0],
                  [-s2,  0, c2]]) # around   y

      m3 = np.array([[  1,  0,  0],
                  [  0, c3,-s3],
                  [  0, s3, c3]]) # around   x

   elif Y:
      m1 = np.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = np.array([[ c2,  0, s2],
                  [  0,  1,  0],
                  [-s2,  0, c2]]) # around   y

      m3 = np.array([[ c3,-s3,  0],
                  [ s3, c3,  0],
                  [  0,  0,  1]]) # around   z

   else:
      m1 = np.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = np.array([[  1,  0,  0],
                  [  0, c2,-s2],
                  [  0, s2, c2]]) # around   x

      m3 = np.array([[ c3,-s3,  0],
                  [ s3, c3,  0],
                  [  0,  0,  1]]) # around   z

   M = np.dot(m3.T,np.dot(m2.T,m1.T)) 

   return M


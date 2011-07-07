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
import numpy as npy
import warnings

coordname = {'G': 'Galactic', 'E': 'Ecliptic', 'C': 'Equatorial'}

class ConsistencyWarning(Warning):
    """Warns for a problem in the consistency of data
    """
    pass

if __name__ != '__main__':
    warnings.filterwarnings("always", category=ConsistencyWarning, module=__name__)

class Rotator:
    """This class provides tools for spherical rotations. It is meant to be used 
    in the healpy library for plotting, and for this reason reflects the
    convention used in the Healpix IDL library.

    Example:
      >>> r = Rotator(coord=['G','E'])
      >>> theta_ecl, phi_ecl = r(theta_gal, phi_gal)
    or in a shorter way:
      >>> theta_ecl, phi_ecl = Rotator(coord=['G','E'])(theta_gal, phi_gal)
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
        self._matrix = npy.identity(3)
        self._do_rotation = False
        for r,c,i in zip(self._rots, self._coords,self._invs):
            rotmat,do_rot,rotnorm = get_rotation_matrix(r,
                                                        eulertype=self._eultype)
            convmat,do_conv,coordnorm = get_coordconv_matrix(c)
            r = npy.dot(rotmat,convmat)
            if i: r = r.T
            self._matrix = npy.dot(self._matrix, r)
            self._do_rotation = self._do_rotation or (do_rot or do_conv)

    def _is_coords_consistent(self):
        c,i = zip(self._coords,self._invs)[0]
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
        v = [npy.allclose(x,y,rtol=0,atol=1e-15) for x,y in zip(self._rots,a._rots)]
        return ( npy.array(v).all() and
                 (self._coords == a._coords) and
                 (self._invs == a._invs) )
    
    def __call__(self,*args,**kwds):
        """Use the rotator to rotate either spherical coordinates (theta, phi)
        or a vector (x,y,z). You can use lonla keyword to use longitude, latitude
        (in degree) instead of theta, phi (in radian). In this case, returns 
        longitude, latitude in degree.
        Accepted forms:
        >>> r = Rotator()
        >>> r(x,y,z)  # x,y,z either scalars or arrays
        >>> r(theta,phi) # theta, phi scalars or arrays 
        >>> r(lon,lat,lonlat=True)  # lon, lat scalars or arrays
        >>> r(vec) # vec 1-D array with 3 elements, or 2-D array 3xN
        >>> r(direction) # direction 1-D array with 2 elements, or 2xN array
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
    
    def get_matrix(self):
        return npy.matrix(self._matrix)
    mat = property(get_matrix,doc='Return a matrix representing the rotation')

    def get_coordin(self):
        if not self.consistent: return None
        c,i = zip(self._coords,self._invs)[-1]
        return c[i]
    coordin = property(get_coordin, doc="the input coordinate system")

    def get_coordout(self):
        if not self.consistent: return None
        c,i = zip(self._coords,self._invs)[0]
        return c[not i]
    coordout = property(get_coordout, doc="the output coordinate system")

    def get_coordin_str(self):
        return coordname.get(self.coordin,'')
    coordinstr = property(get_coordin_str, doc="the input coordinate system in str")

    def get_coordout_str(self):
        return coordname.get(self.coordout,'')
    coordoutstr = property(get_coordout_str, doc="the output coordinate system in str")

    def get_rots(self):
        return self._rots
    rots = property(get_rots, doc="the sequence of rots defining")
    
    def get_coords(self):
        return self._coords
    coords = property(get_coords, doc="the sequence of coords")
    
    def do_rot(self,i):
        return not npy.allclose(self.rots[i],npy.zeros(3),rtol=0.,atol=1.e-15)

    def angle_ref(self,*args,**kwds):
        """Compute the angle between transverse reference direction of initial and final frames
        For example, if angle of polarisation is psi in initial frame, it will be psi+angle_ref in final
        frame.
        Input:
          - direction or vector (see Rotator.__call__)
        Keywords:
          - lonlat: if True, assume input is longitude,latitude in degrees. Otherwise,
                    theta,phi in radian. Default: False
          - inv: if True, use the inverse transforms. Default: False
        Return:
          - angle in radian (a scalar or an array if input is a sequence of direction/vector)
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
                v = ang2vec(arg[0],arg[1],lonlat=lonlat)
            else:
                v = arg
        elif len(args) == 2:
            v = ang2vec(args[0],args[1],lonlat=lonlat)
        elif len(args) == 3:
            v = args
        else:
            raise TypeError('Either 1, 2 or 3 arguments accepted')
        vp = R(v,inv=inv)
        north_pole = R([0.,0.,1.],inv=inv)
        sinalpha = north_pole[0]*vp[1]-north_pole[1]*vp[0]
        cosalpha = north_pole[2] - vp[2]*npy.dot(north_pole,vp)
        return npy.arctan2(sinalpha,cosalpha)

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
    """Rotate a vector (or a list of vectors) using the euler matrix
    given as first argument.

    Usage: vec_rot = rotateVector(rotmat,vec,vy=None,vz=None,do_rot=True)

    - rotmat : the 3x3 rotation matrix
    - vec[0], vec[1], vec[2] : vx,vy,vz (can be vectors)
      or: vec, vy, vz : vx,vy,vz if vy and vz are given.
    - do_rot: really perform the rotation if True, do nothing if false

    Return: vx,vy,vz
    """
    if vy is None and vz is None:
       if do_rot: return npy.tensordot(rotmat,vec,axes=(1,0))
       else: return vec
    elif vy is not None and vz is not None:
       if do_rot: return npy.tensordot(rotmat,npy.array([vec,vy,vz]),axes=(1,0))
       else: return vec,vy,vz
    else:
       raise TypeError("You must give either vec only or vec, vy "
                       "and vz parameters")

def rotateDirection(rotmat,theta,phi=None,do_rot=True,lonlat=False):
   """Rotate the direction pointed by theta,phi using the rotation matrix
   given as first argument.
   
   Usage: dir_rot = rotateDirection(rotmat,theta,phi=None,do_rot=True)
   
   - rotmat : the 3x3 rotation matrix
   - theta[0],theta[1] : theta, phi (can be vectors)
     or: theta, phi : theta, phi if phi is given.
   - do_rot: really perform the rotation if True, do nothing if false
   
   Return: theta_rot,phi_rot
   """
   vx,vy,vz=rotateVector(rotmat,ang2vec(theta,phi,lonlat=lonlat),do_rot=do_rot)
   return vec2ang(vx,vy,vz,lonlat=lonlat)

def vec2ang(vec,vy=None,vz=None,lonlat=False):
   """Transform a vector to a direction given by theta,phi.
   """
   if vy is None and vz is None:
      vx,vy,vz = vec
   elif vy is not None and vz is not None:
      vx=vec
   else:
      raise TypeError("You must either give both vy and vz or none of them")
   r = npy.sqrt(vx**2+vy**2+vz**2)
   theta = npy.arccos(vz/r)
   phi = npy.arctan2(vy,vx)
   if lonlat:
       return npy.asarray([npy.degrees(phi),90-npy.degrees(theta)])
   else:
       return npy.asarray([theta,phi])

def ang2vec(theta,phi=None,lonlat=False):
   """Transform a direction theta,phi to a unit vector.
   """
   if phi is None:
      theta,phi=theta
   if lonlat:
       lon,lat=theta,phi
       theta,phi = npy.pi/2.-npy.radians(lat), npy.radians(lon)
   ct,st,cp,sp = npy.cos(theta),npy.sin(theta),npy.cos(phi),npy.sin(phi)
   return npy.asarray([st*cp,st*sp,ct])

def angdist(dir1,dir2,lonlat=False):
    """Return the angular distance between dir1 and dir2.
    """
    if hasattr(lonlat,'__len__') and len(lonlat) == 2:
        lonlat1,lonlat2 = lonlat
    else:
        lonlat1=lonlat2=lonlat
    if len(dir1) == 2: # theta,phi or lonlat, convert to vec
        vec1 = npy.asarray(ang2vec(dir1,lonlat=lonlat1))
    else:
        vec1 = npy.asarray(dir1)
    if vec1.ndim == 1:
        vec1 = npy.expand_dims(vec1,-1)
    if len(dir2) == 2:
        vec2 = npy.asarray(ang2vec(dir2,lonlat=lonlat1)).T
    else:
        vec2 = npy.asarray(dir2)
    if vec2.ndim == 1:
        vec2 = npy.expand_dims(vec2,-1)
    # compute scalar product
    pscal = (vec1*vec2).sum(axis=0)
    return npy.arccos(pscal)


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
    if type(c) is not str:
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
   if deg: convert=npy.pi/180.
   else: convert=1.
   if rot is None:
      rot=npy.zeros(3)
   else:
      rot=npy.array(rot,npy.float64).flatten()*convert
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
   if not npy.allclose(rot,npy.zeros(3),rtol=0.,atol=1.e-15):
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
      matconv = npy.identity(3)
      do_conv = False        
   else:
      eps = 23.452294 - 0.0130125 - 1.63889E-6 + 5.02778E-7
      eps = eps * npy.pi / 180.
      
      # ecliptic to galactic
      e2g = npy.array([[-0.054882486, -0.993821033, -0.096476249],
                   [ 0.494116468, -0.110993846,  0.862281440],
                   [-0.867661702, -0.000346354,  0.497154957]])
      
      # ecliptic to equatorial
      e2q = npy.array([[1.,     0.    ,      0.         ],
                   [0., npy.cos( eps ), -1. * npy.sin( eps )],
                   [0., npy.sin( eps ),    npy.cos( eps )   ]])
      
      # galactic to ecliptic
      g2e = npy.linalg.inv(e2g)
      
      # galactic to equatorial                   
      g2q = npy.dot(e2q , g2e)
      
      # equatorial to ecliptic
      q2e = npy.linalg.inv(e2q)
      
      # equatorial to galactic
      q2g = npy.dot(e2g , q2e)
   
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

   PI = npy.pi
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
   sb = npy.sin(b)
   cb = npy.cos(b)
   cbsa = cb * npy.sin(a)
   b  = -stheta[i] * cbsa + ctheta[i] * sb
   #bo    = math.asin(where(b<1.0, b, 1.0)*deg_to_rad)
   bo    = npy.arcsin(b)*deg_to_rad
   #
   a = npy.arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * npy.cos(a) )
   ao = npy.fmod( (a+psi[i]+fourpi), twopi) * deg_to_rad
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
      convert = npy.pi/180.
      
   c1 = npy.cos(a1*convert)
   s1 = npy.sin(a1*convert)
   c2 = npy.cos(a2*convert)
   s2 = npy.sin(a2*convert)
   c3 = npy.cos(a3*convert)
   s3 = npy.sin(a3*convert)
        
   if ZYX:
      m1 = npy.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = npy.array([[ c2,  0, s2],
                  [  0,  1,  0],
                  [-s2,  0, c2]]) # around   y

      m3 = npy.array([[  1,  0,  0],
                  [  0, c3,-s3],
                  [  0, s3, c3]]) # around   x

   elif Y:
      m1 = npy.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = npy.array([[ c2,  0, s2],
                  [  0,  1,  0],
                  [-s2,  0, c2]]) # around   y

      m3 = npy.array([[ c3,-s3,  0],
                  [ s3, c3,  0],
                  [  0,  0,  1]]) # around   z

   else:
      m1 = npy.array([[ c1,-s1,  0],
                  [ s1, c1,  0],
                  [  0,  0,  1]]) # around   z

      m2 = npy.array([[  1,  0,  0],
                  [  0, c2,-s2],
                  [  0, s2, c2]]) # around   x

      m3 = npy.array([[ c3,-s3,  0],
                  [ s3, c3,  0],
                  [  0,  0,  1]]) # around   z

   M = npy.dot(m3.T,npy.dot(m2.T,m1.T)) 

   return M


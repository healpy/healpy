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
"""A wrapper for the PSHT library.
Beware, this one only works with healpix maps, and double precision maps
and alms.
"""
cimport numpy as c_numpy
import numpy as np
nm = np
np = np

# Numpy must be initialized
c_numpy.import_array()

cdef extern from "stdlib.h":
 void* malloc(long)
 void* calloc(long,long)
 void free(void*)
 ctypedef long size_t

cdef extern from "string.h":
 void *memcpy(void *dst,  void *src, size_t l)

cimport python
cdef extern from "psht.h":
 ctypedef struct pshtd_joblist
 ctypedef struct psht_geom_info
 ctypedef python.Py_complex pshtd_cmplx
 ctypedef struct psht_alm_info

 void 	pshtd_make_joblist (pshtd_joblist **joblist)
 void 	pshtd_clear_joblist (pshtd_joblist *joblist)
 void 	pshtd_destroy_joblist (pshtd_joblist *joblist)
 void 	pshtd_add_job_alm2map (pshtd_joblist *joblist,  pshtd_cmplx *alm, double *map, int add_output)
 void 	pshtd_add_job_map2alm (pshtd_joblist *joblist,  double *map, pshtd_cmplx *alm, int add_output)
 void 	pshtd_add_job_alm2map_pol (pshtd_joblist *joblist,  pshtd_cmplx *almT,  pshtd_cmplx *almG,  pshtd_cmplx *almC, double *mapT, double *mapQ, double *mapU, int add_output)
 void 	pshtd_add_job_map2alm_pol (pshtd_joblist *joblist,  double *mapT,  double *mapQ,  double *mapU, pshtd_cmplx *almT, pshtd_cmplx *almG, pshtd_cmplx *almC, int add_output)
 void 	pshtd_add_job_alm2map_spin (pshtd_joblist *joblist,  pshtd_cmplx *alm1,  pshtd_cmplx *alm2, double *map1, double *map2, int spin, int add_output)
 void 	pshtd_add_job_map2alm_spin (pshtd_joblist *joblist,  double *map1,  double *map2, pshtd_cmplx *alm1, pshtd_cmplx *alm2, int spin, int add_output)
 void 	pshtd_add_job_alm2map_deriv1 (pshtd_joblist *joblist,  pshtd_cmplx *alm, double *mapdtheta, double *mapdphi, int add_output)
 void 	pshtd_execute_jobs (pshtd_joblist *joblist,  psht_geom_info *geom_info,  psht_alm_info *alm_info)
 void 	psht_make_geom_info (int nrings,  int *nph,  int *ofs,  int *stride,  double *phi0,  double *theta,  double *weight, psht_geom_info **geom_info)
 void 	psht_destroy_geom_info (psht_geom_info *info)
 void 	psht_make_healpix_geom_info (int nside, int stride, psht_geom_info **geom_info)
 void 	psht_make_weighted_healpix_geom_info (int nside, int stride,  double *weight, psht_geom_info **geom_info)
 void 	psht_make_gauss_geom_info (int nrings, int nphi, int stride, psht_geom_info **geom_info)
 void 	psht_make_ecp_geom_info (int nrings, int nphi, double phi0, int stride, psht_geom_info **geom_info)
 void  psht_make_triangular_alm_info (int lmax, int mmax, int stride, psht_alm_info **alm_info)
 void 	psht_destroy_alm_info (psht_alm_info *info)

class pshtError(Exception):
 pass

def _howManyMaps(rmap):
 if isinstance(rmap,np.ndarray):
   if len(rmap.shape)==1:
     return 1
   else:
     return rmap.shape[0]
 if isinstance(rmap,(list,tuple)):
   if isinstance(rmap[0],(list,tuple,np.ndarray)):
     return len(rmap)
   else:
     return 1

cdef class job:
 """A wrapper for the PSHT job list.
Beware, this one only works with healpix maps, and double precision maps and alms.
job(int nside, int lmax=-1, int mmax=-1, int stridemap = 1, int stridealm = 1)
Create a new job list. 
 nside : nside
 lmax : if ==-1, set to 3 * nside - 1
 mmax : if ==-1, set to lmax
"""
 cdef pshtd_joblist *jb
 cdef psht_geom_info *geom
 cdef psht_alm_info *almi
 cdef object storein, storeout
 cdef int lmax,nside,mmax

 def __init__(self,int nside, int lmax=-1, int mmax=-1, int stridemap = 1, int stridealm = 1):
   if lmax == -1:
     lmax = 3*nside-1
   if mmax == -1:
     mmax = lmax 
   psht_make_healpix_geom_info (nside, stridemap, &self.geom)
   psht_make_triangular_alm_info (lmax, mmax, stridealm, &self.almi)
   pshtd_make_joblist (&self.jb)
   self.storein = []
   self.storeout = []
   self.lmax = lmax
   self.mmax = mmax
   self.nside = nside

 def _mapsize(self):
   return 12*self.nside**2
 def _almsize(self):
   return ((self.lmax+1)*(self.lmax+2))/2

 def _newmap(self):
   return np.zeros(self._mapsize())
 def _newalm(self):
   return np.zeros(self._almsize(),dtype=np.complex)

 def _testmapalm(self,rmap,alm):
   if len(rmap)!=self._mapsize():
     raise pshtError("Wrong size for map (expected %d, got %d)"%(self._mapsize(),len(rmap)))
   if len(alm)!=self._almsize():
     raise pshtError("Wrong size for alm (expected %d, got %d)"%(self._almsize(),len(alm)))

 def _testdo(self):
   if len(self.storein)!=0:
     raise pshtError("Can't do !")

 def add_alm2map(self, alm,rmap = None):
   """add_alm2map(alm,rmap = None)
Add a alm2map transform to the job list.
 alm : the input alm. Can be either a single alm, or a table of three T,E,B alms.
 rmap : if != None, inplace add the result of the alm2map transform to this map (or table of maps).
"""
   if _howManyMaps(alm)==3:
     self._add_alm2map_pol(alm,rmap)
   else:
     self._add_alm2map(alm,rmap)

 def do_alm2map(self,alm,rmap=None):
   """do_alm2map(alm,rmap = None)
Add a alm2map transform to the job list and execute the joblist. 
Fail if the joblist is not empty. 
 alm : the input alm. Can be either a single alm, or a table of three T,E,B alms.
 rmap : if != None, inplace add the result of the alm2map transform to this map (or table of maps).
"""
   self._testdo()
   self.add_alm2map(alm,rmap)
   return self.execute()[0]

 def _add_alm2map(self, alm,rmap = None):
   cdef double *_map
   cdef pshtd_cmplx *_alm

   alm_proxy=c_numpy.PyArray_ContiguousFromAny(alm,c_numpy.NPY_COMPLEX128,1,1)
   if rmap !=None:
     rmap_proxy=c_numpy.PyArray_ContiguousFromAny(rmap,c_numpy.NPY_DOUBLE,1,1)
     addo = 1
   else:
     rmap_proxy = self._newmap()
     addo = 0

   self._testmapalm(rmap_proxy,alm_proxy)

   _alm = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm_proxy)
   _map = <double*> c_numpy.PyArray_DATA(rmap_proxy)
   pshtd_add_job_alm2map (self.jb, _alm, _map, addo)
   self.storein += [alm_proxy]
   self.storeout += [rmap_proxy]

 def _add_alm2map_pol(self, alm,rmap = None):
   cdef double *_mapT,*_mapQ,*_mapU
   cdef pshtd_cmplx *_almT,*_almG,*_almC

   almT_proxy=c_numpy.PyArray_ContiguousFromAny(alm[0],c_numpy.NPY_COMPLEX128,1,1)
   almG_proxy=c_numpy.PyArray_ContiguousFromAny(alm[1],c_numpy.NPY_COMPLEX128,1,1)
   almC_proxy=c_numpy.PyArray_ContiguousFromAny(alm[2],c_numpy.NPY_COMPLEX128,1,1)
   if rmap !=None:
     if _howManyMaps(rmap)!=3:
       raise pshtError("not enough maps")
     rmapT_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[0],c_numpy.NPY_DOUBLE,1,1)
     rmapQ_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[1],c_numpy.NPY_DOUBLE,1,1)
     rmapU_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[2],c_numpy.NPY_DOUBLE,1,1)
     addo = 1
   else:
     rmapT_proxy = self._newmap()
     rmapQ_proxy = self._newmap()
     rmapU_proxy = self._newmap()
     addo = 0

   self._testmapalm(rmapT_proxy,almT_proxy)
   self._testmapalm(rmapQ_proxy,almG_proxy)
   self._testmapalm(rmapU_proxy,almC_proxy)

   _almT = <pshtd_cmplx*> c_numpy.PyArray_DATA(almT_proxy)
   _almG = <pshtd_cmplx*> c_numpy.PyArray_DATA(almG_proxy)
   _almC = <pshtd_cmplx*> c_numpy.PyArray_DATA(almC_proxy)
   _mapT = <double*> c_numpy.PyArray_DATA(rmapT_proxy)
   _mapQ = <double*> c_numpy.PyArray_DATA(rmapQ_proxy)
   _mapU = <double*> c_numpy.PyArray_DATA(rmapU_proxy)

   pshtd_add_job_alm2map_pol (self.jb, _almT,_almG,_almC, _mapT,_mapQ,_mapU, addo)
   self.storein += [(almT_proxy,almG_proxy,almC_proxy)]
   self.storeout += [(rmapT_proxy,rmapQ_proxy,rmapU_proxy)]

 def add_map2alm(self, rmap ,alm = None):
   """add_map2alm(rmap ,alm = None)
Add a alm2map transform to the job list.
 rmap : the input map. Can be either a single map, or a table of three T,Q,U maps.
 alm : if != None, inplace add the result of the map2alm transform to this alm (or table of alms).
"""
   if _howManyMaps(alm)==3:
     self._add_map2alm_pol(rmap,alm)
   else:
     self._add_map2alm(rmap,alm)

 def do_map2alm(self,rmap,alm=None):
   """do_map2alm(rmap,alm=None):
Add a alm2map transform to the job list and execute the joblist. 
Fail if the joblist is not empty.
 rmap : the input map. Can be either a single map, or a table of three T,Q,U maps.
 alm : if != None, inplace add the result of the map2alm transform to this alm (or table of alms).
"""
   self._testdo()
   self.add_map2alm(rmap,alm)
   return self.execute()[0]

 def _add_map2alm(self, rmap,alm = None):
   cdef double *_map
   cdef pshtd_cmplx *_alm

   rmap_proxy=c_numpy.PyArray_ContiguousFromAny(rmap,c_numpy.NPY_DOUBLE,1,1)

   if alm !=None:
     alm_proxy=c_numpy.PyArray_ContiguousFromAny(alm,c_numpy.NPY_COMPLEX128,1,1)
     addo = 1
   else:
     alm_proxy = self._newalm()
     addo = 0

   self._testmapalm(rmap_proxy,alm_proxy)

   _alm = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm_proxy)
   _map = <double*> c_numpy.PyArray_DATA(rmap_proxy)
   pshtd_add_job_map2alm (self.jb, _map, _alm, addo)
   self.storeout += [alm_proxy]
   self.storein += [rmap_proxy]

 def _add_map2alm_pol(self, rmap,alm = None):
   cdef double *_mapT,*_mapQ,*_mapU
   cdef pshtd_cmplx *_almT,*_almG,*_almC

   rmapT_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[0],c_numpy.NPY_DOUBLE,1,1)
   rmapQ_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[1],c_numpy.NPY_DOUBLE,1,1)
   rmapU_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[2],c_numpy.NPY_DOUBLE,1,1)

   if alm !=None:
     if _howManyMaps(alm)!=3:
       raise pshtError("not enough alm")
     almT_proxy=c_numpy.PyArray_ContiguousFromAny(alm[0],c_numpy.NPY_COMPLEX128,1,1)
     almG_proxy=c_numpy.PyArray_ContiguousFromAny(alm[1],c_numpy.NPY_COMPLEX128,1,1)
     almC_proxy=c_numpy.PyArray_ContiguousFromAny(alm[2],c_numpy.NPY_COMPLEX128,1,1)
     addo = 1
   else:
     almT_proxy = self._newalm()
     almG_proxy = self._newalm()
     almC_proxy = self._newalm()
     addo = 0

   self._testmapalm(rmapT_proxy,almT_proxy)
   self._testmapalm(rmapQ_proxy,almG_proxy)
   self._testmapalm(rmapU_proxy,almC_proxy)

   _almT = <pshtd_cmplx*> c_numpy.PyArray_DATA(almT_proxy)
   _almG = <pshtd_cmplx*> c_numpy.PyArray_DATA(almG_proxy)
   _almC = <pshtd_cmplx*> c_numpy.PyArray_DATA(almC_proxy)
   _mapT = <double*> c_numpy.PyArray_DATA(rmapT_proxy)
   _mapQ = <double*> c_numpy.PyArray_DATA(rmapQ_proxy)
   _mapU = <double*> c_numpy.PyArray_DATA(rmapU_proxy)

   pshtd_add_job_map2alm_pol (self.jb,  _mapT,_mapQ,_mapU,_almT,_almG,_almC, addo)
   self.storeout += [(almT_proxy,almG_proxy,almC_proxy)]
   self.storein += [(rmapT_proxy,rmapQ_proxy,rmapU_proxy)]

 def add_alm2map_spin(self, alm,int spin,rmap = None):
   """add_alm2map_spin(alm,int spin,rmap = None)
Add a alm2map_spin transform to the job list.
 alm : the input alms as a table of two alms
 spin : the spin of the alms (can be 1 2 or 3)
 ramp : if != None, inplace add the result of the alm2map_spin transform to this table of maps.
"""
   cdef double *_map1,*_map2
   cdef pshtd_cmplx *_alm1,*_alm2

   if _howManyMaps(alm)!=2:
     raise pshtError("not enough alm")

   alm1_proxy=c_numpy.PyArray_ContiguousFromAny(alm[0],c_numpy.NPY_COMPLEX128,1,1)
   alm2_proxy=c_numpy.PyArray_ContiguousFromAny(alm[1],c_numpy.NPY_COMPLEX128,1,1)

   if rmap !=None:
     if _howManyMaps(rmap)!=2:
       raise pshtError("not enough map")
     rmap1_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[0],c_numpy.NPY_DOUBLE,1,1)
     rmap2_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[1],c_numpy.NPY_DOUBLE,1,1)
     addo = 1
   else:
     rmap1_proxy = self._newmap()
     rmap2_proxy = self._newmap()
     addo = 0

   self._testmapalm(rmap1_proxy,alm1_proxy)
   self._testmapalm(rmap2_proxy,alm2_proxy)

   _alm1 = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm1_proxy)
   _alm2 = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm2_proxy)
   _map1 = <double*> c_numpy.PyArray_DATA(rmap1_proxy)
   _map2 = <double*> c_numpy.PyArray_DATA(rmap2_proxy)

   pshtd_add_job_alm2map_spin(self.jb, _alm1,_alm2, _map1,_map2,spin, addo)
   self.storein += [(alm1_proxy,alm2_proxy)]
   self.storeout += [(rmap1_proxy,rmap2_proxy)]

 def add_map2alm_spin(self, rmap,int spin,alm = None):
   """add_map2alm_spin(rmap,int spin,alm = None)
Add a map2alm_spin transform to the job list.
 rmap : the input maps as a table of two maps
 spin : the spin of the alms (can be 1 2 or 3)
 ramp : if != None, inplace add the result of the map2alm_spin transform to this table of alms.
"""
   cdef double *_map1,*_map2
   cdef pshtd_cmplx *_alm1,*_alm2

   if _howManyMaps(rmap)!=2:
     raise pshtError("not enough map")

   rmap1_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[0],c_numpy.NPY_DOUBLE,1,1)
   rmap2_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[1],c_numpy.NPY_DOUBLE,1,1)

   if alm !=None:
     if _howManyMaps(alm)!=2:
       raise pshtError("not enough alm")
     alm1_proxy=c_numpy.PyArray_ContiguousFromAny(alm[0],c_numpy.NPY_COMPLEX128,1,1)
     alm2_proxy=c_numpy.PyArray_ContiguousFromAny(alm[1],c_numpy.NPY_COMPLEX128,1,1)
     addo = 1
   else:
     alm1_proxy = self._newalm()
     alm2_proxy = self._newalm()
     addo = 0

   self._testmapalm(rmap1_proxy,alm1_proxy)
   self._testmapalm(rmap2_proxy,alm2_proxy)

   _alm1 = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm1_proxy)
   _alm2 = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm2_proxy)
   _map1 = <double*> c_numpy.PyArray_DATA(rmap1_proxy)
   _map2 = <double*> c_numpy.PyArray_DATA(rmap2_proxy)

   pshtd_add_job_map2alm_spin(self.jb, _map1,_map2,_alm1,_alm2, spin, addo)
   self.storeout += [(alm1_proxy,alm2_proxy)]
   self.storein += [(rmap1_proxy,rmap2_proxy)]

 def add_alm2map_der1(self, alm,rmap = None):
   """add_alm2map_der1(alm,rmap = None)
Add a alm2map_der1 transform to the job list.
 alm : the input alm
 ramp : if != None, inplace add the result of the alm2map_der transform to this table of maps.
"""    
   cdef double *_map1,*_map2
   cdef pshtd_cmplx *_alm

   alm_proxy=c_numpy.PyArray_ContiguousFromAny(alm,c_numpy.NPY_COMPLEX128,1,1)

   if rmap !=None:
     if _howManyMaps(rmap)!=2:
       raise pshtError("not enough map")
     rmap1_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[0],c_numpy.NPY_DOUBLE,1,1)
     rmap2_proxy=c_numpy.PyArray_ContiguousFromAny(rmap[1],c_numpy.NPY_DOUBLE,1,1)
     addo = 1
   else:
     rmap1_proxy = self._newmap()
     rmap2_proxy = self._newmap()
     addo = 0

   self._testmapalm(rmap1_proxy,alm_proxy)
   self._testmapalm(rmap2_proxy,alm_proxy)

   _alm = <pshtd_cmplx*> c_numpy.PyArray_DATA(alm_proxy)
   _map1 = <double*> c_numpy.PyArray_DATA(rmap1_proxy)
   _map2 = <double*> c_numpy.PyArray_DATA(rmap2_proxy)

   pshtd_add_job_alm2map_deriv1(self.jb, _alm, _map1,_map2,addo)
   self.storein += [alm_proxy]
   self.storeout += [(rmap1_proxy,rmap2_proxy)]

 def execute(self):
   """execute()
Execute all the jobs in the joblist. Empty the joblist. 
Returns a list of the results"""
   pshtd_execute_jobs (self.jb, self.geom, self.almi)
   pshtd_clear_joblist(self.jb)
   self.storein = []
   res = self.storeout
   self.storeout=[]
   return res

 def __dealloc__(self):
   self.storein = []
   self.storeout = []
   pshtd_destroy_joblist(self.jb)
   psht_destroy_geom_info(self.geom)
   psht_destroy_alm_info (self.almi)

def __tlm(lmax,mmax):
 if lmax == None:
   lmax = -1
 if mmax == None:
   mmax = -1
 return lmax,mmax

# def alm2map(alm,nside,lmax=None,mmax=None):
#   """alm2map(alm,nside,lmax=None,mmax=None)
# Computes an Healpix map and its first derivatives given the alm, using PSHT.

#      The alm are given as a complex array. You can specify lmax
#      and mmax, or they will be computed from array size (assuming
#      lmax==mmax).

#      Parameters:
#      - alm: a complex array of alm. Size must be of the form
#             size=mmax(lmax-mmax+1)/2+lmax
#      - nside: the nside of the output map.
#      - lmax: explicitly define lmax (needed if mmax!=lmax)
#      - mmax: explicitly define mmax (needed if mmax!=lmax)

#      Return: an Healpix map in RING scheme at nside.
#      """
#   lmax,mmax = __tlm(lmax,mmax)
#   if lmax==-1:
#     if _howManyMaps(alm)==3:
#       lmax = __getlmax(len(alm[0]),mmax)
#     else:
#       lmax = __getlmax(len(alm),mmax)

#   jb = job(nside,lmax,mmax)
#   rmap = jb.do_alm2map(alm)
#   return rmap

# def map2alm(rmap,lmax=None,mmax=None,iter=1, use_weights=False, regression=True):
#   """map2alm(rmap,lmax=None,mmax=None,iter=1, use_weights=False, regression=True)
# Computes the alm of an Healpix map using PSHT.

#   Input:
#     - m: a ndarray (not polarized) or a list of 3 ndarray (polarized)
#   Parameters:
#     - lmax : maximum l of the power spectrum. Default: 3*nside-1
#     - mmax : maximum m of the alm. Default: lmax
#     - iter : number of iteration (default: 1)
#     - use_weights: UNUSED
#     - regression: if True, subtract map average before computing alm. Default: True.
#   Return:
#     - alm as one ndarray or a tuple of 3 ndarrays
#   """
#   lmax,mmax = __tlm(lmax,mmax)
#   pol = False
#   if _howManyMaps(rmap)==3:
#     nside=int(np.sqrt(len(rmap[0])/12.))
#     avg = np.mean(rmap[0])
#     rmap[0] = rmap[0] - avg
#     pol = True
#   else:
#     if regression:
#       avg = np.mean(rmap)
#       rmap = rmap - avg
#     nside=int(np.sqrt(len(rmap)/12.))

#   jb = job(nside,lmax,mmax)
#   alm = jb.do_map2alm(rmap)

#   for i  in range(iter):
#     back = jb.do_alm2map(alm)
#     if pol:
#       delta = [mm - bb for mm,bb in zip(rmap,back)]
#     else:
#       delta = rmap - back
#     alm = jb.do_map2alm(delta,alm)

#   if regression:
#     if pol:
#       alm[0][0] += avg*np.sqrt(4*np.pi)
#     else:
#       alm[0] += avg*np.sqrt(4*np.pi)
#   return alm

def alm2map_spin(alm,nside,spin,lmax=None,mmax=None):
 """alm2map_spin(alm,nside,spin,lmax=None,mmax=None)
Computes two healpix maps from two spinned alms using PSHT.

 Input:
   - rmap: a list of two alms
 Parameters:
   - nside : nside of the result maps
   - spin : spin of the alms (one of 1 2 or 3)
   - lmax : maximum l of the power spectrum. Default: 3*nside-1
   - mmax : maximum m of the alm. Default: lmax
 Return:
   - a list of two maps
 """
 lmax,mmax = __tlm(lmax,mmax)
 if lmax==-1:
   lmax = __getlmax(len(alm[0]),mmax)

 jb = job(nside,lmax,mmax)
 jb.add_alm2map_spin(alm, spin)
 res = jb.execute()
 return res[0]

def map2alm_spin(rmap,spin,lmax=None,mmax=None):
 """map2alm_spin(rmap,spin,lmax=None,mmax=None)
Computes the spinned alm of two Healpix maps using PSHT.

 Input:
   - rmap: a list of two maps
 Parameters:
   - spin : spin of the alms (one of 1 2 or 3)
   - lmax : maximum l of the power spectrum. Default: 3*nside-1
   - mmax : maximum m of the alm. Default: lmax
 Return:
   - a list of two alms
 """
 lmax,mmax = __tlm(lmax,mmax)
 nside=int(np.sqrt(len(rmap[0])/12.))
 jb = job(nside,lmax,mmax)
 jb.add_map2alm_spin(rmap, spin)
 res = jb.execute()
 return res[0]

def alm2map_der1(alm,nside,lmax=None,mmax=None):
 """alm2map_der1(alm,nside,lmax=None,mmax=None)
Computes two healpix maps corresponding to the derivatives of an alm using PSHT.

 Input:
   - alm : the alm
 Parameters:
   - nside : nside of the result maps
   - lmax : maximum l of the power spectrum. Default: 3*nside-1
   - mmax : maximum m of the alm. Default: lmax
 Return:
   - a list of two maps
 """
 lmax,mmax = __tlm(lmax,mmax)
 if lmax==-1:
   lmax = __getlmax(len(alm),mmax)
 jb = job(nside,lmax,mmax)
 jb.add_alm2map_der1(alm)
 res = jb.execute()
 return res[0]

def __getlmax_from_hp(s,mmax=-1):
   if mmax >= 0:
       x=(2*s+mmax**2-mmax-2)/(2*mmax+2)
   else:
       x=(-3+np.sqrt(1+8*s))/2
   if x != np.floor(x):
       return -1
   else:
       return int(x)

try:
 import healpy as _hp
 __getlmax = _hp.Alm.getlmax
except ImportError:  
 __getlmax = __getlmax_from_hp

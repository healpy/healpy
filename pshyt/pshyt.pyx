#cython: embedsignature=True
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
#cython: embedsignature=True

cimport numpy as c_numpy
import numpy as np

import weakref

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

    void   pshtd_make_joblist (pshtd_joblist **joblist)
    void   pshtd_clear_joblist (pshtd_joblist *joblist)
    void   pshtd_destroy_joblist (pshtd_joblist *joblist)
    void   pshtd_add_job_alm2map (pshtd_joblist *joblist,  pshtd_cmplx *alm, double *map, int add_output)
    void   pshtd_add_job_map2alm (pshtd_joblist *joblist,  double *map, pshtd_cmplx *alm, int add_output)
    void   pshtd_add_job_alm2map_pol (pshtd_joblist *joblist,  pshtd_cmplx *almT,  pshtd_cmplx *almG,  pshtd_cmplx *almC, double *mapT, double *mapQ, double *mapU, int add_output)
    void   pshtd_add_job_map2alm_pol (pshtd_joblist *joblist,  double *mapT,  double *mapQ,  double *mapU, pshtd_cmplx *almT, pshtd_cmplx *almG, pshtd_cmplx *almC, int add_output)
    void   pshtd_add_job_alm2map_spin (pshtd_joblist *joblist,  pshtd_cmplx *alm1,  pshtd_cmplx *alm2, double *map1, double *map2, int spin, int add_output)
    void   pshtd_add_job_map2alm_spin (pshtd_joblist *joblist,  double *map1,  double *map2, pshtd_cmplx *alm1, pshtd_cmplx *alm2, int spin, int add_output)
    void   pshtd_add_job_alm2map_deriv1 (pshtd_joblist *joblist,  pshtd_cmplx *alm, double *mapdtheta, double *mapdphi, int add_output)
    void   pshtd_execute_jobs (pshtd_joblist *joblist,  psht_geom_info *geom_info,  psht_alm_info *alm_info)
    void   psht_make_geom_info (int nrings,  int *nph,  int *ofs,  int *stride,  double *phi0,  double *theta,  double *weight, psht_geom_info **geom_info)
    void   psht_destroy_geom_info (psht_geom_info *info)
    void   psht_make_healpix_geom_info (int nside, int stride, psht_geom_info **geom_info)
    void   psht_make_weighted_healpix_geom_info (int nside, int stride,  double *weight, psht_geom_info **geom_info)
    void   psht_make_gauss_geom_info (int nrings, int nphi, int stride, psht_geom_info **geom_info)
    void   psht_make_ecp_geom_info (int nrings, int nphi, double phi0, int stride, psht_geom_info **geom_info)
    void  psht_make_triangular_alm_info (int lmax, int mmax, int stride, psht_alm_info **alm_info)
    void   psht_destroy_alm_info (psht_alm_info *info)

cdef extern from "ylmgen_c.h":
    ctypedef double ylmgen_dbl2[2]
    ctypedef struct Ylmgen_C:
        int lmax, mmax, smax, nth, ith
        double *ylm
        ylmgen_dbl2 **lambda_wx
 
    void Ylmgen_init (Ylmgen_C *gen, int l_max, int m_max, int s_max, int spinrec,
                      double epsilon)
    void Ylmgen_set_theta (Ylmgen_C *gen, double *theta, int nth)
    void Ylmgen_destroy (Ylmgen_C *gen)
    void Ylmgen_prepare (Ylmgen_C *gen, int ith, int m)
    void Ylmgen_recalc_Ylm (Ylmgen_C *gen)
    void Ylmgen_recalc_lambda_wx (Ylmgen_C *gen, int spin) 
 

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
        rmap : if != None, inplace add the result of the alm2map_spin transform to this table of maps.
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
        rmap : if != None, inplace add the result of the map2alm_spin transform to this table of alms.
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
        rmap : if != None, inplace add the result of the alm2map_der transform to this table of maps.
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

cdef class Ylmgen:
    """Class wrapping the Ylmgen_C structure and functions
    to compute scalar and tensor Spherical Harmonics Ylm(theta, phi=0)
    and Wlm/Xlm(spin, theta, phi=0).
    """
    cdef Ylmgen_C * _ylmgen
    cdef object _ylm
    cdef object _lambda_wx

    def __init__(self, lmax, mmax=None, smax=0, spinrec=False,
                 epsilon=1.e-15):
        """Create a Ylmgen object, which can generate Ylm, Wlm, Xlm values.

        Parameters
        ----------
        lmax : int, scalar
          The maximum l for which Ylm/Xlm/Xlm will be computed.
        mmax : None or int, scalar, optional
          The maximum m for which Ylm/Wlm/Xlm will be computable.
          If None, mmax will be set to lmax. Default: None.
        smax : int, scalar, optional
          The maximum spin for which Ylm/Wlm/Xlm will be computables.
          Default : 0. (No spin harmonics)
        spinrec : bool, scalar
          If True, use recursion from spin=0 to compute spin=1, 2.
          Otherwise, use Wigner d-symbol. Default: False
        epsilon : double, scalar, optional
          The values of Ylm below epsilon may be considered as null in
          recursion formulae. Default: 1e-15

        Returns
        -------
        A new Ylmgen object

        Example
        -------
        >>> ylmgen = Ylmgen(10, smax=2, spinrec=True)
        >>> ylmgen.ylm(0, np.pi / 2)
        array([  2.82094792e-01,   2.99182751e-17,  -3.15391565e-01,
                -6.85513802e-17,   3.17356641e-01,   1.07417129e-16,
                -3.17846011e-01,  -1.46342212e-16,   3.18036967e-01,
                 1.85290566e-16,  -3.18130494e-01])
        >>> ylmgen.lambda_wx(2, 0, np.pi / 2)
        array([[  0.00000000e+00,   0.00000000e+00],
               [  0.00000000e+00,   0.00000000e+00],
               [  1.89234939e+00,   0.00000000e+00],
               [  6.85513802e-16,  -0.00000000e+00],
               [ -6.34713281e+00,  -0.00000000e+00],
               [ -3.00767960e-15,   0.00000000e+00],
               [  1.33495325e+01,   0.00000000e+00],
               [  7.90247945e-15,  -0.00000000e+00],
               [ -2.28986616e+01,  -0.00000000e+00],
               [ -1.63055698e-14,   0.00000000e+00],
               [  3.49943543e+01,   0.00000000e+00]])
        """
        cdef c_numpy.intp_t dims[1]
        lmax = int(lmax)
        smax = int(smax)
        spinrec = bool(spinrec)
        epsilon = float(epsilon)
        if mmax is None:
            mmax = lmax
        else:
            mmax = int(mmax)
        if lmax < 0:
            raise ValueError("lmax must be >=0")
        if mmax < 0 or mmax > lmax:
            raise ValueError("mmax must be between 0 and lmax")
        if smax < 0 or smax > 2:
            raise ValueError('only spin 0, 1, 2 are supported')
        self._ylmgen = <Ylmgen_C*>malloc(sizeof(Ylmgen_C))
        Ylmgen_init(self._ylmgen, lmax, mmax, smax, spinrec, epsilon)
        dims[0] = lmax + 1
        self._ylm = c_numpy.PyArray_SimpleNewFromData(1, dims,
                                                      c_numpy.NPY_DOUBLE,
                                                      <void*>self._ylmgen.ylm)
        self._lambda_wx = [None] * (smax + 1)

    cpdef set_theta(self, theta):
        """Initialize a table of theta, for which Ylm may be generated.

        Parameters
        ----------
        theta : float, array-like
          A 1-dim array or list of theta values

        See Also
        --------
        prepare, recalc_Ylm, recalc_lambda_wx.
        The methods ylm and lambda_wx automatically call this method.
        """
        cdef c_numpy.ndarray[c_numpy.float64_t, ndim = 1] theta_
        theta_ = np.ascontiguousarray(theta, dtype = np.float64)
        nth = theta_.size
        Ylmgen_set_theta(self._ylmgen, <double*>theta_.data, nth)

    cpdef prepare(self, int ith, int m):
        """Prepare the computation for theta[ith] and m.

        Parameters
        ----------
        ith : int, scalar
          The index in the theta array of the theta value for which to compute
          Ylm and/or Wlm/Xlm. See set_theta.
        m : int, scalar
          The value of m.

        See Also
        --------
        set_theta, recalc_Ylm, recalc_lambda_wx.
        The methods ylm and lambda_wx automatically call this method.
        """
        Ylmgen_prepare(self._ylmgen, ith, m)

    cpdef recalc_Ylm(self):
        """Recompute Ylm if needed.

        See Also
        --------
        set_theta, prepare, recalc_lambda_wx.
        The methods ylm and lambda_wx automatically call this method.
        """
        Ylmgen_recalc_Ylm(self._ylmgen)

    cpdef recalc_lambda_wx(self, spin):
        """Recompute Wlm/Xlm if needed.

        See Also
        --------
        set_theta, prepare, recalc_Ylm.
        The methods ylm and lambda_wx automatically call this method.
        """
        if spin > self._ylmgen.smax or spin <= 0:
            raise ValueError('spin must be between 0 and %d' % (self._ylmgen.smax))
        Ylmgen_recalc_lambda_wx(self._ylmgen, spin)

    cpdef get_ylm(self):
        """Return a copy of the internal buffer containing computed Ylm values.
        """
        return self._ylm.copy()
         
    cpdef ylm(self, int m, double theta):
        """Computes the Ylm(theta, phi) for l=0..lmax, m, theta, phi=0.
        Beware that values for l < m have no meaning and should not be used. 

        Parameters
        ----------
        m : int, scalar
          The value of m.
        theta : double, scalar
          The value of theta within [0,pi]

        Returns
        -------
        A copy of the internal buffer containing computed Ylm values, as a 
        one-dimensional numpy array of length lmax + 1.
        The index of the array correspond to l.
        Beware that values for l < m have no meaning and should not be used. 
        """
        cdef int lmax
        cdef double theta_arr[1]
        lmax = self._ylmgen.lmax
        theta_arr[0] = theta
        Ylmgen_set_theta(self._ylmgen, theta_arr, 1)
        Ylmgen_prepare(self._ylmgen, 0, m)
        Ylmgen_recalc_Ylm(self._ylmgen)
        return self._ylm.copy()

    cpdef get_lambda_wx(self, int spin):
        """Return a copy of the internal buffer containing computed Wlm
        and Xlm values. If they were never computed, returns None.
        """
        if spin > self._ylmgen.smax or spin <= 0:
            raise ValueError('spin must be between 1 and %d' % (self._ylmgen.smax))
        if self._ylmgen.lambda_wx[spin] == NULL:
            return None
        cdef c_numpy.ndarray[c_numpy.float64_t, ndim=2, mode='c'] arr
        cdef c_numpy.intp_t dims[2]
        if self._lambda_wx[spin] is None:
            dims[0] = self._ylmgen.lmax + 1
            dims[1] = 2
            arr = c_numpy.PyArray_SimpleNewFromData(
                            2, dims, c_numpy.NPY_FLOAT64,
                            <void*>self._ylmgen.lambda_wx[spin])
            self._lambda_wx[spin] = arr
        return self._lambda_wx[spin].copy()

    cpdef lambda_wx(self, int spin, int m, double theta):
        """Computes the Ylm(theta, phi) for l=0..lmax, m, theta, phi=0.
        Beware that values for l < m have no meaning and should not be used. 

        Parameters
        ----------
        spin : int, scalar
          The spin value (1 or 2).
        m : int, scalar
          The value of m.
        theta : double, scalar
          The value of theta within [0,pi]

        Returns
        -------
        A copy of the internal buffer containing computed Wlm/Xlm values, as a 
        two-dimensional numpy array of shape (lmax + 1, 2).
        The first index of the array correspond to l.
        The second index of the array correspond to W (index 0) or X (index 1).
        Beware that values for l < m have no meaning and should not be used.
        """
        cdef int lmax
        cdef double theta_arr[1]
        lmax = self._ylmgen.lmax
        theta_arr[0] = theta
        Ylmgen_set_theta(self._ylmgen, theta_arr, 1)
        Ylmgen_prepare(self._ylmgen, 0, m)
        Ylmgen_recalc_lambda_wx(self._ylmgen, spin)
        return self.get_lambda_wx(spin)

    def __dealloc__(self):
        for spin in range(self._ylmgen.smax):
            del self._lambda_wx[spin]
        Ylmgen_destroy(self._ylmgen)
        free(self._ylmgen)


try:
    import healpy as _hp
    __getlmax = _hp.Alm.getlmax
except ImportError:  
    __getlmax = __getlmax_from_hp

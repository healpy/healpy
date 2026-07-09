/*
 *  This file is part of Healpy.
 *
 *  Healpy is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpy is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpy; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about Healpy, see http://code.google.com/p/healpy
 */
/*

   This module provides Healpix functions to Python.
   It uses the healpix_cxx library.

*/

#include <Python.h>

#include "arr.h"
#include "healpix_base.h"
#include "healpix_map.h"
#include "_healpy_utils.h"

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "numpy/ndarrayobject.h"

/*
   ang2pix
*/
template<Healpix_Ordering_Scheme scheme>static void
  ufunc_ang2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp n=dimensions[0];

  npy_intp is1=steps[0],is2=steps[1],is3=steps[2], os=steps[3];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(npy_intp i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = pointing(*(double *)ip2,*(double *)ip3);
        ptg.normalize();
        *(int64 *)op = hb.ang2pix(ptg);
      } catch(PlanckError &e) {
        *(int64 *)op = -1;
      }
  }
}

/*
   pix2ang
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2ang(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = hb.pix2ang(*(int64 *)ip2);
        *(double *)op1 = ptg.theta;
        *(double *)op2 = ptg.phi;
      } catch (PlanckError & e) {
        *(double *)op1 = NAN;
        *(double *)op2 = NAN;
      }
    }
}

/*
   xyf2pix
*/
template<Healpix_Ordering_Scheme scheme>static void
  ufunc_xyf2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp n=dimensions[0];

  npy_intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3], os=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op=args[4];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(npy_intp i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op+=os)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        *(int64 *)op = hb.xyf2pix((int)*(int64 *)ip2,(int)*(int64 *)ip3,(int)*(int64 *)ip4);
      } catch(PlanckError &e) {
        *(int64 *)op = -1;
      }
  }
}

/*
   pix2xyf
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2xyf(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        int x, y, f;
        hb.pix2xyf(*(int64 *)ip2, x, y, f);
        *(int64 *)op1 = x;
        *(int64 *)op2 = y;
        *(int64 *)op3 = f;
      } catch (PlanckError & e) {
        *(int64 *)op1 = -1;
        *(int64 *)op2 = -1;
        *(int64 *)op3 = -1;
      }
    }
}

/*
   ring2nest
*/
static void
ufunc_ring2nest(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],os=steps[2];
  char *ip1=args[0], *ip2=args[1], *op=args[2];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, RING); }
      try {
        *(int64 *)op = hb.ring2nest(*(int64 *)ip2);
      } catch(PlanckError & e) {
        *(int64 *)op = -1;
      }
    }
}

/*
   nest2ring
*/
static void
 ufunc_nest2ring (char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],os=steps[2];
  char *ip1=args[0], *ip2=args[1], *op=args[2];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, NEST); }
      try {
        *(int64 *)op = hb.nest2ring(*(int64 *)ip2);
      } catch(PlanckError & e) {
        *(int64 *)op = -1;
      }
    }
}


/*
  pix2vec
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2vec(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        vec3 v = hb.pix2vec(*(int64 *)ip2);
        *(double *)op1 = v.x;
        *(double *)op2 = v.y;
        *(double *)op3 = v.z;
      } catch (PlanckError & e) {
        *(double *)op1 = NAN;
        *(double *)op2 = NAN;
        *(double *)op3 = NAN;
      }
    }
}

/*
  vec2pix
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_vec2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      vec3 v (*(double *)ip2,*(double *)ip3,*(double *)ip4);
      try {
        int64 ipix = hb.vec2pix(v);
        *(int64 *)op1 = ipix;
      } catch (PlanckError &e) {
        *(int64 *)op1 = -1;
      }
    }
}


/*
  get_interpol
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_get_interpol(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],is3=steps[2],
    os1=steps[3],os2=steps[4],os3=steps[5],os4=steps[6],
    os5=steps[7],os6=steps[8],os7=steps[9],os8=steps[10];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2],
    *op1=args[3],*op2=args[4],*op3=args[5],*op4=args[6],
    *op5=args[7],*op6=args[8],*op7=args[9],*op8=args[10];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3,
        op1+=os1,op2+=os2,op3+=os3,op4+=os4,
        op5+=os5,op6+=os6,op7+=os7,op8+=os8 )
    {
      fix_arr<int64,4> pix;
      fix_arr<double,4> wgt;
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = pointing(*(double*)ip2, *(double*)ip3);
        ptg.normalize();
        hb.get_interpol(ptg, pix, wgt);
        *(int64*)op1 = (int64)pix[0];
        *(int64*)op2 = (int64)pix[1];
        *(int64*)op3 = (int64)pix[2];
        *(int64*)op4 = (int64)pix[3];
        *(double*)op5 = wgt[0];
        *(double*)op6 = wgt[1];
        *(double*)op7 = wgt[2];
        *(double*)op8 = wgt[3];
      } catch (PlanckError &e) {
        *(int64*)op1 = -1;
        *(int64*)op2 = -1;
        *(int64*)op3 = -1;
        *(int64*)op4 = -1;
        *(double*)op5 = NAN;
        *(double*)op6 = NAN;
        *(double*)op7 = NAN;
        *(double*)op8 = NAN;
      }
    }
}


/*
  get_neighbors
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_get_neighbors(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],is2=steps[1],
    os1=steps[2],os2=steps[3],os3=steps[4],os4=steps[5],
    os5=steps[6],os6=steps[7],os7=steps[8],os8=steps[9];
  char *ip1=args[0], *ip2=args[1],
    *op1=args[2],*op2=args[3],*op3=args[4],*op4=args[5],
    *op5=args[6],*op6=args[7],*op7=args[8],*op8=args[9];

  Healpix_Base2 hb;
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2,
        op1+=os1,op2+=os2,op3+=os3,op4+=os4,
        op5+=os5,op6+=os6,op7+=os7,op8+=os8 )
    {
      fix_arr<int64,8> pix;
      hb.SetNside(*(int64*)ip1, scheme);
      try {
        hb.neighbors(*(int64*)ip2, pix);
        *(int64*)op1 = (int64)pix[0];
        *(int64*)op2 = (int64)pix[1];
        *(int64*)op3 = (int64)pix[2];
        *(int64*)op4 = (int64)pix[3];
        *(int64*)op5 = (int64)pix[4];
        *(int64*)op6 = (int64)pix[5];
        *(int64*)op7 = (int64)pix[6];
        *(int64*)op8 = (int64)pix[7];
      } catch (PlanckError & e) {
        *(int64*)op1 = -1;
        *(int64*)op2 = -1;
        *(int64*)op3 = -1;
        *(int64*)op4 = -1;
        *(int64*)op5 = -1;
        *(int64*)op6 = -1;
        *(int64*)op7 = -1;
        *(int64*)op8 = -1;
      }
    }
}

/*
  max_pixrad
*/
static void
  ufunc_max_pixrad(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  register npy_intp i, n=dimensions[0];
  register npy_intp is1=steps[0],os1=steps[1];
  char *ip1=args[0], *op1=args[1];

  Healpix_Base2 hb;
  int64 oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, op1+=os1)
    {
      int64 nside = *(int64*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, RING);
          /* RING and NEST should give the same result but use RING because NEST only allows power of 2 nside */
        }
      double max_pixrad = hb.max_pixrad();
      *(double *)op1 = max_pixrad;
    }
}


static char *docstring =
  "This module contains basic ufunc related to healpix pixelisation\n"
  "scheme, such as ang2pix, ring<->nest swapping, etc.\n"
  "\n"
  "Available ufunc: _ang2pix_ring, _ang2pix_nest, _pix2ang_ring,\n"
  "                 _pix2ang_nest, _ring2nest, _nest2ring,\n"
  "                 _get_interpol_ring, _get_interpol_nest.";

/* to define the ufunc */
static PyUFuncGenericFunction ang2pix_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_ang2pix<RING>
};
static PyUFuncGenericFunction ang2pix_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_ang2pix<NEST>
};
static PyUFuncGenericFunction pix2ang_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2ang<RING>
};
static PyUFuncGenericFunction pix2ang_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2ang<NEST>
};
static PyUFuncGenericFunction xyf2pix_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_xyf2pix<RING>
};
static PyUFuncGenericFunction xyf2pix_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_xyf2pix<NEST>
};
static PyUFuncGenericFunction pix2xyf_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2xyf<RING>
};
static PyUFuncGenericFunction pix2xyf_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2xyf<NEST>
};
static PyUFuncGenericFunction vec2pix_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_vec2pix<RING>
};
static PyUFuncGenericFunction vec2pix_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_vec2pix<NEST>
};
static PyUFuncGenericFunction pix2vec_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2vec<RING>
};
static PyUFuncGenericFunction pix2vec_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_pix2vec<NEST>
};
static PyUFuncGenericFunction ring2nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_ring2nest
};
static PyUFuncGenericFunction nest2ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_nest2ring
};
static PyUFuncGenericFunction get_interpol_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_get_interpol<RING>
};
static PyUFuncGenericFunction get_interpol_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_get_interpol<NEST>
};
static PyUFuncGenericFunction get_neighbors_ring_functions[] = {
  (PyUFuncGenericFunction) ufunc_get_neighbors<RING>
};
static PyUFuncGenericFunction get_neighbors_nest_functions[] = {
  (PyUFuncGenericFunction) ufunc_get_neighbors<NEST>
};
static PyUFuncGenericFunction max_pixrad_functions[] = {
  (PyUFuncGenericFunction) ufunc_max_pixrad
};


static void *const blank_data[] = { (void *)NULL };

static const char ang2pix_signatures[] = {
  NPY_INT64, NPY_DOUBLE, NPY_DOUBLE, NPY_INT64
};
static const char pix2ang_signatures[] = {
  NPY_INT64, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE
};
static const char xyf2pix_signatures[] = {
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64
};
static const char pix2xyf_signatures[] = {
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64
};
static const char pix2vec_signatures[] = {
  NPY_INT64, NPY_INT64, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE
};
static const char vec2pix_signatures[] = {
  NPY_INT64, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INT64
};
static const char ring2nest_signatures[] = {
  NPY_INT64, NPY_INT64, NPY_INT64
};
static const char get_interpol_signatures[] = {
  NPY_INT64, NPY_DOUBLE, NPY_DOUBLE,
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64,
  NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE
};
static const char get_neighbors_ring_signatures[] = {
  NPY_INT64, NPY_INT64, // input
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64, // output
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64 // output
};
static const char get_neighbors_nest_signatures[] = {
  NPY_INT64, NPY_INT64, // input
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64, // output
  NPY_INT64, NPY_INT64, NPY_INT64, NPY_INT64 // output
};
static const char max_pixrad_signatures[] = {
  NPY_INT64, NPY_DOUBLE
};

static int m_exec(PyObject *module) {
  if (PyModule_AddObjectRef(module, "_ang2pix_ring", PyUFunc_FromFuncAndData(
      ang2pix_ring_functions, blank_data,
      ang2pix_signatures, 1,
      3, 1, PyUFunc_None, "_ang2pix_ring",
      "nside,theta,phi [rad] -> ipix (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_ang2pix_nest", PyUFunc_FromFuncAndData(
      ang2pix_nest_functions, blank_data,
      ang2pix_signatures, 1,
      3, 1, PyUFunc_None, "_ang2pix_nest",
      "nside,theta,phi [rad] -> ipix (NEST)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2ang_ring", PyUFunc_FromFuncAndData(
      pix2ang_ring_functions, blank_data,
      pix2ang_signatures, 1,
      2, 2, PyUFunc_None, "_pix2ang_ring",
      "nside,ipix -> theta,phi [rad] (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2ang_nest", PyUFunc_FromFuncAndData(
      pix2ang_nest_functions, blank_data,
      pix2ang_signatures, 1,
      2, 2, PyUFunc_None, "_pix2ang_nest",
      "nside,ipix -> theta,phi [rad] (NEST)",0)) < 0)
    return -1;

  //=========

  if (PyModule_AddObjectRef(module, "_xyf2pix_ring", PyUFunc_FromFuncAndData(
      xyf2pix_ring_functions, blank_data,
      xyf2pix_signatures, 1,
      4, 1, PyUFunc_None, "_xyf2pix_ring",
      "nside,x,y,face -> ipix (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_xyf2pix_nest", PyUFunc_FromFuncAndData(
      xyf2pix_nest_functions, blank_data,
      xyf2pix_signatures, 1,
      4, 1, PyUFunc_None, "_xyf2pix_nest",
      "nside,x,y,face -> ipix (NEST)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2xyf_ring", PyUFunc_FromFuncAndData(
      pix2xyf_ring_functions, blank_data,
      pix2xyf_signatures, 1,
      2, 3, PyUFunc_None, "_pix2xyf_ring",
      "nside,ipix -> x,y,face (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2xyf_nest", PyUFunc_FromFuncAndData(
      pix2xyf_nest_functions, blank_data,
      pix2xyf_signatures, 1,
      2, 3, PyUFunc_None, "_pix2xyf_nest",
      "nside,ipix -> x,y,face (NEST)",0)) < 0)
    return -1;

  //=========

  if (PyModule_AddObjectRef(module, "_vec2pix_ring", PyUFunc_FromFuncAndData(
      vec2pix_ring_functions, blank_data,
      vec2pix_signatures, 1,
      4, 1, PyUFunc_None, "_vec2pix_ring",
      "nside,x,y,z -> ipix (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_vec2pix_nest", PyUFunc_FromFuncAndData(
      vec2pix_nest_functions, blank_data,
      vec2pix_signatures, 1,
      4, 1, PyUFunc_None, "_vec2pix_nest",
      "nside,x,y,z -> ipix (NEST)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2vec_ring", PyUFunc_FromFuncAndData(
      pix2vec_ring_functions, blank_data,
      pix2vec_signatures, 1,
      2, 3, PyUFunc_None, "_pix2vec_ring",
      "nside,ipix -> x,y,z (RING)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_pix2vec_nest", PyUFunc_FromFuncAndData(
      pix2vec_nest_functions, blank_data,
      pix2vec_signatures, 1,
      2, 3, PyUFunc_None, "_pix2vec_nest",
      "nside,ipix -> x,y,z (NEST)",0)) < 0)
    return -1;

  //=============

  if (PyModule_AddObjectRef(module, "_ring2nest", PyUFunc_FromFuncAndData(
      ring2nest_functions, blank_data,
      ring2nest_signatures, 1,
      2, 1, PyUFunc_None, "_ring2nest",
      "ipix(ring) -> ipix(nest)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_nest2ring", PyUFunc_FromFuncAndData(
      nest2ring_functions, blank_data,
      ring2nest_signatures, 1,
      2, 1, PyUFunc_None, "_nest2ring",
      "ipix(nest) -> ipix(ring)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_get_interpol_ring", PyUFunc_FromFuncAndData(
      get_interpol_ring_functions, blank_data,
      get_interpol_signatures, 1,
      3, 8, PyUFunc_None, "_get_interpol_ring",
      "nside,theta,phi->4 nearest pixels+4weights",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_get_interpol_nest", PyUFunc_FromFuncAndData(
      get_interpol_nest_functions, blank_data,
      get_interpol_signatures, 1,
      3, 8, PyUFunc_None, "_get_interpol_nest",
      "nside,theta,phi->4 nearest pixels+4weights",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_get_neighbors_ring", PyUFunc_FromFuncAndData(
      get_neighbors_ring_functions, blank_data,
      get_neighbors_ring_signatures, 1,
      2, 8, PyUFunc_None, "_get_neigbors_ring",
      "nside, ipix [rad] -> 8 neighbors",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "_get_neighbors_nest", PyUFunc_FromFuncAndData(
      get_neighbors_nest_functions, blank_data,
      get_neighbors_nest_signatures, 1,
      2, 8, PyUFunc_None, "_get_neigbors_nest",
      "nside, ipix [rad] -> 8 neighbors",0)) < 0)
    return -1;


  if (PyModule_AddObjectRef(module, "_max_pixrad", PyUFunc_FromFuncAndData(
      max_pixrad_functions, blank_data,
      max_pixrad_signatures, 1,
      1, 1, PyUFunc_None, "max_pixrad",
      "nside -> max_distance to pixel corners from center)",0)) < 0)
    return -1;

  if (PyModule_AddObjectRef(module, "UNSEEN", PyFloat_FromDouble(Healpix_undef)) < 0)
    return -1;

  return 0;
}

static PyModuleDef_Slot moduledef_slots[] = {
  {Py_mod_exec, (void *) m_exec},
  {0, NULL}
};

static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_healpy_pixel_lib",
  docstring, 0, NULL, moduledef_slots,
  NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit__healpy_pixel_lib(void)
{
  import_array();
  import_ufunc();
  return PyModuleDef_Init(&moduledef);
}

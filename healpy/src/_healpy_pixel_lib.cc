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
#include "numpy/noprefix.h"

/*
   ang2pix
*/
template<Healpix_Ordering_Scheme scheme>static void
  ufunc_ang2pix(char **args, intp *dimensions, intp *steps, void *func)
{
  intp n=dimensions[0];

  intp is1=steps[0],is2=steps[1],is3=steps[2], os=steps[3];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(intp i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = pointing(*(double *)ip2,*(double *)ip3);
        ptg.normalize();
        *(long *)op = hb.ang2pix(ptg);
      } catch(PlanckError &e) {
        *(long *)op = -1;
      }
  }
}

/*
   pix2ang
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2ang(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = hb.pix2ang(*(long *)ip2);
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
  ufunc_xyf2pix(char **args, intp *dimensions, intp *steps, void *func)
{
  intp n=dimensions[0];

  intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3], os=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op=args[4];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(intp i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op+=os)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        *(long *)op = hb.xyf2pix((int)*(long *)ip2,(int)*(long *)ip3,(int)*(long *)ip4);
      } catch(PlanckError &e) {
        *(long *)op = -1;
      }
  }
}

/*
   pix2xyf
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2xyf(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        int x, y, f;
        hb.pix2xyf(*(long *)ip2, x, y, f);
        *(long *)op1 = x;
        *(long *)op2 = y;
        *(long *)op3 = f;
      } catch (PlanckError & e) {
        *(long *)op1 = -1;
        *(long *)op2 = -1;
        *(long *)op3 = -1;
      }
    }
}

/*
   ring2nest
*/
static void
ufunc_ring2nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os=steps[2];
  char *ip1=args[0], *ip2=args[1], *op=args[2];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, RING); }
      try {
        *(long *)op = hb.ring2nest(*(long *)ip2);
      } catch(PlanckError & e) {
        *(long *)op = -1;
      }
    }
}

/*
   nest2ring
*/
static void
 ufunc_nest2ring (char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os=steps[2];
  char *ip1=args[0], *ip2=args[1], *op=args[2];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, NEST); }
      try {
        *(long *)op = hb.nest2ring(*(long *)ip2);
      } catch(PlanckError & e) {
        *(long *)op = -1;
      }
    }
}


/*
  pix2vec
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2vec(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        vec3 v = hb.pix2vec(*(long *)ip2);
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
  ufunc_vec2pix(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      vec3 v (*(double *)ip2,*(double *)ip3,*(double *)ip4);
      try {
        long ipix = hb.vec2pix(v);
        *(long *)op1 = ipix;
      } catch (PlanckError &e) {
        *(long *)op1 = -1;
      }
    }
}


/*
  get_interpol
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_get_interpol(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],
    os1=steps[3],os2=steps[4],os3=steps[5],os4=steps[6],
    os5=steps[7],os6=steps[8],os7=steps[9],os8=steps[10];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2],
    *op1=args[3],*op2=args[4],*op3=args[5],*op4=args[6],
    *op5=args[7],*op6=args[8],*op7=args[9],*op8=args[10];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3,
        op1+=os1,op2+=os2,op3+=os3,op4+=os4,
        op5+=os5,op6+=os6,op7+=os7,op8+=os8 )
    {
      fix_arr<int64,4> pix;
      fix_arr<double,4> wgt;
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      try {
        pointing ptg = pointing(*(double*)ip2, *(double*)ip3);
        ptg.normalize();
        hb.get_interpol(ptg, pix, wgt);
        *(long*)op1 = (long)pix[0];
        *(long*)op2 = (long)pix[1];
        *(long*)op3 = (long)pix[2];
        *(long*)op4 = (long)pix[3];
        *(double*)op5 = wgt[0];
        *(double*)op6 = wgt[1];
        *(double*)op7 = wgt[2];
        *(double*)op8 = wgt[3];
      } catch (PlanckError &e) {
        *(long*)op1 = -1;
        *(long*)op2 = -1;
        *(long*)op3 = -1;
        *(long*)op4 = -1;
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
  ufunc_get_neighbors(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],
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
      hb.SetNside(*(long*)ip1, scheme);
      try {
        hb.neighbors(*(long*)ip2, pix);
        *(long*)op1 = (long)pix[0];
        *(long*)op2 = (long)pix[1];
        *(long*)op3 = (long)pix[2];
        *(long*)op4 = (long)pix[3];
        *(long*)op5 = (long)pix[4];
        *(long*)op6 = (long)pix[5];
        *(long*)op7 = (long)pix[6];
        *(long*)op8 = (long)pix[7];
      } catch (PlanckError & e) {
        *(long*)op1 = -1;
        *(long*)op2 = -1;
        *(long*)op3 = -1;
        *(long*)op4 = -1;
        *(long*)op5 = -1;
        *(long*)op6 = -1;
        *(long*)op7 = -1;
        *(long*)op8 = -1;
      }
    }
}

/*
  max_pixrad
*/
static void
  ufunc_max_pixrad(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],os1=steps[1];
  char *ip1=args[0], *op1=args[1];

  Healpix_Base2 hb;
  long oldnside=-1;

  for(i=0; i<n; i++, ip1+=is1, op1+=os1)
    {
      long nside = *(long*)ip1;
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, NEST);
        	/* ring and nest should give the same result */ 
        }
      double max_pixrad = hb.max_pixrad();
      *(double *)op1 = max_pixrad;
    }
}


static char *docstring = CP_(
  "This module contains basic ufunc related to healpix pixelisation\n"
  "scheme, such as ang2pix, ring<->nest swapping, etc.\n"
  "\n"
  "Available ufunc: _ang2pix_ring, _ang2pix_nest, _pix2ang_ring,\n"
  "                 _pix2ang_nest, _ring2nest, _nest2ring,\n"
  "                 _get_interpol_ring, _get_interpol_nest.");

/* to define the ufunc */
static PyUFuncGenericFunction ang2pix_ring_functions[] = {
  ufunc_ang2pix<RING>
};
static PyUFuncGenericFunction ang2pix_nest_functions[] = {
  ufunc_ang2pix<NEST>
};
static PyUFuncGenericFunction pix2ang_ring_functions[] = {
  ufunc_pix2ang<RING>
};
static PyUFuncGenericFunction pix2ang_nest_functions[] = {
  ufunc_pix2ang<NEST>
};
static PyUFuncGenericFunction xyf2pix_ring_functions[] = {
  ufunc_xyf2pix<RING>
};
static PyUFuncGenericFunction xyf2pix_nest_functions[] = {
  ufunc_xyf2pix<NEST>
};
static PyUFuncGenericFunction pix2xyf_ring_functions[] = {
  ufunc_pix2xyf<RING>
};
static PyUFuncGenericFunction pix2xyf_nest_functions[] = {
  ufunc_pix2xyf<NEST>
};
static PyUFuncGenericFunction vec2pix_ring_functions[] = {
  ufunc_vec2pix<RING>
};
static PyUFuncGenericFunction vec2pix_nest_functions[] = {
  ufunc_vec2pix<NEST>
};
static PyUFuncGenericFunction pix2vec_ring_functions[] = {
  ufunc_pix2vec<RING>
};
static PyUFuncGenericFunction pix2vec_nest_functions[] = {
  ufunc_pix2vec<NEST>
};
static PyUFuncGenericFunction ring2nest_functions[] = {
  ufunc_ring2nest
};
static PyUFuncGenericFunction nest2ring_functions[] = {
  ufunc_nest2ring
};
static PyUFuncGenericFunction get_interpol_ring_functions[] = {
  ufunc_get_interpol<RING>
};
static PyUFuncGenericFunction get_interpol_nest_functions[] = {
  ufunc_get_interpol<NEST>
};
static PyUFuncGenericFunction get_neighbors_ring_functions[] = {
  ufunc_get_neighbors<RING>
};
static PyUFuncGenericFunction get_neighbors_nest_functions[] = {
  ufunc_get_neighbors<NEST>
};
static PyUFuncGenericFunction max_pixrad_functions[] = {
  ufunc_max_pixrad
};


static void * blank_data[] = { (void *)NULL };

static char ang2pix_signatures[] = {
  PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_LONG
};
static char pix2ang_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE
};
static char xyf2pix_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG
};
static char pix2xyf_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG
};
static char pix2vec_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE
};
static char vec2pix_signatures[] = {
  PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_LONG
};
static char ring2nest_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_LONG
};
static char get_interpol_signatures[] = {
  PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE,
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG,
  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE
};
static char get_neighbors_ring_signatures[] = {
  PyArray_LONG, PyArray_LONG, // input
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG, // output
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG // output
};
static char get_neighbors_nest_signatures[] = {
  PyArray_LONG, PyArray_LONG, // input
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG, // output
  PyArray_LONG, PyArray_LONG, PyArray_LONG, PyArray_LONG // output
};
static char max_pixrad_signatures[] = {
  PyArray_LONG, PyArray_DOUBLE  
};

#if PY_MAJOR_VERSION >= 3
static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_healpy_pixel_lib",
  NULL, -1, NULL
};
#endif

#if PY_MAJOR_VERSION < 3
#define FREE_MODULE_AND_FAIL do { return; } while(0)
#else
#define FREE_MODULE_AND_FAIL do { Py_DECREF(m); return NULL; } while(0)
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION < 3
init_healpy_pixel_lib(void)
#else
PyInit__healpy_pixel_lib(void)
#endif
{
  PyObject *m;

  import_array();
  import_ufunc();

#if PY_MAJOR_VERSION < 3
  m = Py_InitModule3("_healpy_pixel_lib", NULL, docstring);
	if (!m) return;
#else
	m = PyModule_Create(&moduledef);
	if (!m) return NULL;
#endif

  if (PyModule_AddObject(m, "_ang2pix_ring", PyUFunc_FromFuncAndData(
      ang2pix_ring_functions, blank_data,
      ang2pix_signatures, 1,
      3, 1, PyUFunc_None, CP_("_ang2pix_ring"),
      CP_("nside,theta,phi [rad] -> ipix (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_ang2pix_nest", PyUFunc_FromFuncAndData(
      ang2pix_nest_functions, blank_data,
      ang2pix_signatures, 1,
      3, 1, PyUFunc_None, CP_("_ang2pix_nest"),
      CP_("nside,theta,phi [rad] -> ipix (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2ang_ring", PyUFunc_FromFuncAndData(
      pix2ang_ring_functions, blank_data,
      pix2ang_signatures, 1,
      2, 2, PyUFunc_None, CP_("_pix2ang_ring"),
      CP_("nside,ipix -> theta,phi [rad] (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2ang_nest", PyUFunc_FromFuncAndData(
      pix2ang_nest_functions, blank_data,
      pix2ang_signatures, 1,
      2, 2, PyUFunc_None, CP_("_pix2ang_nest"),
      CP_("nside,ipix -> theta,phi [rad] (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  //=========

  if (PyModule_AddObject(m, "_xyf2pix_ring", PyUFunc_FromFuncAndData(
      xyf2pix_ring_functions, blank_data,
      xyf2pix_signatures, 1,
      4, 1, PyUFunc_None, CP_("_xyf2pix_ring"),
      CP_("nside,x,y,face -> ipix (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_xyf2pix_nest", PyUFunc_FromFuncAndData(
      xyf2pix_nest_functions, blank_data,
      xyf2pix_signatures, 1,
      4, 1, PyUFunc_None, CP_("_xyf2pix_nest"),
      CP_("nside,x,y,face -> ipix (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2xyf_ring", PyUFunc_FromFuncAndData(
      pix2xyf_ring_functions, blank_data,
      pix2xyf_signatures, 1,
      2, 3, PyUFunc_None, CP_("_pix2xyf_ring"),
      CP_("nside,ipix -> x,y,face (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2xyf_nest", PyUFunc_FromFuncAndData(
      pix2xyf_nest_functions, blank_data,
      pix2xyf_signatures, 1,
      2, 3, PyUFunc_None, CP_("_pix2xyf_nest"),
      CP_("nside,ipix -> x,y,face (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  //=========

  if (PyModule_AddObject(m, "_vec2pix_ring", PyUFunc_FromFuncAndData(
      vec2pix_ring_functions, blank_data,
      vec2pix_signatures, 1,
      4, 1, PyUFunc_None, CP_("_vec2pix_ring"),
      CP_("nside,x,y,z -> ipix (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_vec2pix_nest", PyUFunc_FromFuncAndData(
      vec2pix_nest_functions, blank_data,
      vec2pix_signatures, 1,
      4, 1, PyUFunc_None, CP_("_vec2pix_nest"),
      CP_("nside,x,y,z -> ipix (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2vec_ring", PyUFunc_FromFuncAndData(
      pix2vec_ring_functions, blank_data,
      pix2vec_signatures, 1,
      2, 3, PyUFunc_None, CP_("_pix2vec_ring"),
      CP_("nside,ipix -> x,y,z (RING)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_pix2vec_nest", PyUFunc_FromFuncAndData(
      pix2vec_nest_functions, blank_data,
      pix2vec_signatures, 1,
      2, 3, PyUFunc_None, CP_("_pix2vec_nest"),
      CP_("nside,ipix -> x,y,z (NEST)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  //=============

  if (PyModule_AddObject(m, "_ring2nest", PyUFunc_FromFuncAndData(
      ring2nest_functions, blank_data,
      ring2nest_signatures, 1,
      2, 1, PyUFunc_None, CP_("_ring2nest"),
      CP_("ipix(ring) -> ipix(nest)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_nest2ring", PyUFunc_FromFuncAndData(
      nest2ring_functions, blank_data,
      ring2nest_signatures, 1,
      2, 1, PyUFunc_None, CP_("_nest2ring"),
      CP_("ipix(nest) -> ipix(ring)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_get_interpol_ring", PyUFunc_FromFuncAndData(
      get_interpol_ring_functions, blank_data,
      get_interpol_signatures, 1,
      3, 8, PyUFunc_None, CP_("_get_interpol_ring"),
      CP_("nside,theta,phi->4 nearest pixels+4weights"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_get_interpol_nest", PyUFunc_FromFuncAndData(
      get_interpol_nest_functions, blank_data,
      get_interpol_signatures, 1,
      3, 8, PyUFunc_None, CP_("_get_interpol_nest"),
      CP_("nside,theta,phi->4 nearest pixels+4weights"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_get_neighbors_ring", PyUFunc_FromFuncAndData(
      get_neighbors_ring_functions, blank_data,
      get_neighbors_ring_signatures, 1,
      2, 8, PyUFunc_None, CP_("_get_neigbors_ring"),
      CP_("nside, ipix [rad] -> 8 neighbors"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "_get_neighbors_nest", PyUFunc_FromFuncAndData(
      get_neighbors_nest_functions, blank_data,
      get_neighbors_nest_signatures, 1,
      2, 8, PyUFunc_None, CP_("_get_neigbors_nest"),
      CP_("nside, ipix [rad] -> 8 neighbors"),0)) < 0)
    FREE_MODULE_AND_FAIL;


  if (PyModule_AddObject(m, "_max_pixrad", PyUFunc_FromFuncAndData(
      max_pixrad_functions, blank_data,
      max_pixrad_signatures, 1,
      1, 1, PyUFunc_None, CP_("max_pixrad"),
      CP_("nside -> max_distance to pixel corners from center)"),0)) < 0)
    FREE_MODULE_AND_FAIL;

  if (PyModule_AddObject(m, "UNSEEN", PyFloat_FromDouble(Healpix_undef)) < 0)
    FREE_MODULE_AND_FAIL;

#if PY_MAJOR_VERSION >= 3
  return m;
#endif
}

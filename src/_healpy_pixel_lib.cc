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
   OpenMP

   Each element of these ufunc loops is independent, so loop() splits long
   inputs across OpenMP threads (OMP_NUM_THREADS controls how many).  Every
   thread keeps its own Healpix_Base2 and nside cache.  Inputs shorter than
   HEALPY_OMP_MIN_N run on the calling thread: for them starting the team
   costs more than the work.

   GNU libgomp cannot start a team in a child created by fork() once the
   parent has used one (the child waits forever for the parent's pool), and
   multiprocessing forks by default on Linux before Python 3.14.  The loops
   therefore only use threads in the process that imported this module; a
   forked child runs them serially.

   When healpy is built without OpenMP the loops are serial.
*/
#ifdef _OPENMP
#include <omp.h>
#ifndef _WIN32
#include <unistd.h>
static pid_t import_pid;
#endif
#endif

#define HEALPY_OMP_MIN_N 1000

static bool use_threads(npy_intp n)
{
#ifdef _OPENMP
  if (n < HEALPY_OMP_MIN_N || omp_get_max_threads() < 2)
    return false;
#ifndef _WIN32
  if (getpid() != import_pid)
    return false;
#endif
  return true;
#else
  (void)n;
  return false;
#endif
}

/*
   Call body(hb, i) for i in [0, n), with hb set to the nside of element i
   (read from nside_ptr with stride nside_step).
*/
template<Healpix_Ordering_Scheme scheme, typename Body> static void
  loop(npy_intp n, const char *nside_ptr, npy_intp nside_step, Body body)
{
#ifdef _OPENMP
  if (use_threads(n))
    {
#pragma omp parallel
      {
        Healpix_Base2 hb;
        int64 oldnside=-1;
#pragma omp for schedule(static)
        for(npy_intp i=0; i<n; i++)
          {
            int64 nside = *(const int64 *)(nside_ptr+i*nside_step);
            if (nside!=oldnside)
              { oldnside=nside; hb.SetNside(nside, scheme); }
            body(hb, i);
          }
      }
      return;
    }
#endif
  Healpix_Base2 hb;
  int64 oldnside=-1;
  for(npy_intp i=0; i<n; i++)
    {
      int64 nside = *(const int64 *)(nside_ptr+i*nside_step);
      if (nside!=oldnside)
        { oldnside=nside; hb.SetNside(nside, scheme); }
      body(hb, i);
    }
}

/*
   ang2pix
*/
template<Healpix_Ordering_Scheme scheme>static void
  ufunc_ang2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],is3=steps[2], os=steps[3];
  char *ip2=args[1], *ip3=args[2], *op=args[3];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        pointing ptg = pointing(*(double *)(ip2+i*is2),*(double *)(ip3+i*is3));
        ptg.normalize();
        *(int64 *)(op+i*os) = hb.ang2pix(ptg);
      } catch(PlanckError &e) {
        *(int64 *)(op+i*os) = -1;
      }
    });
}

/*
   pix2ang
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2ang(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],os1=steps[2],os2=steps[3];
  char *ip2=args[1], *op1=args[2], *op2=args[3];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        pointing ptg = hb.pix2ang(*(int64 *)(ip2+i*is2));
        *(double *)(op1+i*os1) = ptg.theta;
        *(double *)(op2+i*os2) = ptg.phi;
      } catch (PlanckError & e) {
        *(double *)(op1+i*os1) = NAN;
        *(double *)(op2+i*os2) = NAN;
      }
    });
}

/*
   xyf2pix
*/
template<Healpix_Ordering_Scheme scheme>static void
  ufunc_xyf2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],is3=steps[2],is4=steps[3], os=steps[4];
  char *ip2=args[1], *ip3=args[2], *ip4=args[3], *op=args[4];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        *(int64 *)(op+i*os) = hb.xyf2pix((int)*(int64 *)(ip2+i*is2),
          (int)*(int64 *)(ip3+i*is3),(int)*(int64 *)(ip4+i*is4));
      } catch(PlanckError &e) {
        *(int64 *)(op+i*os) = -1;
      }
    });
}

/*
   pix2xyf
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2xyf(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        int x, y, f;
        hb.pix2xyf(*(int64 *)(ip2+i*is2), x, y, f);
        *(int64 *)(op1+i*os1) = x;
        *(int64 *)(op2+i*os2) = y;
        *(int64 *)(op3+i*os3) = f;
      } catch (PlanckError & e) {
        *(int64 *)(op1+i*os1) = -1;
        *(int64 *)(op2+i*os2) = -1;
        *(int64 *)(op3+i*os3) = -1;
      }
    });
}

/*
   ring2nest
*/
static void
ufunc_ring2nest(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],os=steps[2];
  char *ip2=args[1], *op=args[2];

  loop<RING>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        *(int64 *)(op+i*os) = hb.ring2nest(*(int64 *)(ip2+i*is2));
      } catch(PlanckError & e) {
        *(int64 *)(op+i*os) = -1;
      }
    });
}

/*
   nest2ring
*/
static void
 ufunc_nest2ring (char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],os=steps[2];
  char *ip2=args[1], *op=args[2];

  loop<NEST>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        *(int64 *)(op+i*os) = hb.nest2ring(*(int64 *)(ip2+i*is2));
      } catch(PlanckError & e) {
        *(int64 *)(op+i*os) = -1;
      }
    });
}


/*
  pix2vec
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_pix2vec(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      try {
        vec3 v = hb.pix2vec(*(int64 *)(ip2+i*is2));
        *(double *)(op1+i*os1) = v.x;
        *(double *)(op2+i*os2) = v.y;
        *(double *)(op3+i*os3) = v.z;
      } catch (PlanckError & e) {
        *(double *)(op1+i*os1) = NAN;
        *(double *)(op2+i*os2) = NAN;
        *(double *)(op3+i*os3) = NAN;
      }
    });
}

/*
  vec2pix
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_vec2pix(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4];
  char *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      vec3 v (*(double *)(ip2+i*is2),*(double *)(ip3+i*is3),*(double *)(ip4+i*is4));
      try {
        int64 ipix = hb.vec2pix(v);
        *(int64 *)(op1+i*os1) = ipix;
      } catch (PlanckError &e) {
        *(int64 *)(op1+i*os1) = -1;
      }
    });
}


/*
  get_interpol
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_get_interpol(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],is3=steps[2],
    os1=steps[3],os2=steps[4],os3=steps[5],os4=steps[6],
    os5=steps[7],os6=steps[8],os7=steps[9],os8=steps[10];
  char *ip2=args[1], *ip3=args[2],
    *op1=args[3],*op2=args[4],*op3=args[5],*op4=args[6],
    *op5=args[7],*op6=args[8],*op7=args[9],*op8=args[10];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      fix_arr<int64,4> pix;
      fix_arr<double,4> wgt;
      try {
        pointing ptg = pointing(*(double*)(ip2+i*is2), *(double*)(ip3+i*is3));
        ptg.normalize();
        hb.get_interpol(ptg, pix, wgt);
        *(int64*)(op1+i*os1) = (int64)pix[0];
        *(int64*)(op2+i*os2) = (int64)pix[1];
        *(int64*)(op3+i*os3) = (int64)pix[2];
        *(int64*)(op4+i*os4) = (int64)pix[3];
        *(double*)(op5+i*os5) = wgt[0];
        *(double*)(op6+i*os6) = wgt[1];
        *(double*)(op7+i*os7) = wgt[2];
        *(double*)(op8+i*os8) = wgt[3];
      } catch (PlanckError &e) {
        *(int64*)(op1+i*os1) = -1;
        *(int64*)(op2+i*os2) = -1;
        *(int64*)(op3+i*os3) = -1;
        *(int64*)(op4+i*os4) = -1;
        *(double*)(op5+i*os5) = NAN;
        *(double*)(op6+i*os6) = NAN;
        *(double*)(op7+i*os7) = NAN;
        *(double*)(op8+i*os8) = NAN;
      }
    });
}


/*
  get_neighbors
*/
template<Healpix_Ordering_Scheme scheme> static void
  ufunc_get_neighbors(char **args, const npy_intp *dimensions, const npy_intp *steps, void *func)
{
  npy_intp is2=steps[1],
    os1=steps[2],os2=steps[3],os3=steps[4],os4=steps[5],
    os5=steps[6],os6=steps[7],os7=steps[8],os8=steps[9];
  char *ip2=args[1],
    *op1=args[2],*op2=args[3],*op3=args[4],*op4=args[5],
    *op5=args[6],*op6=args[7],*op7=args[8],*op8=args[9];

  loop<scheme>(dimensions[0], args[0], steps[0], [&](const Healpix_Base2 &hb, npy_intp i)
    {
      fix_arr<int64,8> pix;
      try {
        hb.neighbors(*(int64*)(ip2+i*is2), pix);
        *(int64*)(op1+i*os1) = (int64)pix[0];
        *(int64*)(op2+i*os2) = (int64)pix[1];
        *(int64*)(op3+i*os3) = (int64)pix[2];
        *(int64*)(op4+i*os4) = (int64)pix[3];
        *(int64*)(op5+i*os5) = (int64)pix[4];
        *(int64*)(op6+i*os6) = (int64)pix[5];
        *(int64*)(op7+i*os7) = (int64)pix[6];
        *(int64*)(op8+i*os8) = (int64)pix[7];
      } catch (PlanckError & e) {
        *(int64*)(op1+i*os1) = -1;
        *(int64*)(op2+i*os2) = -1;
        *(int64*)(op3+i*os3) = -1;
        *(int64*)(op4+i*os4) = -1;
        *(int64*)(op5+i*os5) = -1;
        *(int64*)(op6+i*os6) = -1;
        *(int64*)(op7+i*os7) = -1;
        *(int64*)(op8+i*os8) = -1;
      }
    });
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


static char docstring[] =
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
#if defined(_OPENMP) && !defined(_WIN32)
  import_pid = getpid();
#endif

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

static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_healpy_pixel_lib",
  docstring, 0, NULL, (PyModuleDef_Slot []) {
        {Py_mod_exec, (void *) m_exec},
        {/* terminal element, all NULL */}
  }
};

PyMODINIT_FUNC
PyInit__healpy_pixel_lib(void)
{
  import_array();
  import_ufunc();
  return PyModuleDef_Init(&moduledef);
}

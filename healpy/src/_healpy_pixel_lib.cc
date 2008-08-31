/*

   This module provides Healpix functions to Python.
   It uses the healpix_cxx library.

*/

#include <Python.h>

#include "arr.h"
#include "healpix_base.h"

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "numpy/noprefix.h"

static void
ufunc_ang2pix_ring(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_ang2pix_nest(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_pix2ang_ring(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_pix2ang_nest(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_ring2nest(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_nest2ring(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_get_interpol_ring(char **args, intp *dimension, intp *steps, void *funv);
static void
ufunc_get_interpol_nest(char **args, intp *dimension, intp *steps, void *funv);


static void
ufunc_pix2vec_ring(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_pix2vec_nest(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_vec2pix_ring(char **args, intp *dimensions, intp *steps, void *func);
static void
ufunc_vec2pix_nest(char **args, intp *dimensions, intp *steps, void *func);


static char *docstring = 
  "This module conatains basic ufunc related to healpix pixellisation\n"
  "scheme, such as ang2pix, ring<->nest swapping, etc.\n"
  "\n"
  "Available ufunc: _ang2pix_ring, _ang2pix_nest, _pix2ang_ring,\n"
  "                 _pix2ang_nest, _ring2nest, _nest2ring,\n"
  "                 _get_interpol_ring, _get_interpol_nest.";

/* to define the ufunc */
static PyUFuncGenericFunction ang2pix_ring_functions[] = {
  ufunc_ang2pix_ring
};
static PyUFuncGenericFunction ang2pix_nest_functions[] = {
  ufunc_ang2pix_nest
};
static PyUFuncGenericFunction pix2ang_ring_functions[] = {
  ufunc_pix2ang_ring
};
static PyUFuncGenericFunction pix2ang_nest_functions[] = {
  ufunc_pix2ang_nest
};
static PyUFuncGenericFunction vec2pix_ring_functions[] = {
  ufunc_vec2pix_ring
};
static PyUFuncGenericFunction vec2pix_nest_functions[] = {
  ufunc_vec2pix_nest
};
static PyUFuncGenericFunction pix2vec_ring_functions[] = {
  ufunc_pix2vec_ring
};
static PyUFuncGenericFunction pix2vec_nest_functions[] = {
  ufunc_pix2vec_nest
};
static PyUFuncGenericFunction ring2nest_functions[] = {
  ufunc_ring2nest
};
static PyUFuncGenericFunction nest2ring_functions[] = {
  ufunc_nest2ring
};
static PyUFuncGenericFunction get_interpol_ring_functions[] = {
  ufunc_get_interpol_ring
};
static PyUFuncGenericFunction get_interpol_nest_functions[] = {
  ufunc_get_interpol_nest
};


static void * blank_data[] = { (void *)NULL };

static char ang2pix_signatures[] = {
  PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_LONG
};
static char pix2ang_signatures[] = {
  PyArray_LONG, PyArray_LONG, PyArray_DOUBLE, PyArray_DOUBLE
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

PyMODINIT_FUNC
init_healpy_pixel_lib(void)
{
  PyObject *m, *d, *f;

  m = Py_InitModule3("_healpy_pixel_lib", NULL, docstring);
  import_array();
  import_ufunc();

  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);
  
  f = PyUFunc_FromFuncAndData(ang2pix_ring_functions, blank_data,
			      ang2pix_signatures, 1,
			      3, 1, PyUFunc_None, "_ang2pix_ring",
			      "nside,theta,phi -> ipix (RING)",0);

  PyDict_SetItemString(d, "_ang2pix_ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(ang2pix_nest_functions, blank_data,
			      ang2pix_signatures, 1,
			      3, 1, PyUFunc_None, "_ang2pix_nest",
			      "nside,theta,phi -> ipix (NEST)",0);

  PyDict_SetItemString(d, "_ang2pix_nest", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(pix2ang_ring_functions, blank_data,
			      pix2ang_signatures, 1,
			      2, 2, PyUFunc_None, "_pix2ang_ring",
			      "nside,ipix -> theta,phi (RING)",0);

  PyDict_SetItemString(d, "_pix2ang_ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(pix2ang_nest_functions, blank_data,
			      pix2ang_signatures, 1,
			      2, 2, PyUFunc_None, "_pix2ang_nest",
			      "nside,ipix -> theta,phi (NEST)",0);

  PyDict_SetItemString(d, "_pix2ang_nest", f);
  Py_DECREF(f);
  
  //=========

  f = PyUFunc_FromFuncAndData(vec2pix_ring_functions, blank_data,
			      vec2pix_signatures, 1,
			      4, 1, PyUFunc_None, "_vec2pix_ring",
			      "nside,x,y,z -> ipix (RING)",0);

  PyDict_SetItemString(d, "_vec2pix_ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(vec2pix_nest_functions, blank_data,
			      vec2pix_signatures, 1,
			      4, 1, PyUFunc_None, "_vec2pix_nest",
			      "nside,x,y,z -> ipix (NEST)",0);

  PyDict_SetItemString(d, "_vec2pix_nest", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(pix2vec_ring_functions, blank_data,
			      pix2vec_signatures, 1,
			      2, 3, PyUFunc_None, "_pix2vec_ring",
			      "nside,ipix -> x,y,z (RING)",0);

  PyDict_SetItemString(d, "_pix2vec_ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(pix2vec_nest_functions, blank_data,
			      pix2vec_signatures, 1,
			      2, 3, PyUFunc_None, "_pix2vec_nest",
			      "nside,ipix -> x,y,z (NEST)",0);

  PyDict_SetItemString(d, "_pix2vec_nest", f);
  Py_DECREF(f);
  
  //=============
  
  f = PyUFunc_FromFuncAndData(ring2nest_functions, blank_data,
			      ring2nest_signatures, 1,
			      2, 1, PyUFunc_None, "_ring2nest",
			      "ipix(ring) -> ipix(nest)",0);

  PyDict_SetItemString(d, "_ring2nest", f);
  Py_DECREF(f);

  f = PyUFunc_FromFuncAndData(nest2ring_functions, blank_data,
			      ring2nest_signatures, 1,
			      2, 1, PyUFunc_None, "_nest2ring",
			      "ipix(nest) -> ipix(ring)",0);

  PyDict_SetItemString(d, "_nest2ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(get_interpol_ring_functions, blank_data,
			      get_interpol_signatures, 1,
			      3, 8, PyUFunc_None, "_get_interpol_ring",
			      "nside,theta,phi->4 nearest pixels+4weights",0);

  PyDict_SetItemString(d, "_get_interpol_ring", f);
  Py_DECREF(f);
  
  f = PyUFunc_FromFuncAndData(get_interpol_nest_functions, blank_data,
			      get_interpol_signatures, 1,
			      3, 8, PyUFunc_None, "_get_interpol_nest",
			      "nside,theta,phi->4 nearest pixels+4weights",0);

  PyDict_SetItemString(d, "_get_interpol_nest", f);
  Py_DECREF(f);

  
  f = PyFloat_FromDouble(Healpix_undef);

  PyDict_SetItemString(d, "UNSEEN", f);
  Py_DECREF(f);
  
  return;
}

/* 
   ang2pix_ring
*/
static void
ufunc_ang2pix_ring(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],
    os=steps[3];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) 
    {
      hb.SetNside(*(long*)ip1,RING);
      *(long *)op = hb.ang2pix(pointing(*(double *)ip2,
					*(double *)ip3));
    }
}

/* 
   ang2pix_nest
*/
static void
ufunc_ang2pix_nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],
    os=steps[3];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) 
    {
      hb.SetNside(*(long*)ip1,NEST);
      *(long *)op = hb.ang2pix(pointing(*(double *)ip2,
					*(double *)ip3));
    }
}

/*
   pix2ang_ring
*/
static void
ufunc_pix2ang_ring(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2) 
    {
      pointing ptg;
      hb.SetNside(*(long*)ip1,RING);
      ptg = hb.pix2ang(*(long *)ip2);
      *(double *)op1 = ptg.theta;
      *(double *)op2 = ptg.phi;
    }
}

/*
   pix2ang_nest
*/
static void
ufunc_pix2ang_nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2) 
    {
      pointing ptg;
      hb.SetNside(*(long*)ip1,NEST);
      ptg = hb.pix2ang(*(long *)ip2);
      *(double *)op1 = ptg.theta;
      *(double *)op2 = ptg.phi;
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
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) 
    {
      hb.SetNside(*(long*)ip1,RING);
      *(long *)op = hb.ring2nest(*(long *)ip2);
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
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) 
    {
      hb.SetNside(*(long*)ip1,NEST);
      *(long *)op = hb.nest2ring(*(long *)ip2);
    }
}

/*
  pix2vec_ring
*/
static void
ufunc_pix2vec_ring(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3) 
    {
      vec3 v;
      hb.SetNside(*(long*)ip1,RING);
      v = hb.pix2vec(*(long *)ip2);
      *(double *)op1 = v.x;
      *(double *)op2 = v.y;
      *(double *)op3 = v.z;
    }
}


/*
  pix2vec_nest
*/
static void
ufunc_pix2vec_nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4];
  char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3) 
    {
      vec3 v;
      hb.SetNside(*(long*)ip1,NEST);
      v = hb.pix2vec(*(long *)ip2);
      *(double *)op1 = v.x;
      *(double *)op2 = v.y;
      *(double *)op3 = v.z;
    }
}

/*
  vec2pix_ring
*/
static void
ufunc_vec2pix_ring(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1) 
    {
      vec3 v;
      long ipix;
      hb.SetNside(*(long*)ip1,RING);
      v.x = *(double *)ip2;
      v.y = *(double *)ip3;
      v.z = *(double *)ip4;
      ipix = hb.vec2pix(v);
      *(long *)op1 = ipix;
    }
}

/*
  vec2pix_nest
*/
static void
ufunc_vec2pix_nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1) 
    {
      vec3 v;
      long ipix;
      hb.SetNside(*(long*)ip1,NEST);
      v.x = *(double *)ip2;
      v.y = *(double *)ip3;
      v.z = *(double *)ip4;
      ipix = hb.vec2pix(v);
      *(long *)op1 = ipix;
    }
}


/*
  get_interpol_ring
*/
static void
ufunc_get_interpol_ring(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],
    os1=steps[3],os2=steps[4],os3=steps[5],os4=steps[6],
    os5=steps[7],os6=steps[8],os7=steps[9],os8=steps[10];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2],
    *op1=args[3],*op2=args[4],*op3=args[5],*op4=args[6],
    *op5=args[7],*op6=args[8],*op7=args[9],*op8=args[10];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3,
	op1+=os1,op2+=os2,op3+=os3,op4+=os4,
	op5+=os5,op6+=os6,op7+=os7,op8+=os8 ) 
    {
      fix_arr<int,4> pix;
      fix_arr<double,4> wgt;
      hb.SetNside(*(long*)ip1,RING);
      hb.get_interpol(pointing(*(double*)ip2, *(double*)ip3),
		       pix, wgt);
      *(long*)op1 = (long)pix[0];
      *(long*)op2 = (long)pix[1];
      *(long*)op3 = (long)pix[2];
      *(long*)op4 = (long)pix[3];
      *(double*)op5 = wgt[0];
      *(double*)op6 = wgt[1];
      *(double*)op7 = wgt[2];
      *(double*)op8 = wgt[3];
    }
}

/*
  get_interpol_nest
*/
static void
ufunc_get_interpol_nest(char **args, intp *dimensions, intp *steps, void *func)
{
  register intp i, n=dimensions[0];
  register intp is1=steps[0],is2=steps[1],is3=steps[2],
    os1=steps[3],os2=steps[4],os3=steps[5],os4=steps[6],
    os5=steps[7],os6=steps[8],os7=steps[9],os8=steps[10];
  char *ip1=args[0], *ip2=args[1], *ip3=args[2],
    *op1=args[3],*op2=args[4],*op3=args[5],*op4=args[6],
    *op5=args[7],*op6=args[8],*op7=args[9],*op8=args[10];
  
  Healpix_Base hb;
  
  for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3,
	op1+=os1,op2+=os2,op3+=os3,op4+=os4,
	op5+=os5,op6+=os6,op7+=os7,op8+=os8 ) 
    {
      fix_arr<int,4> pix;
      fix_arr<double,4> wgt;
      hb.SetNside(*(long*)ip1,NEST);
      hb.get_interpol(pointing(*(double*)ip2, *(double*)ip3),
		      pix, wgt);
      *(long*)op1 = (long)pix[0];
      *(long*)op2 = (long)pix[1];
      *(long*)op3 = (long)pix[2];
      *(long*)op4 = (long)pix[3];
      *(double*)op5 = wgt[0];
      *(double*)op6 = wgt[1];
      *(double*)op7 = wgt[2];
      *(double*)op8 = wgt[3];
    }
}

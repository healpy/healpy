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

#include "numpy/arrayobject.h"

#include <string>
//#include <string.h>
#include <iostream>

#include "arr.h"
#include "alm.h"
#include "healpix_map.h"
#include "xcomplex.h"
#include "alm_healpix_tools.h"
#include "alm_map_tools.h"
#include "powspec.h"
#include "alm_powspec_tools.h"
#include "ylmgen.h"
#include "healpix_data_io.h"

#define IS_DEBUG_ON 0

/* Some helpful macro */
#define XMALLOC(X,Y,Z) if( !(X = (Y*)malloc(Z*sizeof(Y))) ) { PyErr_NoMemory(); goto fail;}
#define XNEW(X,Y,Z) if( !(X = new(Y[Z])) ) { PyErr_NoMemory(); goto fail; }
#define XFREE(X) if( X ) free(X);
#define XDELETE(X) if( X ) delete[] X;
#define DBGPRINTF(X,...) if( IS_DEBUG_ON ) printf(X, ## __VA_ARGS__)

static double sky_signal_direct(Alm<xcomplex<double> >alm, 
				double theta, double phi);

static long npix2nside(long npix);
static long nside2npix(long nside);
static long getn(long s);
//static void getij(long n, long idx, long *i, long *j);
static long getidx(long n, long i, long j);
static void cholesky(int n, double *data, double *res);


static PyObject *healpy_map2alm(PyObject *self, PyObject *args, 
				PyObject *kwds);

static PyObject *healpy_alm2map(PyObject *self, PyObject *args, 
				PyObject *kwds);

static PyObject *healpy_alm2map_der1(PyObject *self, PyObject *args, 
				PyObject *kwds);

static PyObject *healpy_synalm(PyObject *self, PyObject *args, 
			       PyObject *kwds);

static PyObject *healpy_getn(PyObject *self, PyObject *args);

static PyObject *healpy_alm2signal(PyObject *self, PyObject *args, 
				   PyObject *kwds);

static PyObject *healpy_Ylm(PyObject *self, PyObject *args, 
			    PyObject *kwds);

static PyMethodDef SphtMethods[] = {
  {"_map2alm", (PyCFunction)healpy_map2alm, METH_VARARGS | METH_KEYWORDS, 
   "Compute alm or cl from an input map.\n"
   "The input map is assumed to be ordered in RING.\n"
   "anafast(map,lmax=3*nside-1,mmax=lmax,cl=False,\n"
   "        iter=3,use_weights=False,data_path=None,regression=True)"},
  {"_alm2map", (PyCFunction)healpy_alm2map, METH_VARARGS | METH_KEYWORDS, 
   "Compute a map from alm.\n"
   "The output map is ordered in RING scheme.\n"
   "alm2map(alm,nside=64,lmax=-1,mmax=-1)"},
  {"_alm2map_der1", (PyCFunction)healpy_alm2map_der1, METH_VARARGS | METH_KEYWORDS, 
   "Compute a map and derivatives from alm.\n"
   "The output map is ordered in RING scheme.\n"
   "alm2map_der1(alm,nside=64,lmax=-1,mmax=-1)"},
  {"_synalm", (PyCFunction)healpy_synalm, METH_VARARGS | METH_KEYWORDS,
   "Compute alm's given cl's and unit variance random arrays.\n"},
  {"_getn", healpy_getn, METH_VARARGS,
   "Compute number n such that n(n+1)/2 is equal to the argument.\n"},
  {"_alm2signal", (PyCFunction)healpy_alm2signal, METH_VARARGS | METH_KEYWORDS,
   "Compute sum(Alm*Ylm) for a direction (theta,phi).\n"},
  {"_getylm", (PyCFunction)healpy_Ylm, METH_VARARGS | METH_KEYWORDS,
   "Compute Ylm(theta,0) for m and l=m..lmax. Return a lmax+1 long array\n"
   "with values for l<m set to zero.\n"
   "call: _getylm(lmax,mmax,m,theta)"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC
init_healpy_sph_transform_lib(void)
{
  PyObject *m;
  m =  Py_InitModule("_healpy_sph_transform_lib", SphtMethods);
  
  import_array();
}

/***********************************************************************
    map2alm

       input: map, lmax=3*nside-1, mmax=lmax, cl=False 
              iter=3, use_weights=False

       output: alm (or cl if cl=True)
*/
static PyObject *healpy_map2alm(PyObject *self, PyObject *args, 
				PyObject *kwds)
{
  PyArrayObject *mapIin = NULL;
  PyArrayObject *mapQin = NULL;
  PyArrayObject *mapUin = NULL;
  int lmax=-1, mmax=-1;
  int nside=-1;
  int npix=-1;
  int num_iter=3;
  int docl=0;
  int use_weights=0;
  char * datapath=NULL;
  int polarisation = 0; /* not polarised by default */
  int regression=1;

  static char* kwlist[] = {"","lmax", "mmax","cl","iter", 
			   "use_weights", "data_path", "regression", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iiiiisi", kwlist,
				   &PyArray_Type, &mapIin,
				   &lmax, &mmax, &docl,
				   &num_iter,&use_weights,&datapath,&regression))
    {
      PyErr_Clear(); /* I want to try the other calling way */

      PyObject *t = NULL;
      if( !PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiiisi", kwlist,
				   &t,
				   &lmax, &mmax, &docl,
				       &num_iter,&use_weights,&datapath,&regression) )	
	return NULL;
      else
	{
	  if( PySequence_Size(t) != 3 )
	    {
	      PyErr_SetString(PyExc_TypeError, 
			      "First argument must be a sequence with "
			      "three elements.");
	      return NULL;
	    }
	  PyObject *o1 = NULL;
	  PyObject *o2 = NULL;
	  PyObject *o3 = NULL;
	  o1 = PySequence_GetItem(t, 0);
	  o2 = PySequence_GetItem(t, 1);
	  o3 = PySequence_GetItem(t, 2);
	  /* I decrease reference counts here,
	     because PySequence_GetItem increase 
	     reference count, and I just want to 
	     borrow a reference the time of this 
	     function. */
	  Py_XDECREF(o1);
	  Py_XDECREF(o2);
	  Py_XDECREF(o3);
	  if( ! ( PyArray_Check(o1) && PyArray_Check(o2) 
		  && PyArray_Check(o3) ) )
	    {
	      PyErr_SetString(PyExc_TypeError, 
			      "First argument must be a sequence with "
			      "three arrays");
	      return NULL;
	    }
	  else
	    {
	      mapIin = (PyArrayObject*) o1;
	      mapQin = (PyArrayObject*) o2;
	      mapUin = (PyArrayObject*) o3;
	    }
	}
	polarisation = 1;  /* we have three maps : polarisation! */
    }

  /* Check array is contiguous */
  if( !(mapIin->flags & NPY_C_CONTIGUOUS) ) 
    {
      PyErr_SetString(PyExc_ValueError,
		      "Array must be C contiguous for this operation.");
      return NULL;      
    }

  if( polarisation )
    {
      if( !(mapQin->flags & NPY_C_CONTIGUOUS) ) 
	{
	  PyErr_SetString(PyExc_ValueError,
			  "Array must be C contiguous for this operation.");
	  return NULL;      
	}
      if( !(mapUin->flags & NPY_C_CONTIGUOUS) ) 
	{
	  PyErr_SetString(PyExc_ValueError,
			  "Array must be C contiguous for this operation.");
	  return NULL;      
	}
    }

  /* Check type of data : must be double ('d') */
  if( mapIin->descr->type != 'd' )
    {
      PyErr_SetString(PyExc_TypeError,
		      "Type must be Float64 for this function");
      return NULL;
    }

  if( polarisation )
    {
      if( mapQin->descr->type != 'd' )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Type must be Float64 for this function");
	  return NULL;
	}
      if( mapUin->descr->type != 'd' )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Type must be Float64 for this function");
	  return NULL;
	}
    }

  /* Check number of dimension : must be 1 */
  if( mapIin->nd != 1 )
    {
      PyErr_SetString(PyExc_TypeError,
		      "Array must be 1D.");
      return NULL;
    }
  else
    npix = mapIin->dimensions[0];

  if( polarisation )
    {
      if( mapQin->nd != 1 )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Array must be 1D.");
	  return NULL;
	}
      else
	if( mapQin->dimensions[0] != npix)
	  {
	    PyErr_SetString(PyExc_TypeError,
			    "All maps must have same dimension.");
	    return NULL;
	  }
      if( mapUin->nd != 1 )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Array must be 1D.");
	  return NULL;
	}
      else
	if( mapUin->dimensions[0] != npix)
	  {
	    PyErr_SetString(PyExc_TypeError,
			    "All maps must have same dimension.");
	    return NULL;
	  }
    }
  
  /* Check that the number of pixel is ok (12*nside^2) */
  nside = npix2nside(npix);
  if( nside < 0 )
    {
      PyErr_SetString(PyExc_ValueError,
		      "Number of pixel not valid for healpix map.");
      return NULL;
    }

  /* lmax and mmax */
  if( lmax < 0 )
    lmax = 3*nside-1;
  if( mmax <0 || mmax > lmax )
    mmax = lmax;

  Healpix_Map<double> mapI;
  {
    arr<double> arr_map((double*)mapIin->data, npix);
    mapI.Set(arr_map, RING);
  }

  Healpix_Map<double> mapQ;
  if( polarisation ) 
    {
      arr<double> arr_map((double*)mapQin->data, npix);
      mapQ.Set(arr_map, RING);
    }

  Healpix_Map<double> mapU;
  if( polarisation )
    {
      arr<double> arr_map((double*)mapUin->data, npix);
      mapU.Set(arr_map, RING);
    }

  npy_intp szalm = Alm<xcomplex<double> >::Num_Alms(lmax,mmax);

  PyArrayObject *almIout = NULL;
  almIout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szalm, 
					      PyArray_CDOUBLE);
  if( !almIout ) return NULL;

  PyArrayObject *almGout = NULL;
  if( polarisation )
    {
      almGout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szalm, 
						  PyArray_CDOUBLE);
      if( !almGout ) {
	Py_DECREF(almIout);
	return NULL;
      }
    }

  PyArrayObject *almCout = NULL;
  if( polarisation )
    {
      almCout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szalm, 
						  PyArray_CDOUBLE);
      if( !almCout ) {
	Py_DECREF(almIout);
	Py_DECREF(almGout);
	return NULL;
      }
    }

  Alm< xcomplex<double> > almIalm;
  {
    arr< xcomplex<double> > alm_arr((xcomplex<double>*)almIout->data, szalm);
    almIalm.Set(alm_arr, lmax, mmax);
  }

  Alm< xcomplex<double> > almGalm;
  if( polarisation )
    {
      arr< xcomplex<double> > alm_arr((xcomplex<double>*)almGout->data, szalm);
      almGalm.Set(alm_arr, lmax, mmax);
    }

  Alm< xcomplex<double> > almCalm;
  if( polarisation )
    {
      arr< xcomplex<double> > alm_arr((xcomplex<double>*)almCout->data, szalm);
      almCalm.Set(alm_arr, lmax, mmax);
    }

  arr<double> weight;

  if( use_weights )
    {
      read_weight_ring(datapath, nside, weight);
      for (int m=0; m<weight.size(); ++m) weight[m]+=1;
    }
  else
    {
      weight.alloc(2*nside);
      weight.fill(1.);
    }

  double avg = 0.0;
  if( regression ) {
    avg = mapI.average();
    mapI.add(-avg);
  }

  if( !polarisation )
    map2alm_iter(mapI,almIalm,num_iter,weight);
  else
    map2alm_pol_iter(mapI, mapQ, mapU, almIalm, almGalm, almCalm, num_iter,
		     weight);

  if( regression ) {
    almIalm(0,0) += avg*sqrt(fourpi);
    mapI.add(avg);
  }

  if( !docl )
    {
      if( !polarisation )
	return Py_BuildValue("N",almIout);
      else
	return Py_BuildValue("NNN", almIout, almGout, almCout);
    }
  else
    {
      if( !polarisation )
	{
	  PowSpec powspec;
	  extract_powspec (almIalm,powspec);
	  npy_intp szcl = (npy_intp)(powspec.Lmax()+1);
	  
	  PyArrayObject *ctt=NULL;
	  
	  ctt = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szcl, 
						  PyArray_DOUBLE);
	  if( !ctt ) 
	    return NULL;
	  
	  for( int l=0; l<szcl; l++ )
	    *((double*)PyArray_GETPTR1(ctt,l)) =  powspec.tt(l);
	  return Py_BuildValue("NN",ctt,almIout);
	}
      else
	{
	  PowSpec powspec;
	  extract_powspec(almIalm, almGalm,almCalm,powspec);
	  
	  npy_intp szcl = (npy_intp)(powspec.Lmax()+1);
	  
	  PyArrayObject *ctt=NULL;
	  PyArrayObject *cee=NULL;
	  PyArrayObject *cbb=NULL;
	  PyArrayObject *cte=NULL;

	  ctt = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szcl, 
						  PyArray_DOUBLE);
	  if( !ctt ) 
	    return NULL;
	  cee = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szcl, 
						  PyArray_DOUBLE);
	  if( !cee ) 
	    return NULL;
	  cbb = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szcl, 
						  PyArray_DOUBLE);
	  if( !cbb ) 
	    return NULL;
	  cte = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szcl, 
						  PyArray_DOUBLE);
	  if( !cte ) 
	    return NULL;

	  for( int l=0; l<szcl; l++ )
	    {
	      *((double*)PyArray_GETPTR1(ctt,l)) =  powspec.tt(l);
	      *((double*)PyArray_GETPTR1(cee,l)) =  powspec.gg(l);
	      *((double*)PyArray_GETPTR1(cbb,l)) =  powspec.cc(l);
	      *((double*)PyArray_GETPTR1(cte,l)) =  powspec.tg(l);
	    }
	  return Py_BuildValue("(NNNN)(NNN)",ctt,cee,cbb,cte,
			       almIout, almGout, almCout);
	}
    }
}

/***********************************************************************
    alm2map

       input: cl or alm, nside, iter=3, use_weights=True

       output: map in RING scheme
*/
static PyObject *healpy_alm2map(PyObject *self, PyObject *args, 
				PyObject *kwds)
{
  PyArrayObject *almIin = NULL;
  PyArrayObject *almGin = NULL;
  PyArrayObject *almCin = NULL;
  int nside = 64;
  int lmax = -1;
  int mmax = -1;
  int polarisation = 0;
  
  static char* kwlist[] = {"","nside", "lmax", "mmax", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iii", kwlist,
				   &PyArray_Type, &almIin,
				   &nside,
				   &lmax, 
				   &mmax))
				   
    {
      PyErr_Clear(); /* I want to try the other calling way */

      PyObject *t = NULL;
      if( !PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", kwlist,
				       &t,
				       &nside,
				       &lmax, 
				       &mmax) )
	return NULL;
      else
	{
	  if( PySequence_Size(t) != 3 )
	    {
	      PyErr_SetString(PyExc_TypeError, 
			      "First argument must be a sequence with "
			      "three elements.");
	      return NULL;
	    }
	  PyObject *o1 = NULL;
	  PyObject *o2 = NULL;
	  PyObject *o3 = NULL;
	  o1 = PySequence_GetItem(t, 0);
	  o2 = PySequence_GetItem(t, 1);
	  o3 = PySequence_GetItem(t, 2);
	  /* I decrease reference counts here,
	     because PySequence_GetItem increase 
	     reference count, and I just want to 
	     borrow a reference the time of this 
	     function. */
	  Py_XDECREF(o1);
	  Py_XDECREF(o2);
	  Py_XDECREF(o3);
	  if( ! ( PyArray_Check(o1) && PyArray_Check(o2) 
		  && PyArray_Check(o3) ) )
	    {
	      PyErr_SetString(PyExc_TypeError, 
			      "First argument must be a sequence with "
			      "three arrays");
	      return NULL;
	    }
	  else
	    {
	      almIin = (PyArrayObject*) o1;
	      almGin = (PyArrayObject*) o2;
	      almCin = (PyArrayObject*) o3;
	    }
	}
	polarisation = 1;  /* we have three maps : polarisation! */
    }
    
  /* Check array is contiguous */
  if( !(almIin->flags & NPY_C_CONTIGUOUS) ) 
    {
      PyErr_SetString(PyExc_ValueError,
		      "Array must be C contiguous for this operation.");
      return NULL;      
    }
  if( polarisation ) 
    {
      if( !(almGin->flags & NPY_C_CONTIGUOUS) ) 
	{
	  PyErr_SetString(PyExc_ValueError,
			  "Array must be C contiguous for this operation.");
	  return NULL;      
	}
      if( !(almCin->flags & NPY_C_CONTIGUOUS) ) 
	{
	  PyErr_SetString(PyExc_ValueError,
			  "Array must be C contiguous for this operation.");
	  return NULL;      
	}
    }
  
  /* Check type of data : must be double, real ('d') or complex ('D') */
  if( almIin->descr->type != 'D' )
    {
      PyErr_SetString(PyExc_TypeError,
		      "Type must be Complex for this function");
      return NULL;
    }
  if( polarisation )
    {
      if( almIin->descr->type != 'D' )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Type must be Complex for this function");
	  return NULL;
	}
      if( almIin->descr->type != 'D' )
	{
	  PyErr_SetString(PyExc_TypeError,
			  "Type must be Complex for this function");
	  return NULL;
	}
    }

  /* Check number of dimension : must be 1 */
  if( almIin->nd != 1 )
    {
      PyErr_SetString(PyExc_ValueError,
		      "The map must be a 1D array");
      return NULL;
    }
  if( polarisation )
    {
      if( almGin->nd != 1 )
	{
	  PyErr_SetString(PyExc_ValueError,
			  "The map must be a 1D array");
	  return NULL;
	}
      if( almCin->nd != 1 )
	{
	  PyErr_SetString(PyExc_ValueError,
			  "The map must be a 1D array");
	  return NULL;
	}
    }

  /* Need to have lmax and mmax defined */
  if( lmax < 0 )
    {
      /* Check that the dimension is compatible with lmax=mmax */
      long imax = almIin->dimensions[0] - 1;
      double ell;
      ell = (-3.+sqrt(9.+8.*imax))/2.;
      if( ell != floor(ell) )
	{
	  PyErr_SetString(PyExc_ValueError, "Wrong alm size "
			  "(or give lmax and mmax).\n");
	  return NULL;
	}
      lmax=(int)floor(ell);
      mmax = lmax;
    }
  if( mmax < 0 || mmax > lmax)
    mmax = lmax;

  /* Check lmax and mmax are ok compared to alm.size */
  int szalm = Alm< xcomplex<double> >::Num_Alms(lmax,mmax);
  if( almIin->dimensions[0] != szalm )
    {
      PyErr_SetString(PyExc_ValueError, "Wrong alm size.\n");
      return NULL;
    }
  if( polarisation )
    {
      if( almIin->dimensions[0] != szalm )
	{
	  PyErr_SetString(PyExc_ValueError, "Wrong alm size.\n");
	  return NULL;
	}
      if( almIin->dimensions[0] != szalm )
	{
	  PyErr_SetString(PyExc_ValueError, "Wrong alm size.\n");
	  return NULL;
	}
    }

  /* Now we can build an Alm and give it to alm2map_iter */
  Alm< xcomplex<double> > almIalm;
  {
    arr< xcomplex<double> > alm_arr((xcomplex<double>*)almIin->data, szalm);
    almIalm.Set(alm_arr, lmax, mmax);
  }
  Alm< xcomplex<double> > almGalm;
  if( polarisation ) 
    {
      arr< xcomplex<double> > alm_arr((xcomplex<double>*)almGin->data, szalm);
      almGalm.Set(alm_arr, lmax, mmax);
    }
  Alm< xcomplex<double> > almCalm;
  if( polarisation )
    {
      arr< xcomplex<double> > alm_arr((xcomplex<double>*)almCin->data, szalm);
      almCalm.Set(alm_arr, lmax, mmax);
    }
  
  /* We must prepare the map */
  
  npy_intp npix = nside2npix(nside);
  PyArrayObject *mapIout = NULL;
  mapIout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
					      PyArray_DOUBLE);
  if( !mapIout ) 
    return NULL;

  PyArrayObject *mapQout = NULL;
  if( polarisation )
    {
      mapQout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
						  PyArray_DOUBLE);
      if( !mapQout ) 
	return NULL;
    }

  PyArrayObject *mapUout = NULL;
  if( polarisation )
    {
      mapUout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
						  PyArray_DOUBLE);
      if( !mapUout ) 
	return NULL;
    }

  Healpix_Map<double> mapI;
  {
    arr<double> arr_map((double*)mapIout->data, npix);
    mapI.Set(arr_map, RING);
  }

  Healpix_Map<double> mapQ;
  if( polarisation ) 
    {
      arr<double> arr_map((double*)mapQout->data, npix);
      mapQ.Set(arr_map, RING);
    }
  Healpix_Map<double> mapU;
  if( polarisation )
    {
      arr<double> arr_map((double*)mapUout->data, npix);
      mapU.Set(arr_map, RING);
    }
  
  /* We now call alm2map */

  if( !polarisation )
    {
      double offset = almIalm(0,0).real()/sqrt(fourpi);
      xcomplex<double> almI00 = almIalm(0,0);
      almIalm(0,0) = 0;
      alm2map(almIalm,mapI);
      mapI.add(offset);
      almIalm(0,0) = almI00;
    }
  else
    {
      double offset = almIalm(0,0).real()/sqrt(fourpi);
      xcomplex<double> almI00 = almIalm(0,0);
      almIalm(0,0) = 0;
      alm2map_pol(almIalm,almGalm,almCalm,mapI,mapQ,mapU);
      mapI.add(offset);      
      almIalm(0,0) = almI00;
    }

  /* We now return the map */
  if( !polarisation )
    return Py_BuildValue("N",mapIout);
  else
    return Py_BuildValue("NNN",mapIout,mapQout,mapUout);
}

/***********************************************************************
    alm2map_der1

       input: alm, nside, 

       output: map in RING scheme
*/
static PyObject *healpy_alm2map_der1(PyObject *self, PyObject *args, 
        PyObject *kwds) {
  PyArrayObject *almIin = NULL;
  int nside = 64;
  int lmax = -1;
  int mmax = -1;
  
  static char* kwlist[] = {"","nside", "lmax", "mmax", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iii", kwlist,
           &PyArray_Type, &almIin,
           &nside,
           &lmax, 
           &mmax)) {
    return NULL;
  } 
    
  /* Check array is contiguous */
  if( !(almIin->flags & NPY_C_CONTIGUOUS) ) {
      PyErr_SetString(PyExc_ValueError,
          "Array must be C contiguous for this operation.");
      return NULL;      
  }
  
  /* Check type of data : must be double, real ('d') or complex ('D') */
  if( almIin->descr->type != 'D' )  {
      PyErr_SetString(PyExc_TypeError,
          "Type must be Complex for this function");
      return NULL;
  }

  /* Check number of dimension : must be 1 */
  if( almIin->nd != 1 ) {
      PyErr_SetString(PyExc_ValueError,
          "The map must be a 1D array");
      return NULL;
    }
  
  /* Need to have lmax and mmax defined */
  if( lmax < 0 ) {
      /* Check that the dimension is compatible with lmax=mmax */
      long imax = almIin->dimensions[0] - 1;
      double ell;
      ell = (-3.+sqrt(9.+8.*imax))/2.;
      if( ell != floor(ell) ) {
        PyErr_SetString(PyExc_ValueError, "Wrong alm size "
          "(or give lmax and mmax).\n");
        return NULL;
      }
      lmax=(int)floor(ell);
      mmax = lmax;
  }
  if( mmax < 0 || mmax > lmax) {
    mmax = lmax;
  }
  
  /* Check lmax and mmax are ok compared to alm.size */
  int szalm = Alm< xcomplex<double> >::Num_Alms(lmax,mmax);
  if( almIin->dimensions[0] != szalm ) {
      PyErr_SetString(PyExc_ValueError, "Wrong alm size.\n");
      return NULL;
  }
  
  /* Now we can build an Alm and give it to alm2map_iter */
  Alm< xcomplex<double> > almIalm;
  {
    arr< xcomplex<double> > alm_arr((xcomplex<double>*)almIin->data, szalm);
    almIalm.Set(alm_arr, lmax, mmax);
  }
  
  /* We must prepare the map */
  
  npy_intp npix = nside2npix(nside);

  PyArrayObject *mapIout = NULL;
  mapIout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
                PyArray_DOUBLE);
  if( !mapIout ) 
    return NULL;
  Healpix_Map<double> mapI;
  {
    arr<double> arr_map((double*)mapIout->data, npix);
    mapI.Set(arr_map, RING);
  }

  PyArrayObject *mapDtheta = NULL;
  mapDtheta = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
                PyArray_DOUBLE);
  if( !mapDtheta ) 
    return NULL;
  Healpix_Map<double> mapDt;
  {
    arr<double> arr_map((double*)mapDtheta->data, npix);
    mapDt.Set(arr_map, RING);
  }

  PyArrayObject *mapDphi = NULL;
  mapDphi = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&npix, 
                PyArray_DOUBLE);
  if( !mapDphi ) 
    return NULL;
  Healpix_Map<double> mapDp;
  {
    arr<double> arr_map((double*)mapDphi->data, npix);
    mapDp.Set(arr_map, RING);
  }


  /* We now call alm2map_der1 */

  double offset = almIalm(0,0).real()/sqrt(fourpi);
  xcomplex<double> almI00 = almIalm(0,0);
  almIalm(0,0) = 0;
  alm2map_der1(almIalm,mapI,mapDt,mapDp);
  mapI.add(offset);
  almIalm(0,0) = almI00;

  return Py_BuildValue("NNN",mapIout,mapDtheta,mapDphi);
}

long nside2npix(long nside)
{
  return 12*nside*nside;
}

long npix2nside(long npix)
{
  long nside;
  long npix2;
  nside = (long)floor(sqrt(npix/12.));
  npix2 = 12*nside*nside;
  if( npix2 == npix )
    return nside;
  else
    return -1;
}

/*  Functions needed to create alm from cl
    in the polarised case.
    The idea is as follow:
      - this function take as argument a sequence of n(n+1)/2 arrays,
        where n is the number of components 
      - these arrays represents the symetric correlation matrices
        of n components
      - this functions use the cholesky decomposition to compute
        the 
        
 */
/***********************************************************************
    synalm

       input: a sequence of Cl(i,j) where 0<=i,j<n denoting the cross-correlation
              only i>=j are given.

       output: alm
*/
static PyObject *healpy_synalm(PyObject *self, PyObject *args, 
			       PyObject *kwds)
{
  int lmax=-1, mmax=-1;
  int ncl, nalm;
  const double sqrt_two = sqrt(2.);

  /* Take also a sequence of unit variance random vectors for the alm. */

  static char* kwlist[] = {"","", "", "", NULL};

  PyObject *t = NULL;
  PyObject *u = NULL;
  DBGPRINTF("Parsing keyword\n");
  if( !PyArg_ParseTupleAndKeywords(args, kwds, "OOii", kwlist,
				   &t, 
				   &u,
				   &lmax, &mmax) )
    return NULL;

  DBGPRINTF("Checking sequence\n");
  if( (!PySequence_Check(t)) || (!PySequence_Check(u)) )
    return NULL;

  ncl = PySequence_Size(t);
  nalm = getn(ncl);
  DBGPRINTF("Sequence size: ncl=%d, nalm=%d\n", ncl, nalm);
  if( nalm <= 0 )
    {
      PyErr_SetString(PyExc_TypeError, 
		      "First argument must be a sequence with "
		      "n(n+1)/2 elements, and second argument "
		      "a sequence with n elements.");
      return NULL;
    }
  if( PySequence_Size(u) != nalm )
    {
      PyErr_SetString(PyExc_TypeError, 
		      "First argument must be a sequence with "
		      "n(n+1)/2 elements, and second argument "
		      "a sequence with n elements.");
      return NULL;
    }

  DBGPRINTF("Allocating memory\n");
  /* Memory allocation */
  PyArrayObject **cls = NULL;
  PyArrayObject **alms = NULL;
  Alm< xcomplex<double> > *almalms = NULL;
  double *mat = NULL;
  double *res = NULL;
  XMALLOC(cls, PyArrayObject*, ncl);
  XMALLOC(alms, PyArrayObject*, nalm);
  XNEW(almalms, Alm< xcomplex<double> >, nalm);
  XMALLOC(mat, double, ncl);
  XMALLOC(res, double, ncl);

  /*  From now on, I should do a 'goto fail' to return from the function
      in order to free allocated memory */
  /* Get the cls objects.
     If an object is None, set the array to NULL 
  */
  for( int i=0; i<ncl; i++ )
    {
      DBGPRINTF("Get item cl %d/%d\n", i+1, ncl);
      PyObject *o;
      o = PySequence_GetItem(t,i);
      /* I decrease reference counts here,
	 because PySequence_GetItem increase 
	 reference count, and I just want to 
	 borrow a reference the time of this 
	 function. */	  
      Py_XDECREF(o);
      if( o == Py_None )
	{
	  cls[i] = NULL;
	  DBGPRINTF("Cls[%d] is None\n", i);
	}
      else if( ! PyArray_Check(o) )
	{
	  PyErr_SetString(PyExc_TypeError, 
			  "First argument must be a sequence of "
			  "arrays");
	  goto fail;
	}
      else
	cls[i] = (PyArrayObject*) o;
    }
  for( int i=0; i<nalm; i++ )
    {
      PyObject *o;
      DBGPRINTF("Get item alm %d/%d\n", i+1, nalm);
      o = PySequence_GetItem(u,i);
      /* I decrease reference counts here,
	 because PySequence_GetItem increase 
	 reference count, and I just want to 
	 borrow a reference the time of this 
	 function. */	  
      Py_XDECREF(o);
      if( ! PyArray_Check(o) )
	{
	  PyErr_SetString(PyExc_TypeError, 
			  "First argument must be a sequence of "
			  "arrays");
	  goto fail;
	}
      alms[i] = (PyArrayObject*) o;
    }
  if( lmax<0 )
    {
      PyErr_SetString(PyExc_ValueError, 
		      "lmax must be positive.");
      goto fail;
    }
  if( mmax <0 || mmax >lmax )
    mmax=lmax;

  DBGPRINTF("lmax=%d, mmax=%d\n", lmax, mmax);

  /* Now, I check the arrays cls and alms are 1D and complex for alms */ 
  DBGPRINTF("Check dimension and size of cls\n");
  for( int i=0; i<ncl; i++ )
    {
      if( cls[i] == NULL ) 
	continue;
      if( (cls[i]->nd != 1) 
	  || ((cls[i]->descr->type != 'd') && (cls[i]->descr->type != 'f')) )
	{
	  PyErr_SetString(PyExc_TypeError,
		      "Type of cls must be float64 and arrays must be 1D.");
	  goto fail;
	}
    }
  DBGPRINTF("Check dimension and size of alms\n");
  for( int i=0; i<nalm; i++ )
    {
      if( (alms[i]->nd != 1) || (alms[i]->descr->type != 'D') )
	{
	  PyErr_SetString(PyExc_TypeError,
		      "Type of alms must be complex128 and arrays must be 1D.");
	  goto fail;
	}
    }
  
  /* Now, I check that all alms have the same size and that it is compatible with
     lmax and mmax */
  DBGPRINTF("Check alms have identical size\n");
  int szalm;
  szalm = -1;
  for( int i=0; i<nalm; i++ )
    {
      if( i==0 )
	szalm = alms[i]->dimensions[0];
      else if( alms[i]->dimensions[0] != szalm )
	{
	  PyErr_SetString(PyExc_ValueError,
			  "All alms arrays must have same size.");
	  goto fail;
	}      
    }
  if( szalm != Alm< xcomplex<double> >::Num_Alms(lmax,mmax) )
    {
      PyErr_SetString(PyExc_ValueError,
		      "lmax and mmax are not compatible with size of alms.");
      goto fail;
    }
  DBGPRINTF("Alms have all size %d\n", szalm);

  /* Set the objects Alm */
  DBGPRINTF("Set alm objects\n");
  for( int i=0; i<nalm; i++)
    {
      DBGPRINTF("Setting almalms[%d]\n", i);
      arr< xcomplex<double> > * alm_arr;
      alm_arr = new arr< xcomplex<double> >((xcomplex<double>*)alms[i]->data, szalm);
      DBGPRINTF("Set...\n");
      almalms[i].Set(*alm_arr, lmax, mmax);
      delete alm_arr;
    }


  /* Now, I can loop over l,
     for each l, I make the Cholesky decomposition of the correlation matrix
     given by the cls[*][l] 
  */
  DBGPRINTF("Start loop over l\n");
  for( int l=0; l<=lmax; l++ )
    {      
      DBGPRINTF("l=%d\n", l);
      /* fill the matrix of cls */
      for( int i=0; i<ncl; i++ )
	{
	  if( cls[i] == NULL )
	    mat[i] = 0.0;
	  else if( cls[i]->dimensions[0] < l )
	    mat[i] = 0.0;
	  else
	    {
	      if( cls[i]->descr->type == 'f' )
		mat[i] = (double)(*((float*)PyArray_GETPTR1(cls[i],l)));
	      else
		mat[i] = *((double*)PyArray_GETPTR1(cls[i],l));
	    }
	}

      /* Make the Cholesky decomposition */
      cholesky(nalm, mat, res);

      if( l == 400 )
	{
	  DBGPRINTF("matrice: ");
	  for( int i=0; i<ncl; i++ )
	    DBGPRINTF("%d: %lg  ", i, mat[i]);
	  DBGPRINTF("\n");
	  DBGPRINTF("cholesky: ");
	  for( int i=0; i<ncl; i++ )
	    DBGPRINTF("%d: %lg  ", i, res[i]);
	  DBGPRINTF("\n");
	}

      /* Apply the matrix to each m */

      /* m=0 */
      DBGPRINTF("   m=%d: ", 0);
      for( int i=nalm-1; i>=0; i-- )
	{
	  double x;
	  x = 0.0;
	  almalms[i](l,0).im = 0.0;
	  for( int j=0; j<=i; j++ )
	    x += res[getidx(nalm,i,j)]*almalms[j](l,0).re;
	  almalms[i](l,0).re = x;	  
	  DBGPRINTF(" %lg %lg ;", almalms[i](l,0).re, almalms[i](l,0).im);
	}
      DBGPRINTF("\n");

      /* m > 1 */
      for( int m=1; m<=l; m++ )
	{
	  DBGPRINTF("   m=%d: ", m);
	  for( int i=nalm-1; i>=0; i-- )
	    {
	      double xr, xi;
	      xi = xr = 0.0;
	      for( int j=0; j<=i; j++ )
		{
		  xr += res[getidx(nalm,i,j)]*almalms[j](l,m).re;
		  xi += res[getidx(nalm,i,j)]*almalms[j](l,m).im;
		  DBGPRINTF("(res[%d]=%lg, alm=%lg,%lg) %lg %lg", (int)getidx(nalm,i,j), 
			    res[getidx(nalm,i,j)],
			    almalms[j](l,m).re, almalms[j](l,m).im, 
			    xr, xi);
		}
	      almalms[i](l,m).re = xr/sqrt_two;
	      almalms[i](l,m).im = xi/sqrt_two;
	      DBGPRINTF(" xre,xim[%d]: %lg %lg ;", i, 
			almalms[i](l,m).re, almalms[i](l,m).im);
	    }
	  DBGPRINTF("\n");
      }
   }

  /* Should be finished now... */
  XFREE(cls);
  XFREE(alms);
  XDELETE(almalms);
  XFREE(mat);
  XFREE(res);
  Py_INCREF(Py_None);
  return Py_None;
 
  /* To be done before returning */
 fail:
  XFREE(cls);
  XFREE(alms);
  XDELETE(almalms);
  XFREE(mat);
  XFREE(res);

  return NULL;
}

PyObject *healpy_getn(PyObject *self, PyObject *args)
{
  long s;
  if( ! PyArg_ParseTuple(args, "l", &s) )
    {
      PyErr_SetString(PyExc_TypeError, "This function takes an integer as argument.");
      return NULL;
    }
  long n;
  n = getn(s);
  return Py_BuildValue("l",n);
}

/* Helpers functions */

/*
void getij(long n, long idx, long *i, long *j)
{
  *i = (long)ceil(((2*n+1)-sqrt((2*n+1)*(2*n+1)-8*(idx-n)))/2);
  *j = idx - (*i)*(2*n+1- (*i) )/2;
}
*/

long getidx(long n, long i, long j)
{
  long tmp;
  if( i > j )
    {tmp = j; j=i; i=tmp;}
  return i*(2*n-1-i)/2+j;
}


long getn(long s)
{
  long x;
  x = (long)floor((-1+sqrt(1+8*s))/2);
  if( (x*(x+1)/2) != s )
    return -1;
  else
    return x;
}

void cholesky(int n, double *data, double *res)
{
  int i,j,k;
  double sum;

  for( j=0; j<n; j++ )
    {
      for( i=0; i<n; i++ )
	{
	  if( i==j )
	    {
	      sum = data[getidx(n,j,j)];
	      for(k=0; k<j; k++ )
		sum -= res[getidx(n,k,j)]*res[getidx(n,k,j)];
	      if( sum <= 0 )
		res[getidx(n,j,j)] = 0.0;
	      else
		res[getidx(n,j,j)] = sqrt(sum);
	    }
	  else if( i>j)
	    {
	      sum = data[getidx(n,i,j)];
	      for( k=0; k<j; k++ )
		sum -= res[getidx(n,i,k)]*res[getidx(n,j,k)];
	      if( res[getidx(n,j,j)] != 0.0 )
		res[getidx(n,i,j)] = sum/res[getidx(n,j,j)];
	      else
		res[getidx(n,i,j)] = 0.0;
	    }
	}
    }
  return;
}


/***********************************************************************
    alm2signal

       input: 
          alm: a 1D ndarray
          theta:
          phi:    direction in the sky

       parameter: lmax, mmax (default, lmax=mmax compatible with szalm)

       output: the signal at given position
*/
static PyObject *healpy_alm2signal(PyObject *self, PyObject *args, 
				   PyObject *kwds)
{
  int lmax=-1, mmax=-1;
  double theta,phi;

  static char* kwlist[] = {"","", "", "lmax", "mmax", NULL};

  PyArrayObject *alm = NULL;

  if( !PyArg_ParseTupleAndKeywords(args, kwds, "O!dd|ii", kwlist,
				   &PyArray_Type, &alm,
				   &theta, &phi,
				   &lmax, &mmax) )
    return NULL;
  
  /* Check array is contiguous */
  if( !(alm->flags & NPY_C_CONTIGUOUS) ) 
    {
      PyErr_SetString(PyExc_ValueError,
		      "Array must be C contiguous for this operation.");
      return NULL;      
    }
  
  /* Check type of data : must be double complex ('D') */
  if( alm->descr->type != 'D' )
    {
      PyErr_SetString(PyExc_TypeError,
		      "Type must be double complex (complex128) "
		      "for this function");
      return NULL;
    }

  /* Check number of dimension : must be 1 */
  if( alm->nd != 1 )
    {
      PyErr_SetString(PyExc_ValueError,
		      "The map must be a 1D array");
      return NULL;
    }
  /* Need to have lmax and mmax defined */
  if( lmax < 0 )
    {
      /* Check that the dimension is compatible with lmax=mmax */
      long imax = alm->dimensions[0] - 1;
      double ell;
      ell = (-3.+sqrt(9.+8.*imax))/2.;
      if( ell != floor(ell) )
	{
	  PyErr_SetString(PyExc_ValueError, "Wrong alm size "
			  "(or give lmax and mmax).\n");
	  return NULL;
	}
      lmax=(int)floor(ell);
      mmax = lmax;
    }
  if( mmax < 0 || mmax > lmax)
    mmax = lmax;

  /* Check lmax and mmax are ok compared to alm.size */
  int szalm = Alm< xcomplex<double> >::Num_Alms(lmax,mmax);
  if( alm->dimensions[0] != szalm )
    {
      PyErr_SetString(PyExc_ValueError, "Wrong alm size.\n");
      return NULL;
    }

  /* Initialize the Alm object */
  Alm< xcomplex<double> > almobj;
  {
    arr< xcomplex<double> > alm_arr((xcomplex<double>*)alm->data, szalm);
    almobj.Set(alm_arr, lmax, mmax);
  }

  /* call sky_signal_direct */
  double result;

  result = sky_signal_direct(almobj, theta, phi);

  /* return the value */
  return Py_BuildValue("d",result);

}



double sky_signal_direct(Alm<xcomplex<double> >alm, double theta, double phi)
{
  // 
  int lmax = alm.Lmax();
  int mmax = alm.Mmax();
  Ylmgen ylmgen(lmax, mmax);

  double cth, sth;
  cth = cos(theta);
  sth = sin(theta);

  arr<double> ylm(lmax+1);

  // to do as in alm2map_cxx
  double offset = alm(0,0).real()/sqrt(fourpi);
  alm(0,0) = 0;

  // Term m=0
  int firstl = 0;
  ylmgen.get_Ylm(cth, sth, 0, ylm, firstl);
  double a = 0.0;
  for( int l=firstl; l<=lmax; l++ )
    {
      a += alm(l,0).real() * ylm[l];
    }
  
  // Terms m>=1
  double b = 0.0;
  for( int m=1; m<=mmax; m++ )
    {
      ylmgen.get_Ylm(cth, sth, m, ylm, firstl);      
      for( int l=firstl; l<= lmax; l++ )
	{
	  b += ( alm(l,m).real() * cos(m*phi) - 
		 alm(l,m).imag() * sin(m*phi) )*ylm[l];
	}
    }
  
  return (a + 2. * b) + offset;
}



static PyObject *healpy_Ylm(PyObject *self, PyObject *args, 
			    PyObject *kwds)
{
  int lmax, mmax, m;
  double theta, cth, sth;

  static char* kwlist[] = {"", "", "", "", NULL};


  if( !PyArg_ParseTupleAndKeywords(args, kwds, "iiid", kwlist,
				   &lmax, &mmax, &m, &theta) )
    return NULL;

  /* Compute the Ylm */
  
  Ylmgen ylmgen(lmax, mmax);
  arr<double> ylm;
  int firstl;

  cth = cos(theta);
  sth = sin(theta);
  ylmgen.get_Ylm(cth, sth, m, ylm, firstl);

  /* Create a numpy array */
  PyArrayObject *ylm_arr;
  int nd = 1;
  npy_intp dims[1];
  dims[0] = ylm.size();
  
  ylm_arr = (PyArrayObject*)PyArray_SimpleNew(nd, dims, NPY_DOUBLE);
  if( !ylm_arr )
    return NULL;

  /* copy values of ylm to ylm_arr */
  memcpy(PyArray_BYTES(ylm_arr), ylm.begin(), PyArray_DIM(ylm_arr, 0)*sizeof(double));

  /* Set lower values to zero */
  for( int i=0; i<firstl; i++ )
    *(double*)PyArray_GETPTR1(ylm_arr, i) = 0.0;

  /* return the array object */
  return Py_BuildValue("N",ylm_arr);
}

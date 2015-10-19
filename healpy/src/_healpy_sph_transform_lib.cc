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
#include <iostream>

#include "arr.h"
#include "alm.h"
#include "lsconstants.h"
#include "healpix_map.h"
#include "xcomplex.h"
#include "alm_healpix_tools.h"
#include "powspec.h"
#include "alm_powspec_tools.h"
#include "healpix_data_io.h"
#include "_healpy_utils.h"

#define IS_DEBUG_ON 0

/* Some helpful macro */
#define XMALLOC(X,Y,Z) if( !(X = (Y*)malloc(Z*sizeof(Y))) ) { PyErr_NoMemory(); goto fail;}
#define XNEW(X,Y,Z) if( !(X = new Y[Z]) ) { PyErr_NoMemory(); goto fail; }
#define XFREE(X) if( X ) free(X);
#define XDELETE(X) if( X ) delete[] X;
#define DBGPRINTF(X,...) if( IS_DEBUG_ON ) printf(X, ## __VA_ARGS__)


static long nside2npix(long nside)
{
  return 12*nside*nside;
}

static long npix2nside(long npix)
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


/* Helper functions */

/*
static void getij(long n, long idx, long *i, long *j)
{
  *i = (long)ceil(((2*n+1)-sqrt((2*n+1)*(2*n+1)-8*(idx-n)))/2);
  *j = idx - (*i)*(2*n+1- (*i) )/2;
}
*/

static long getidx(long n, long i, long j)
{
  long tmp;
  if( i > j )
    {tmp = j; j=i; i=tmp;}
  return i*(2*n-1-i)/2+j;
}


static long getn(long s)
{
  long x;
  x = (long)floor((-1+sqrt(1+8*s))/2);
  if( (x*(x+1)/2) != s )
    return -1;
  else
    return x;
}

static void cholesky(int n, double *data, double *res)
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
    map2alm

       input: map, lmax=3*nside-1, mmax=lmax, cl=False
              iter=3, use_weights=False

       output: alm (or cl if cl=True)
*/
static PyObject *healpy_map2alm(PyObject *self, PyObject *args,
                                PyObject *kwds)
{
  PyArrayObject *mapIin = NULL, *mapQin = NULL, *mapUin = NULL;
  int lmax=-1, mmax=-1;
  int nside=-1;
  int npix=-1;
  int num_iter=3;
  int docl=0;
  int use_weights=0;
  char * datapath=NULL;
  int polarisation = 0; /* not polarised by default */

  static const char* kwlist[] = {"","lmax", "mmax","cl","iter",
                           "use_weights", "data_path", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iiiiisi", (char **)kwlist,
                                   &PyArray_Type, &mapIin,
                                   &lmax, &mmax, &docl,
                                   &num_iter,&use_weights,&datapath))
    {
      PyErr_Clear(); /* I want to try the other calling way */

      PyObject *t = NULL;
      if( !PyArg_ParseTupleAndKeywords(args, kwds, "O|iiiiisi", (char **)kwlist,
                                   &t,
                                   &lmax, &mmax, &docl,
                                       &num_iter,&use_weights,&datapath))
        return NULL;
      else
        {
          healpyAssertType(PySequence_Size(t)==3,
            "First argument must be a sequence with three elements.");
          PyObject *o1 = PySequence_GetItem(t, 0),
                   *o2 = PySequence_GetItem(t, 1),
                   *o3 = PySequence_GetItem(t, 2);
          /* I decrease reference counts here,
             because PySequence_GetItem increase
             reference count, and I just want to
             borrow a reference the time of this
             function. */
          Py_XDECREF(o1);
          Py_XDECREF(o2);
          Py_XDECREF(o3);
          healpyAssertType(PyArray_Check(o1)&&PyArray_Check(o2)&&PyArray_Check(o3),
            "First argument must be a sequence with three arrays");

          mapIin = (PyArrayObject*) o1;
          mapQin = (PyArrayObject*) o2;
          mapUin = (PyArrayObject*) o3;
        }
        polarisation = 1;  /* we have three maps : polarisation! */
    }

  /* Check array is contiguous */
  healpyAssertValue(mapIin->flags&NPY_C_CONTIGUOUS,
    "Array must be C contiguous for this operation.");

  if( polarisation )
    healpyAssertValue(mapQin->flags&mapUin->flags&NPY_C_CONTIGUOUS,
                      "Array must be C contiguous for this operation.");

  /* Check type of data : must be double ('d') */
  healpyAssertType(mapIin->descr->type == 'd',
    "Type must be Float64 for this function");

  if( polarisation )
    {
    healpyAssertType(mapQin->descr->type == 'd',
      "Type must be Float64 for this function");
    healpyAssertType(mapUin->descr->type == 'd',
      "Type must be Float64 for this function");
    }

  /* Check number of dimension : must be 1 */
  healpyAssertType(mapIin->nd==1,"Array must be 1D.");
  npix = mapIin->dimensions[0];

  if( polarisation )
    {
    healpyAssertType((mapQin->nd==1)&&(mapUin->nd==1),"Array must be 1D.");
    healpyAssertType((mapQin->dimensions[0]==npix)&&
                     (mapQin->dimensions[0]==npix),
                      "All maps must have same dimension.");
    }

  /* Check that the number of pixel is ok (12*nside^2) */
  nside = npix2nside(npix);
  healpyAssertValue(nside>=0,"Number of pixel not valid for healpix map.");

  /* lmax and mmax */
  if( lmax < 0 )
    lmax = 3*nside-1;
  if( mmax <0 || mmax > lmax )
    mmax = lmax;

  Healpix_Map<double> mapI, mapQ, mapU;
  {
  arr<double> arr_map((double*)mapIin->data, npix);
  mapI.Set(arr_map, RING);
  }

  if(polarisation)
    {
    arr<double> arr_map((double*)mapQin->data, npix);
    mapQ.Set(arr_map, RING);
    }

  if( polarisation )
    {
    arr<double> arr_map((double*)mapUin->data, npix);
    mapU.Set(arr_map, RING);
    }

  npy_intp szalm = Alm<xcomplex<double> >::Num_Alms(lmax,mmax);

  PyArrayObject *almIout = (PyArrayObject*)PyArray_SimpleNew
    (1, (npy_intp*)&szalm, PyArray_CDOUBLE);
  if( !almIout ) return NULL;

  PyArrayObject *almGout=NULL, *almCout=NULL;
  if( polarisation )
    {
      almGout = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szalm,
                                                  PyArray_CDOUBLE);
      if( !almGout ) {
        Py_DECREF(almIout);
        return NULL;
      }
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
      for (tsize m=0; m<weight.size(); ++m) weight[m]+=1;
    }
  else
      weight.allocAndFill(2*nside,1.);

  if( !polarisation )
    map2alm_iter(mapI,almIalm,num_iter,weight);
  else
    map2alm_pol_iter(mapI, mapQ, mapU, almIalm, almGalm, almCalm, num_iter,
                     weight);

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

          PyArrayObject *ctt = (PyArrayObject*)PyArray_SimpleNew
            (1, (npy_intp*)&szcl, PyArray_DOUBLE);
          if( !ctt ) return NULL;
          PyArrayObject *cee = (PyArrayObject*)PyArray_SimpleNew
            (1, (npy_intp*)&szcl, PyArray_DOUBLE);
          if( !cee ) return NULL;
          PyArrayObject *cbb = (PyArrayObject*)PyArray_SimpleNew
            (1, (npy_intp*)&szcl, PyArray_DOUBLE);
          if( !cbb ) return NULL;
          PyArrayObject *cte = (PyArrayObject*)PyArray_SimpleNew
            (1, (npy_intp*)&szcl, PyArray_DOUBLE);
          if( !cte ) return NULL;

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
  PyArrayObject *almIin=NULL, *almGin=NULL, *almCin=NULL;
  int nside = 64;
  int lmax = -1;
  int mmax = -1;
  int polarisation = 0;

  static const char* kwlist[] = {"","nside", "lmax", "mmax", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iii", (char **)kwlist,
                                   &PyArray_Type, &almIin,
                                   &nside,
                                   &lmax,
                                   &mmax))

    {
      PyErr_Clear(); /* I want to try the other calling way */

      PyObject *t = NULL;
      if( !PyArg_ParseTupleAndKeywords(args, kwds, "O|iii", (char **)kwlist,
                                       &t,
                                       &nside,
                                       &lmax,
                                       &mmax) )
        return NULL;
      else
        {
          healpyAssertType(PySequence_Size(t)==3,
            "First argument must be a sequence with three elements.");
          PyObject *o1 = PySequence_GetItem(t, 0),
                   *o2 = PySequence_GetItem(t, 1),
                   *o3 = PySequence_GetItem(t, 2);
          /* I decrease reference counts here,
             because PySequence_GetItem increase
             reference count, and I just want to
             borrow a reference the time of this
             function. */
          Py_XDECREF(o1);
          Py_XDECREF(o2);
          Py_XDECREF(o3);
          healpyAssertType(PyArray_Check(o1)&&PyArray_Check(o2)&&PyArray_Check(o3),
            "First argument must be a sequence with three arrays");
          almIin = (PyArrayObject*) o1;
          almGin = (PyArrayObject*) o2;
          almCin = (PyArrayObject*) o3;
        }
        polarisation = 1;  /* we have three maps : polarisation! */
    }

  /* Check array is contiguous */
  healpyAssertValue(almIin->flags&NPY_C_CONTIGUOUS,
                      "Array must be C contiguous for this operation.");
  if( polarisation )
    healpyAssertValue(almGin->flags&almCin->flags&NPY_C_CONTIGUOUS,
                          "Array must be C contiguous for this operation.");

  /* Check type of data : must be double, real ('d') or complex ('D') */
  healpyAssertType(almIin->descr->type == 'D',
                      "Type must be Complex for this function");
  if( polarisation )
    {
      healpyAssertType(almGin->descr->type == 'D',
                          "Type must be Complex for this function");
      healpyAssertType(almCin->descr->type == 'D',
                          "Type must be Complex for this function");
    }

  /* Check number of dimension : must be 1 */
  healpyAssertType(almIin->nd==1,"The a_lm must be a 1D array.");
  if( polarisation )
    {
    healpyAssertType(almGin->nd==1,"The a_lm must be a 1D array.");
    healpyAssertType(almCin->nd==1,"The a_lm must be a 1D array.");
    }

  /* Need to have lmax and mmax defined */
  if( lmax < 0 )
    {
      /* Check that the dimension is compatible with lmax=mmax */
      long imax = almIin->dimensions[0] - 1;
      double ell = (-3.+sqrt(9.+8.*imax))/2.;
      healpyAssertType(ell==floor(ell),
        "Wrong alm size (or give lmax and mmax)");
      lmax=(int)floor(ell);
      mmax = lmax;
    }
  if( mmax < 0 || mmax > lmax)
    mmax = lmax;

  /* Check lmax and mmax are ok compared to alm.size */
  int szalm = Alm< xcomplex<double> >::Num_Alms(lmax,mmax);
  healpyAssertValue(almIin->dimensions[0]==szalm,"Wrong alm size.");
  if( polarisation )
    {
    healpyAssertValue(almGin->dimensions[0]==szalm,"Wrong alm size.");
    healpyAssertValue(almCin->dimensions[0]==szalm,"Wrong alm size.");
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
      mapI.Add(offset);
      almIalm(0,0) = almI00;
    }
  else
    {
      double offset = almIalm(0,0).real()/sqrt(fourpi);
      xcomplex<double> almI00 = almIalm(0,0);
      almIalm(0,0) = 0;
      alm2map_pol(almIalm,almGalm,almCalm,mapI,mapQ,mapU);
      mapI.Add(offset);
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

  static const char* kwlist[] = {"","nside", "lmax", "mmax", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iii", (char **)kwlist,
           &PyArray_Type, &almIin,
           &nside,
           &lmax,
           &mmax)) {
    return NULL;
  }

  /* Check array is contiguous */
  healpyAssertValue(almIin->flags & NPY_C_CONTIGUOUS,
    "Array must be C contiguous for this operation.");

  /* Check type of data : must be double, real ('d') or complex ('D') */
  healpyAssertType(almIin->descr->type == 'D',
    "Type must be Complex for this function");

  /* Check number of dimension : must be 1 */
  healpyAssertValue(almIin->nd==1, "The map must be a 1D array");

  /* Need to have lmax and mmax defined */
  if( lmax < 0 ) {
      /* Check that the dimension is compatible with lmax=mmax */
      long imax = almIin->dimensions[0] - 1;
      double ell = (-3.+sqrt(9.+8.*imax))/2.;
      healpyAssertValue(ell == floor(ell),
        "Wrong alm size (or give lmax and mmax).");
      lmax=(int)floor(ell);
      mmax = lmax;
  }
  if( mmax < 0 || mmax > lmax) {
    mmax = lmax;
  }

  /* Check lmax and mmax are ok compared to alm.size */
  int szalm = Alm< xcomplex<double> >::Num_Alms(lmax,mmax);
  healpyAssertValue(almIin->dimensions[0] == szalm,"Wrong alm size.");

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
  mapI.Add(offset);
  almIalm(0,0) = almI00;

  return Py_BuildValue("NNN",mapIout,mapDtheta,mapDphi);
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

  static const char* kwlist[] = {"","", "", "", NULL};

  PyObject *t = NULL;
  PyObject *u = NULL;
  DBGPRINTF("Parsing keyword\n");
  if( !PyArg_ParseTupleAndKeywords(args, kwds, "OOii", (char **)kwlist,
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
  healpyAssertType((nalm>0)&&(PySequence_Size(u)==nalm),
                      "First argument must be a sequence with "
                      "n(n+1)/2 elements, and second argument "
                      "a sequence with n elements.");

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
          //|| ((cls[i]->descr->type != 'd') && (cls[i]->descr->type != 'f')) )
          || (cls[i]->descr->type != 'd') )
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
  if( szalm != int(Alm< xcomplex<double> >::Num_Alms(lmax,mmax)) )
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
          almalms[i](l,0)=xcomplex<double>(almalms[i](l,0).real(),0.0);
          for( int j=0; j<=i; j++ )
            x += res[getidx(nalm,i,j)]*almalms[j](l,0).real();
          almalms[i](l,0)=xcomplex<double>(x,0.);
          DBGPRINTF(" %lg %lg ;", almalms[i](l,0).real(), almalms[i](l,0).imag());
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
                  xr += res[getidx(nalm,i,j)]*almalms[j](l,m).real();
                  xi += res[getidx(nalm,i,j)]*almalms[j](l,m).imag();
                  DBGPRINTF("(res[%d]=%lg, alm=%lg,%lg) %lg %lg", (int)getidx(nalm,i,j),
                            res[getidx(nalm,i,j)],
                            almalms[j](l,m).real(), almalms[j](l,m).imag(),
                            xr, xi);
                }
              almalms[i](l,m)=xcomplex<double>(xr/sqrt_two,xi/sqrt_two);
              DBGPRINTF(" xre,xim[%d]: %lg %lg ;", i,
                        almalms[i](l,m).real(), almalms[i](l,m).imag());
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
  healpyAssertType (PyArg_ParseTuple(args, "l", &s),
    "This function takes an integer as argument.");
  long n = getn(s);
  return Py_BuildValue("l",n);
}

static PyMethodDef methods[] = {
  {"_map2alm", (PyCFunction)healpy_map2alm, METH_VARARGS | METH_KEYWORDS,
   "Compute alm or cl from an input map.\n"
   "The input map is assumed to be ordered in RING.\n"
   "anafast(map,lmax=3*nside-1,mmax=lmax,cl=False,\n"
   "        iter=3,use_weights=False,data_path=None"},
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
  {NULL, NULL, 0, NULL} /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_healpy_sph_transform_lib",
  NULL, -1, methods
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION < 3
init_healpy_sph_transform_lib(void)
#else
PyInit__healpy_sph_transform_lib(void)
#endif
{
  import_array();

#if PY_MAJOR_VERSION < 3
  Py_InitModule("_healpy_sph_transform_lib", methods);
#else
  return PyModule_Create(&moduledef);
#endif
}

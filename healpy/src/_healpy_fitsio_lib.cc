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

#include <string>
#include <iostream>

#include "numpy/arrayobject.h"

#include "healpix_data_io.h"
#include "arr.h"

static PyObject *healpy_pixwin(PyObject *self, PyObject *args, 
			       PyObject *kwds);

static PyMethodDef HEALFITSMethods[] = {
  {"_pixwin", (PyCFunction)healpy_pixwin, METH_VARARGS | METH_KEYWORDS, 
   "Return the pixel window for some nside\n"
   "_pixwin(nside,data_path,pol=False)"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_healpy_fitsio_lib(void)
{
  PyObject *m;
  m =  Py_InitModule("_healpy_fitsio_lib", HEALFITSMethods);
  
  import_array();
}

/***********************************************************************
    healpy_pixwin

       input: nside, data_path, pol=False

       output: W(l), the pixel window function
*/
static PyObject *healpy_pixwin(PyObject *self, PyObject *args, 
			       PyObject *kwds)
{
  int nside;
  char * datapath=NULL;
  int polarisation = 0; /* not polarised by default */

  static char* kwlist[] = {"","", "pol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "is|i", kwlist,
				   &nside, &datapath,
				   &polarisation))
    return NULL;

  double x = log(nside)/log(2.);
  if( x != floor(x) )
    {
      PyErr_SetString(PyExc_ValueError,
		      "Wrong nside value (must be a power of 2)");
      return NULL;
    }

  /*
  void read_pixwin (const std::string &dir, int nside, arr<double> &temp);
  void read_pixwin (const std::string &dir, int nside, arr<double> &temp,
		    arr<double> &pol);
  */
  arr<double> pw_temp;
  arr<double> pw_pol;

  read_pixwin(datapath, nside, pw_temp, pw_pol);
  
  PyArrayObject *pixwin_temp = NULL;
  PyArrayObject *pixwin_pol = NULL;

  npy_intp szpw;
  
  szpw = pw_temp.size();

  pixwin_temp = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szpw, 
						  PyArray_DOUBLE);
  if( !pixwin_temp )
    return NULL;

  pixwin_pol  = (PyArrayObject*)PyArray_SimpleNew(1, (npy_intp*)&szpw, 
						  PyArray_DOUBLE);
  if( !pixwin_pol )
    return NULL;
  
  for(int i=0; i<szpw; i++ )
    {
      *(double*)PyArray_GETPTR1(pixwin_temp, i) = pw_temp[i];
      *(double*)PyArray_GETPTR1(pixwin_pol, i) = pw_pol[i];
    }

  if( !polarisation )
    {
      Py_DECREF(pixwin_pol);
      return Py_BuildValue("N",pixwin_temp);
    }
  else
    return Py_BuildValue("NN",pixwin_temp,pixwin_pol);

}

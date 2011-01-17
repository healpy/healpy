/*
 *  This file is part of libpsht.
 *
 *  libpsht is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libpsht is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libpsht; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libpsht is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file psht_almhelpers.c
 *  Spherical transform library
 *
 *  Copyright (C) 2008, 2009, 2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "psht_almhelpers.h"
#include "c_utils.h"

void psht_make_triangular_alm_info (int lmax, int mmax, int stride,
  psht_alm_info **alm_info)
  {
  ptrdiff_t m;
  int tval;
  psht_alm_info *info = RALLOC(psht_alm_info,1);
  info->lmax = lmax;
  info->mmax = mmax;
  info->mstart = RALLOC(ptrdiff_t,mmax+1);
  info->stride = stride;
  tval = 2*lmax+1;
  for (m=0; m<=mmax; ++m)
    info->mstart[m] = stride*((m*(tval-m))>>1);
  *alm_info = info;
  }

void psht_make_rectangular_alm_info (int lmax, int mmax, int stride,
  psht_alm_info **alm_info)
  {
  ptrdiff_t m;
  psht_alm_info *info = RALLOC(psht_alm_info,1);
  info->lmax = lmax;
  info->mmax = mmax;
  info->mstart = RALLOC(ptrdiff_t,mmax+1);
  info->stride = stride;
  for (m=0; m<=mmax; ++m)
    info->mstart[m] = stride*m*(lmax+1);
  *alm_info = info;
  }

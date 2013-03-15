/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file sharp_almhelpers.c
 *  Spherical transform library
 *
 *  Copyright (C) 2008-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "sharp_almhelpers.h"
#include "c_utils.h"

void sharp_make_triangular_alm_info (int lmax, int mmax, int stride,
  sharp_alm_info **alm_info)
  {
  sharp_alm_info *info = RALLOC(sharp_alm_info,1);
  info->lmax = lmax;
  info->nm = mmax+1;
  info->mval = RALLOC(int,mmax+1);
  info->mvstart = RALLOC(ptrdiff_t,mmax+1);
  info->stride = stride;
  info->flags = 0;
  ptrdiff_t tval = 2*lmax+1;
  for (ptrdiff_t m=0; m<=mmax; ++m)
    {
    info->mval[m] = m;
    info->mvstart[m] = stride*((m*(tval-m))>>1);
    }
  *alm_info = info;
  }

void sharp_make_rectangular_alm_info (int lmax, int mmax, int stride,
  sharp_alm_info **alm_info)
  {
  sharp_alm_info *info = RALLOC(sharp_alm_info,1);
  info->lmax = lmax;
  info->nm = mmax+1;
  info->mval = RALLOC(int,mmax+1);
  info->mvstart = RALLOC(ptrdiff_t,mmax+1);
  info->stride = stride;
  info->flags = 0;
  for (ptrdiff_t m=0; m<=mmax; ++m)
    {
    info->mval[m] = m;
    info->mvstart[m] = stride*m*(lmax+1);
    }
  *alm_info = info;
  }

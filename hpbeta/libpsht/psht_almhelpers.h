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

/*! \file psht_almhelpers.h
 *  PSHT helper function for the creation of a_lm data structures
 *
 *  Copyright (C) 2008 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_PSHT_ALMHELPERS_H
#define PLANCK_PSHT_ALMHELPERS_H

#include "psht.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Initialises an a_lm data structure according to the scheme used by
    Healpix_cxx.
    \ingroup almgroup */
void psht_make_triangular_alm_info (int lmax, int mmax, int stride,
  psht_alm_info **alm_info);

/*! Initialises an a_lm data structure according to the scheme used by
    Fortran Healpix
    \ingroup almgroup */
void psht_make_rectangular_alm_info (int lmax, int mmax, int stride,
  psht_alm_info **alm_info);

#ifdef __cplusplus
}
#endif

#endif

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

/*! \file psht_geomhelpers.h
 *  PSHT helper function for the creation of grid geometries
 *
 *  Copyright (C) 2006, 2007, 2008 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_PSHT_GEOMHELPERS_H
#define PLANCK_PSHT_GEOMHELPERS_H

#include "psht.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Creates a geometry information describing a HEALPix map with an
    Nside parameter \a nside.
    \ingroup geominfogroup */
void psht_make_healpix_geom_info (int nside, int stride,
  psht_geom_info **geom_info);

/*! Creates a geometry information describing a HEALPix map with an
    Nside parameter \a nside. \a weight contains the relative ring
    weights and must have \a 2*nside entries.
    \ingroup geominfogroup */
void psht_make_weighted_healpix_geom_info (int nside, int stride,
  const double *weight, psht_geom_info **geom_info);

/*! Creates a geometry information describing a Gaussian map with \a nrings
    iso-latitude rings and \a nphi pixels per ring. The azimuth of the first
    pixel in each ring is 0.
    \ingroup geominfogroup */
void psht_make_gauss_geom_info (int nrings, int nphi, int stride,
  psht_geom_info **geom_info);

/*! Creates a geometry information describing an ECP map with \a nrings
    iso-latitude rings and \a nphi pixels per ring. The azimuth of the first
    pixel in each ring is \a phi0 (in radians).
    \ingroup geominfogroup */
void psht_make_ecp_geom_info (int nrings, int nphi, double phi0, int stride,
  psht_geom_info **geom_info);

#ifdef __cplusplus
}
#endif

#endif

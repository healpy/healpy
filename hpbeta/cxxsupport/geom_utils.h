/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file geom_utils.h
 *  Geometric utility functions.
 *
 *  Copyright (C) 2003, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#ifndef PLANCK_GEOM_UTILS_H
#define PLANCK_GEOM_UTILS_H

#include "cxxutils.h"
#include "vec3.h"

/*! Returns the orientation when looking from point \a loc on the unit
    sphere in the direction \a dir. \a loc must be normalized. The result
    ranges from -pi to pi, is 0 for North and pi/2 for West, i.e. the angle
    is given in mathematically positive sense.

    If \a loc is the North or South pole, the returned angle is
    \a atan2(dir.y,dir.x). */
inline double orientation (const vec3 &loc, const vec3 &dir)
  {
// FIXME: here is still optimization potential
  if (loc.x==0 && loc.y==0)
    return (loc.z>0) ? safe_atan2(dir.y,-dir.x) : safe_atan2(dir.y,dir.x);
  vec3 east (-loc.y, loc.x, 0);
  vec3 north = crossprod(loc,east);
  return safe_atan2(-dotprod(dir,east),dotprod(dir,north));
  }

/*! Returns the angle between \a v1 and \a v2 in radians. */
inline double v_angle (const vec3 &v1, const vec3 &v2)
  { return atan2 (crossprod(v1,v2).Length(), dotprod(v1,v2)); }

#endif

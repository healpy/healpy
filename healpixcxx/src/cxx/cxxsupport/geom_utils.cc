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

/*
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "geom_utils.h"
#include "arr.h"

using namespace std;

namespace {

void get_circle (const arr<vec3> &point, tsize q1, tsize q2, vec3 &center,
  double &cosrad)
  {
  center = (point[q1]+point[q2]).Norm();
  cosrad = dotprod(point[q1],center);
  for (tsize i=0; i<q1; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      {
      center=crossprod(point[q1]-point[i],point[q2]-point[i]).Norm();
      cosrad=dotprod(point[i],center);
      if (cosrad<0)
        { center.Flip(); cosrad=-cosrad; }
      }
  }
void get_circle (const arr<vec3> &point, tsize q, vec3 &center,
  double &cosrad)
  {
  center = (point[0]+point[q]).Norm();
  cosrad = dotprod(point[0],center);
  for (tsize i=1; i<q; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      get_circle(point,i,q,center,cosrad);
  }

} // unnamed namespace

void find_enclosing_circle (const arr<vec3> &point, vec3 &center,
  double &cosrad)
  {
  tsize np=point.size();
  planck_assert(np>=3,"too few points");
  center = (point[0]+point[1]).Norm();
  cosrad = dotprod(point[0],center);
  for (tsize i=2; i<np; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      get_circle(point,i,center,cosrad);
  }

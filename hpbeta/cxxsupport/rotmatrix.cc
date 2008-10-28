/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Class for rotation transforms in 3D space
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "rotmatrix.h"
#include "vec3.h"
#include "lsconstants.h"
#include <algorithm>

using namespace std;

rotmatrix::rotmatrix (const vec3 &a, const vec3 &b, const vec3 &c)
  {
  entry[0][0]=a.x; entry[0][1]=b.x; entry[0][2]=c.x;
  entry[1][0]=a.y; entry[1][1]=b.y; entry[1][2]=c.y;
  entry[2][0]=a.z; entry[2][1]=b.z; entry[2][2]=c.z;
  }

void rotmatrix::SetToIdentity ()
  {
  entry[0][0] = entry[1][1] = entry[2][2] = 1.;
  entry[0][1] = entry[1][0] = entry[0][2] =
  entry[2][0] = entry[1][2] = entry[2][1] = 0.;
  }

void rotmatrix::SetToZero ()
  {
  for (int m=0; m<3; ++m)
    entry[m][0] = entry[m][1] = entry[m][2] = 0;
  }

void rotmatrix::Transpose ()
  {
  swap(entry[0][1], entry[1][0]);
  swap(entry[0][2], entry[2][0]);
  swap(entry[1][2], entry[2][1]);
  }

void rotmatrix::toAxisAngle (vec3 &axis, double &angle) const
  {
  double c2 = entry[0][0] + entry[1][1] + entry[2][2] - 1;
  axis.x = entry[2][1] - entry[1][2];
  axis.y = entry[0][2] - entry[2][0];
  axis.z = entry[1][0] - entry[0][1];

  double s2 = axis.Length();

  if (s2>0)
    {
    angle = atan2(s2,c2);
    axis *= 1/s2;
    return;
    }

  if (c2>=2) // angle is 0
    {
    axis = vec3(1,0,0);
    angle = 0;
    return;
    }

  angle = pi;

  int choice = 0; // assume entry[0][0] is the largest
  if ((entry[1][1]>entry[0][0]) && (entry[1][1]>entry[2][2])) choice=1;
  if ((entry[2][2]>entry[0][0]) && (entry[2][2]>entry[1][1])) choice=2;

  if (choice==0)
    {
    axis.x = 0.5*sqrt(entry[0][0]-entry[1][1]-entry[2][2]+1);
    double half_inv = 0.5/axis.x;
    axis.y = half_inv*entry[0][1];
    axis.z = half_inv*entry[0][2];
    return;
    }
  if (choice==1)
    {
    axis.y = 0.5*sqrt(entry[1][1]-entry[0][0]-entry[2][2]+1);
    double half_inv = 0.5/axis.y;
    axis.x = half_inv*entry[0][1];
    axis.z = half_inv*entry[1][2];
    return;
    }

  axis.z = 0.5*sqrt(entry[2][2]-entry[0][0]-entry[1][1]+1);
  double half_inv = 0.5/axis.z;
  axis.x = half_inv*entry[0][2];
  axis.y = half_inv*entry[1][2];
  }

void rotmatrix::Make_CPAC_Euler_Matrix
  (double alpha, double beta, double gamma)
  {
  double ca=cos(alpha), cb=cos(beta), cg=cos(gamma);
  double sa=sin(alpha), sb=sin(beta), sg=sin(gamma);

  entry[0][0]= ca*cb*cg-sa*sg; entry[0][1]=-ca*cb*sg-sa*cg; entry[0][2]= ca*sb;
  entry[1][0]= sa*cb*cg+ca*sg; entry[1][1]=-sa*cb*sg+ca*cg; entry[1][2]= sa*sb;
  entry[2][0]=-sb*cg;          entry[2][1]= sb*sg;          entry[2][2]= cb;
  }

void rotmatrix::Extract_CPAC_Euler_Angles
  (double &alpha, double &beta, double &gamma) const
  {
  double cb = entry[2][2];
  double sb = sqrt(entry[0][2]*entry[0][2] + entry[1][2]*entry[1][2]);
  beta=atan2(sb,cb);
  if (abs(sb)<=1e-6)
    {
    alpha=0;
    if (cb>0)
      gamma=atan2(entry[1][0],entry[0][0]);
    else
      gamma=atan2(entry[0][1],-entry[0][0]);
    }
  else
    {
    alpha=atan2(entry[1][2],entry[0][2]);
    gamma=atan2(entry[2][1],-entry[2][0]);
    }
  }

rotmatrix operator* (const rotmatrix &a, const rotmatrix &b)
  {
  rotmatrix res;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      res.entry[i][j] = a.entry[i][0] * b.entry[0][j]
                      + a.entry[i][1] * b.entry[1][j]
                      + a.entry[i][2] * b.entry[2][j];
  return res;
  }

void TransposeTimes (const rotmatrix &a, const rotmatrix &b, rotmatrix &res)
  {
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      res.entry[i][j] = a.entry[0][i] * b.entry[0][j]
                      + a.entry[1][i] * b.entry[1][j]
                      + a.entry[2][i] * b.entry[2][j];
  }

ostream &operator<< (ostream &os, const rotmatrix &mat)
  {
  for (int i=0;i<3;++i)
    os << mat.entry[i][0] << ' '
       << mat.entry[i][1] << ' '
       << mat.entry[i][2] << endl;
  return os;
  }

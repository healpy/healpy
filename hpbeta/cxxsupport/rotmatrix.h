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

/*! \file rotmatrix.h
 *  Class for rotation transforms in 3D space
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ROTMATRIX_H
#define PLANCK_ROTMATRIX_H

#include <iostream>
#include "cxxutils.h"
#include "vec3.h"

/*! \defgroup rotmatrixgroup Rotation matrices */
/*! \{ */

/*! Class for rotation transforms in 3D space */
class rotmatrix
  {
  public:
    double entry[3][3];

    rotmatrix () {}

    /*! Constructs a rotation matrix from its nine entries */
    rotmatrix (double a00, double a01, double a02,
               double a10, double a11, double a12,
               double a20, double a21, double a22)
      {
      entry[0][0]=a00; entry[0][1]=a01; entry[0][2]=a02;
      entry[1][0]=a10; entry[1][1]=a11; entry[1][2]=a12;
      entry[2][0]=a20; entry[2][1]=a21; entry[2][2]=a22;
      }

    /*! Constructs a rotation matrix so that \a a is the first column,
        \a b is the second column and \a c is the third column.
        \note The vectors \a a, \a b and \a c must form an orthonormal system!
     */
    rotmatrix (const vec3 &a, const vec3 &b, const vec3 &c);

    /*! Sets the matrix to the identity matrix. */
    void SetToIdentity ();
    /*! Sets all matrix elements to zero. */
    void SetToZero ();
    /*! Transposes the matrix. */
    void Transpose ();

    /*! Extracts a unit-length rotation axis \a axis and a rotation angle
        \a angle from the matrix. */
    void toAxisAngle (vec3 &axis, double &angle) const;

    /*! Constructs a matrix which causes a rotation by \a angle around
        \a axis. \a axis must have unit length. */
    void Make_Axis_Rotation_Transform (const vec3 &axis, double angle)
      {
      double sa=sin(angle), ca=cos(angle);
      double ica=1-ca;
      entry[0][0] = axis.x*axis.x*ica + ca;
      entry[1][1] = axis.y*axis.y*ica + ca;
      entry[2][2] = axis.z*axis.z*ica + ca;
      double t1 = axis.x*axis.y*ica, t2 = axis.z*sa;
      entry[1][0] = t1 + t2;
      entry[0][1] = t1 - t2;
      t1 = axis.x*axis.z*ica; t2 = axis.y*sa;
      entry[2][0] = t1 - t2;
      entry[0][2] = t1 + t2;
      t1 = axis.y*axis.z*ica; t2 = axis.x*sa;
      entry[1][2] = t1 - t2;
      entry[2][1] = t1 + t2;
      }

    /*! Creates a rotation matrix \a A, which performs the following operations
        on a vector \a v, when \a Av is calculated:
        -# rotate \a v around the z-axis by \a gamma,
        -# rotate \a v' around the y-axis by \a beta,
        -# rotate \a v'' around the z-axis by \a alpha.

        \note \a alpha, \a beta and \a gamma are given in radians,
              the rotations are right handed.

        \note This transformation rotates the \e vectors, not the coordinate
              axes! */
    void Make_CPAC_Euler_Matrix (double alpha, double beta, double gamma);

    /*! Extracts the Euler angles \a alpha, \a beta and \a gamma from the
        matrix. For their definition see Make_CPAC_Euler_Matrix().

        \note In case of ambiguity \a alpha will be 0. */
    void Extract_CPAC_Euler_Angles
      (double &alpha, double &beta, double &gamma) const;

    /*! Returns the vector \a vec, transformed by the matrix. */
    vec3 Transform (const vec3 &vec) const
      {
      return vec3
        (vec.x*entry[0][0] + vec.y*entry[0][1] + vec.z*entry[0][2],
         vec.x*entry[1][0] + vec.y*entry[1][1] + vec.z*entry[1][2],
         vec.x*entry[2][0] + vec.y*entry[2][1] + vec.z*entry[2][2]);
      }
    /*! Returns the vector \a vec, transformed by the matrix, in \a vec2. */
    void Transform (const vec3 &vec, vec3 &vec2) const
      {
      vec2.x = vec.x*entry[0][0] + vec.y*entry[0][1] + vec.z*entry[0][2];
      vec2.y = vec.x*entry[1][0] + vec.y*entry[1][1] + vec.z*entry[1][2];
      vec2.z = vec.x*entry[2][0] + vec.y*entry[2][1] + vec.z*entry[2][2];
      }
  };

/*! Returns \a a * \a b.
    \relates rotmatrix */
rotmatrix operator* (const rotmatrix &a, const rotmatrix &b);
/*! Returns \a a * \a b in \a res.
    \relates rotmatrix */
inline void matmult (const rotmatrix &a, const rotmatrix &b, rotmatrix &res)
  {
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      res.entry[i][j] = a.entry[i][0] * b.entry[0][j]
                      + a.entry[i][1] * b.entry[1][j]
                      + a.entry[i][2] * b.entry[2][j];
  }

/*! Returns \a a^T * \a b in \a res.
    \relates rotmatrix */
void TransposeTimes (const rotmatrix &a, const rotmatrix &b, rotmatrix &res);

/*! Writes \a mat to \a os.
    \relates rotmatrix */
std::ostream &operator<< (std::ostream &os, const rotmatrix &mat);

/*! \} */

#endif

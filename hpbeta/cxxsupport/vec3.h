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

/*! \file vec3.h
 *  Class representing 3D cartesian vectors
 *
 *  Copyright (C) 2003, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_VEC3_H
#define PLANCK_VEC3_H

#include <cmath>
#include <iostream>
#include "datatypes.h"

/*! \defgroup vec3group 3D vectors */
/*! \{ */

/*! Class representing a 3D cartesian vector. */
template<typename T>class vec3_t
  {
  public:
    T x, /*!< x-coordinate */
      y, /*!< y-coordinate */
      z; /*!< z-coordinate */

    /*! Default constructor. Does not initialize \a x, \a y, and \a z. */
    vec3_t () {}
    /*! Creates a vector with the coordinates \a xc, \a yc, and \a zc. */
    vec3_t (T xc, T yc, T zc)
      : x(xc), y(yc), z(zc) {}
    template<typename T2> explicit vec3_t (const vec3_t<T2> &orig)
      : x(orig.x), y(orig.y), z(orig.z) {}

    /*! Creates a vector from the provided location (z, phi) */
    static vec3_t fromLoc(T _z, T phi) {
        T sinTheta = std::sqrt((1.0-_z)*(1.0+_z));
        return fromLoc(_z, phi, sinTheta, true);
    }
    /*! Creates a vector from the provided location (z, phi, sin(theta)) */
    static vec3_t fromLoc(T _z, T phi, T sinTheta, bool haveSinTheta=true) {
        return haveSinTheta ? vec3_t(sinTheta*std::cos(phi), sinTheta*std::sin(phi), _z) : fromLoc(_z, phi);
    }

    /*! Sets the vector components to \a xc, \a yc, and \a zc. */
    void Set (T xc, T yc, T zc)
      { x=xc; y=yc; z=zc; }
    /*! Creates a unit vector from a z coordinate and an azimuthal angle. */
    void set_z_phi (T z_, T phi)
      {
      using namespace std;
      T sintheta = sqrt((T(1)-z_)*(T(1)+z_));
      x = sintheta*cos(phi);
      y = sintheta*sin(phi);
      z = z_;
      }

    /*! Normalizes the vector to length 1. */
    void Normalize ()
      {
      using namespace std;
      T l = T(1)/sqrt (x*x + y*y + z*z);
      x*=l; y*=l; z*=l;
      }

    vec3_t Norm() const
      {
      vec3_t res(*this);
      res.Normalize();
      return res;
      }

    /*! Returns the length of the vector. */
    T Length () const
      { return sqrt (x*x + y*y + z*z); }

    /*! Returns the squared length of the vector. */
    T SquaredLength () const
      { return (x*x + y*y + z*z); }
    /*! Returns the vector with the signs of all coordinates flipped. */
    const vec3_t operator- () const
      { return vec3_t (-x, -y, -z); }
    /*! Flips the signs of all coordinates. */
    void Flip ()
      { x=-x; y=-y; z=-z; }
    /*! Returns (\a *this + \a vec). */
    const vec3_t operator+ (const vec3_t &vec) const
      { return vec3_t (x+vec.x, y+vec.y, z+vec.z); }
    /*! Adds \a vec to \a *this. */
    vec3_t &operator+= (const vec3_t &vec)
      { x+=vec.x; y+=vec.y; z+=vec.z; return *this; }
    /*! Returns (\a *this - \a vec). */
    const vec3_t operator- (const vec3_t &vec) const
      { return vec3_t (x-vec.x, y-vec.y, z-vec.z); }
    /*! Subtracts \a vec from \a *this. */
    vec3_t &operator-= (const vec3_t &vec)
      { x-=vec.x; y-=vec.y; z-=vec.z; return *this; }
    /*! Returns the vector scaled by \a fact. */
    const vec3_t operator* (T fact) const
      { return vec3_t (x*fact, y*fact, z*fact); }
    /*! Returns the vector scaled by \a 1/fact. */
    const vec3_t operator/ (T fact) const
      { T xfact = T(1)/fact; return vec3_t (x*xfact, y*xfact, z*xfact); }
    /*! Scales the vector by \a fact. */
    vec3_t &operator*= (T fact)
      { x*=fact; y*=fact; z*=fact; return *this; }
  };

/*! Returns the dot product of \a v1 and \a v2.
    \relates vec3_t */
template<typename T> inline T dotprod(const vec3_t<T> &v1, const vec3_t<T> &v2)
  { return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; }

/*! Returns the cross product of \a a and \a b.
    \relates vec3_t */
template<typename T> inline vec3_t<T> crossprod
  (const vec3_t<T> &a, const vec3_t<T> &b)
  { return vec3_t<T>(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }

/*! Writes \a v to \a os.
    \relates vec3_t */
template<typename T> inline std::ostream &operator<<
  (std::ostream &os, const vec3_t<T> &v)
  {
  os << v.x << ", " << v.y << ", " << v.z << std::endl;
  return os;
  }

/*! Specialisation of vec3_t for 64-bit floating point components */
typedef vec3_t<float64> vec3;
/*! Specialisation of vec3_t for 32-bit floating point components */
typedef vec3_t<float32> vec3f;

/*! \} */

#endif

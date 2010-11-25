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

/*! \file safe_cast.h
 *  Numerical cast operator with additional checks that the value is preserved.
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_SAFE_CAST_H
#define PLANCK_SAFE_CAST_H

#include <limits>
#include "error_handling.h"

template<typename T1, typename T2, bool s1, bool s2> struct safe_cast_helper__
  {};

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,true,true>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert(T2(res)==arg, "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,false,false>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert(T2(res)==arg, "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,true,false>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert((res>=0) && (T2(res)==arg),
      "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,false,true>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert((arg>=0) && (T2(res)==arg),
      "safe_cast: value changed during cast");
    return res;
    }
  };

/*! Tries to cast \a arg from its type to a variable of type \c T1.
    If this conversion leads to a change of the actual value (e.g. due to
    overflow or truncation), an exception is thrown. */
template<typename T1, typename T2> inline T1 safe_cast(const T2 &arg)
  {
  return safe_cast_helper__<T1,T2,std::numeric_limits<T1>::is_signed,
    std::numeric_limits<T2>::is_signed>::cast(arg);
  }

#endif

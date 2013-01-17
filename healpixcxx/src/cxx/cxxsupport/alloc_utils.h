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

/*! \file alloc_utils.h
 *  Classes providing raw memory allocation and deallocation support.
 *
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALLOC_UTILS_H
#define PLANCK_ALLOC_UTILS_H

#include <cstdlib>
#include "datatypes.h"

template <typename T> class normalAlloc__
  {
  public:
    T *alloc(tsize sz) const { return (sz>0) ? new T[sz] : 0; }
    void dealloc (T *ptr) const { delete[] ptr; }
  };

template <typename T, int align> class alignAlloc__
  {
  public:
    T *alloc(tsize sz) const
      {
      using namespace std;
      if (sz==0) return 0;
      planck_assert((align&(align-1))==0,"alignment must be power of 2");
      void *res;
/* OSX up to version 10.5 does not define posix_memalign(), but fortunately
   the normal malloc() returns 16 byte aligned storage */
#ifdef __APPLE__
      planck_assert(align<=16, "bad alignment requested");
      res=malloc(sz*sizeof(T));
      planck_assert(res!=0,"error in malloc()");
#else
      planck_assert(posix_memalign(&res,align,sz*sizeof(T))==0,
        "error in posix_memalign()");
#endif
      return static_cast<T *>(res);
      }
    void dealloc(T *ptr) const
      {
      using namespace std;
      free(ptr);
      }
  };

#endif

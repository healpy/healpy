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

/*! \file safe_ptr.h
 *  Pointer wrapper for better exception safety and leakage prevention
 *
 *  Copyright (C) 2005, 2008 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SAFE_PTR_H
#define PLANCK_SAFE_PTR_H

#include "error_handling.h"

template <typename T> class safe_ptr
  {
  private:
    T *p;
    bool set;

    // forbid copying of safe_ptrs, at least until we know how to do it
    safe_ptr (const safe_ptr &) {}
    safe_ptr &operator= (const safe_ptr &) { return *this; }

  public:
    safe_ptr() : p(0), set(false) {}
    safe_ptr (T *p2) : p(p2), set(true) {}
    ~safe_ptr() { delete p; }

    void operator= (T *p2)
      {
      planck_assert (!set, "safe_ptr: already set");
      set = true;
      p=p2;
      }

    void reset()
      {
      delete p;
      p=0;
      set=false;
      }

    operator T*() { return p; }
    operator const T*() const { return p; }
    T *operator->() { return p; }
    const T *operator->() const { return p; }
  };

#endif

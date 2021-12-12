/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2019-2021 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_MISC_UTILS_H
#define DUCC0_MISC_UTILS_H

#include <cstddef>
#include <tuple>

namespace ducc0 {

namespace detail_misc_utils {

using namespace std;

template<typename T> auto calcShare(size_t nshares, size_t myshare,
  const T &begin, const T &end)
  {
  auto nwork = end-begin;
  auto nbase = nwork/nshares;
  auto additional = nwork%nshares;
  auto lo = begin + (myshare*nbase + ((myshare<additional) ? myshare : additional));
  auto hi = lo+nbase+(myshare<additional);
  return make_tuple(lo, hi);
  }

template<typename T> auto calcShare(size_t nshares, size_t myshare, const T &end)
  { return calcShare(nshares, myshare, T(0), end); }

template<typename shp> shp noncritical_shape(const shp &in, size_t elemsz)
  {
  constexpr size_t critstride = 4096; // must be a power of 2
  auto ndim = in.size();
  if (ndim==1) return in;
  shp res(in);
  size_t stride = elemsz;
  for (size_t i=0, xi=ndim-1; i+1<ndim; ++i, --xi)
    {
    size_t tstride = stride*in[xi];
    if ((tstride&(critstride-1))==0)
       res[xi] += 3;
    stride *= res[xi];
    }
  return res;
  }

}

using detail_misc_utils::calcShare;
using detail_misc_utils::noncritical_shape;

}

#endif

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

/*! \file share_utils.h
 *  Various convenience functions for subdividing tasks into chunks
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SHARE_UTILS_H
#define PLANCK_SHARE_UTILS_H

#include "datatypes.h"

/*! Divides the index range [\a glo; \a ghi) into \a nshares approximately
    equal parts, and returns the sub-range [\a lo; \a hi) of the
    part with the number \a myshare (first part has index 0). */
inline void calcShareGeneral (int64 glo, int64 ghi, int64 nshares,
  int64 myshare, int64 &lo, int64 &hi)
  {
  int64 nwork = ghi-glo;
  int64 nbase = nwork/nshares;
  int64 additional = nwork%nshares;
  lo = glo+myshare*nbase + ((myshare<additional) ? myshare : additional);
  hi = lo+nbase+(myshare<additional);
  }

/*! Helper class for dividing a range of work items into chunks of specified
    size. */
class chunkMaker
  {
  private:
    uint64 s_full, s_chunk, offset;

  public:
    /*! Creates an object that produces chunk information for \a s_full_
        work items and a desired chunk size of \a s_chunk_.
        \note Both \a s_chunk_ and \a s_full_ must be larger than 0. */
    chunkMaker (uint64 s_full_, uint64 s_chunk_)
      : s_full(s_full_), s_chunk(s_chunk_), offset(0) {}

    /*! Returns the total number of chunks. */
    uint64 nchunks() const
      { return 1 + (s_full-1)/s_chunk; }

    /*! Returns the start index of the next chunk in \a start, and its size
        in \a size. If all chunks have been processed already, the return
        value is \a false, else \a true. */
    bool getNext (uint64 &start, uint64 &size)
      {
      using namespace std;
      if (offset>=s_full) return false;
      start=offset;
      size=min(s_chunk,s_full-offset);
      offset+=s_chunk;
      return true;
      }
  };

#endif

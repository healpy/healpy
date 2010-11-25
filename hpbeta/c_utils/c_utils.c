/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Convenience functions
 *
 *  Copyright (C) 2008, 2009, 2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include "c_utils.h"

void util_fail_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }
void util_warn_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }

/* This function tries to avoid allocations with a total size close to a high
   power of two (called the "critical stride" here), by adding a few more bytes
   if necssary. This lowers the probability that two arrays differ by a multiple
   of the critical stride in their starting address, which in turn lowers the
   risk of cache line contention. */
static size_t manipsize(size_t sz)
  {
  const size_t critical_stride=4096, cacheline=64, overhead=32;
  if (sz < (critical_stride/2)) return sz;
  if (((sz+overhead)%critical_stride)>(2*cacheline)) return sz;
  return sz+2*cacheline;
  }

#ifdef __SSE__
#include <xmmintrin.h>
void *util_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = _mm_malloc(manipsize(sz),16);
  UTIL_ASSERT(res,"_mm_malloc() failed");
  return res;
  }
void util_free_ (void *ptr)
  { if ((ptr)!=NULL) _mm_free(ptr); }
#else
void *util_malloc_ (size_t sz)
  {
  void *res;
  if (sz==0) return NULL;
  res = malloc(manipsize(sz));
  UTIL_ASSERT(res,"malloc() failed");
  return res;
  }
void util_free_ (void *ptr)
  { if ((ptr)!=NULL) free(ptr); }
#endif

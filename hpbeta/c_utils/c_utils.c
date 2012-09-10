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
 *  Copyright (C) 2008, 2009, 2010, 2011, 2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "c_utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

void util_fail_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }
void util_warn_ (const char *file, int line, const char *func, const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
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

static void OpenMP_status(void)
  {
#ifndef _OPENMP
  printf("OpenMP: not supported by this binary\n");
#else
  int threads = omp_get_max_threads();
  if (threads>1)
    printf("OpenMP active: max. %d threads.\n",threads);
  else
    printf("OpenMP active, but running with 1 thread only.\n");
#endif
  }

static void MPI_status(void)
  {
#ifndef USE_MPI
  printf("MPI: not supported by this binary\n");
#else
  int tasks;
  MPI_Comm_size(MPI_COMM_WORLD,&tasks);
  if (tasks>1)
    printf("MPI active with %d tasks.\n",tasks);
  else
    printf("MPI active, but running with 1 task only.\n");
#endif
  }

static void vec_status(void)
  {
  printf("Vector math: ");
#if(defined(__AVX__))
  printf("AVX\n");
#elif(defined(__SSE2__))
  printf("SSE2\n");
#elif(defined(__SSE__))
  printf("SSE\n");
#else
  printf("not supported by this binary\n");
#endif
  }

void announce_c (const char *name)
  {
  size_t m, nlen=strlen(name);
  printf("\n+-");
  for (m=0; m<nlen; ++m) printf("-");
  printf("-+\n");
  printf("| %s |\n", name);
  printf("+-");
  for (m=0; m<nlen; ++m) printf("-");
  printf("-+\n\n");
  vec_status();
  OpenMP_status();
  MPI_status();
  printf("\n");
  }

void module_startup_c (const char *name, int argc, int argc_expected,
  const char *argv_expected, int verbose)
  {
  if (verbose) announce_c (name);
  if (argc==argc_expected) return;
  if (verbose) fprintf(stderr, "Usage: %s %s\n", name, argv_expected);
  exit(1);
  }

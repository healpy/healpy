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

/*! \file sse_utils.h
 *  SSE/SSE2/SSE3-related functionality
 *
 *  Copyright (C) 2010,2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SSE_UTILS_H
#define PLANCK_SSE_UTILS_H

#if (defined(__SSE__))

#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef __m128 v4sf; /* vector of 4 floats (SSE1) */

typedef union {
  float f[4];
  v4sf v;
} V4SF;

static inline v4sf build_v4sf (float a, float b, float c, float d)
  { return _mm_set_ps(d,c,b,a); }
static inline void read_v4sf (v4sf v, float *a, float *b, float *c, float *d)
  {
  V4SF tmp;
  tmp.v = v;
  if (a) *a=tmp.f[0];
  if (b) *b=tmp.f[1];
  if (c) *c=tmp.f[2];
  if (d) *d=tmp.f[3];
  }

#ifdef __cplusplus
}
#endif

#endif

#if (defined(__SSE2__))

#include <emmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef __m128d v2df; /* vector of 2 doubles (SSE2) */

typedef union {
  double d[2];
  v2df v;
} V2DF;

typedef struct {
  v2df a,b;
} v2df2;
typedef struct {
  V2DF a,b;
} V2DF2;

#define V2DF_SIGNMASK _mm_set1_pd(-0.0)

static inline v2df build_v2df (double a, double b)
  { return _mm_set_pd(b,a); }
static inline void read_v2df (v2df v, double *a, double *b)
  { _mm_store_sd(a,v); _mm_storeh_pd(b,v); }

static inline int v2df_any_gt (v2df a, v2df b)
  {
  return (_mm_movemask_pd(_mm_cmpgt_pd(_mm_andnot_pd(V2DF_SIGNMASK,a),b))!=0);
  }
static inline int v2df_all_ge (v2df a, v2df b)
  {
  return (_mm_movemask_pd(_mm_cmplt_pd(_mm_andnot_pd(V2DF_SIGNMASK,a),b))==0);
  }
static inline V2DF to_V2DF (v2df x)
  { V2DF X; X.v=x; return X; }
static inline V2DF2 to_V2DF2 (v2df2 x)
  { V2DF2 X; X.a.v=x.a; X.b.v=x.b; return X; }
static inline v2df2 to_v2df2 (V2DF2 X)
  { v2df2 x; x.a=X.a.v; x.b=X.b.v; return x; }
static inline v2df2 zero_v2df2(void)
  { v2df2 x; x.a=x.b=_mm_setzero_pd(); return x; }

#ifdef __cplusplus
}
#endif

#endif

#if (defined(__SSE3__))

#include <pmmintrin.h>

#endif

#endif

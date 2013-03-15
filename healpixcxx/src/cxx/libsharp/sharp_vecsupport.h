/*
 *  This file is part of libsharp.
 *
 *  libsharp is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libsharp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libsharp; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libsharp is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*  \file sharp_vecsupport.h
 *  Convenience functions for vector arithmetics
 *
 *  Copyright (C) 2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef SHARP_VECSUPPORT_H
#define SHARP_VECSUPPORT_H

#include <math.h>
#include "sharp_vecutil.h"

typedef double Ts;

#if (VLEN==1)

typedef double Tv;

#define vadd(a,b) ((a)+(b))
#define vaddeq(a,b) ((a)+=(b))
#define vsub(a,b) ((a)-(b))
#define vsubeq(a,b) ((a)-=(b))
#define vmul(a,b) ((a)*(b))
#define vmuleq(a,b) ((a)*=(b))
#define vfmaeq(a,b,c) ((a)+=(b)*(c))
#define vfmseq(a,b,c) ((a)-=(b)*(c))
#define vfmaaeq(a,b,c,d,e) ((a)+=(b)*(c)+(d)*(e))
#define vfmaseq(a,b,c,d,e) ((a)+=(b)*(c)-(d)*(e))
#define vneg(a) (-(a))
#define vload(a) (a)
#define vabs(a) fabs(a)
#define vsqrt(a) sqrt(a)
#define vlt(a,b) (((a)<(b))?1.:0.)
#define vgt(a,b) (((a)>(b))?1.:0.)
#define vge(a,b) (((a)>=(b))?1.:0.)
#define vne(a,b) (((a)!=(b))?1.:0.)
#define vand(a,b) ((((a)*(b))!=0.)?1.:0.)
#define vor(a,b) ((((a)+(b))!=0.)?1.:0.)

static inline Tv vmin (Tv a, Tv b) { return (a<b) ? a : b; }
static inline Tv vmax (Tv a, Tv b) { return (a>b) ? a : b; }

#define vanyTrue(a) ((a)!=0.)
#define vallTrue(a) ((a)!=0.)
#define vblend(m,a,b) (((m)!=0.) ? (a) : (b))
#define vzero 0.
#define vone 1.

#endif

#if (VLEN==2)

#include <emmintrin.h>

#if defined (__SSE3__)
#include <pmmintrin.h>
#endif
#if defined (__SSE4_1__)
#include <smmintrin.h>
#endif

typedef __m128d Tv;

#define vadd(a,b) _mm_add_pd(a,b)
#define vaddeq(a,b) a=_mm_add_pd(a,b)
#define vsub(a,b) _mm_sub_pd(a,b)
#define vsubeq(a,b) a=_mm_sub_pd(a,b)
#define vmul(a,b) _mm_mul_pd(a,b)
#define vmuleq(a,b) a=_mm_mul_pd(a,b)
#define vfmaeq(a,b,c) a=_mm_add_pd(a,_mm_mul_pd(b,c))
#define vfmseq(a,b,c) a=_mm_sub_pd(a,_mm_mul_pd(b,c))
#define vfmaaeq(a,b,c,d,e) \
  a=_mm_add_pd(a,_mm_add_pd(_mm_mul_pd(b,c),_mm_mul_pd(d,e)))
#define vfmaseq(a,b,c,d,e) \
  a=_mm_add_pd(a,_mm_sub_pd(_mm_mul_pd(b,c),_mm_mul_pd(d,e)))
#define vneg(a) _mm_xor_pd(_mm_set1_pd(-0.),a)
#define vload(a) _mm_set1_pd(a)
#define vabs(a) _mm_andnot_pd(_mm_set1_pd(-0.),a)
#define vsqrt(a) _mm_sqrt_pd(a)
#define vlt(a,b) _mm_cmplt_pd(a,b)
#define vgt(a,b) _mm_cmpgt_pd(a,b)
#define vge(a,b) _mm_cmpge_pd(a,b)
#define vne(a,b) _mm_cmpneq_pd(a,b)
#define vand(a,b) _mm_and_pd(a,b)
#define vor(a,b) _mm_or_pd(a,b)
#define vmin(a,b) _mm_min_pd(a,b)
#define vmax(a,b) _mm_max_pd(a,b);
#define vanyTrue(a) (_mm_movemask_pd(a)!=0)
#define vallTrue(a) (_mm_movemask_pd(a)==3)
#if defined(__SSE4_1__)
#define vblend(m,a,b) _mm_blendv_pd(b,a,m)
#else
static inline Tv vblend(Tv m, Tv a, Tv b)
  { return _mm_or_pd(_mm_and_pd(a,m),_mm_andnot_pd(m,b)); }
#endif
#define vzero _mm_setzero_pd()
#define vone _mm_set1_pd(1.)

#endif

#if (VLEN==4)

#include <immintrin.h>
#ifdef __FMA4__
#include <x86intrin.h>
#endif

typedef __m256d Tv;

#define vadd(a,b) _mm256_add_pd(a,b)
#define vaddeq(a,b) a=_mm256_add_pd(a,b)
#define vsub(a,b) _mm256_sub_pd(a,b)
#define vsubeq(a,b) a=_mm256_sub_pd(a,b)
#define vmul(a,b) _mm256_mul_pd(a,b)
#define vmuleq(a,b) a=_mm256_mul_pd(a,b)
#ifdef __FMA4__
#define vfmaeq(a,b,c) a=_mm256_macc_pd(b,c,a)
#define vfmseq(a,b,c) a=_mm256_nmacc_pd(b,c,a)
#define vfmaaeq(a,b,c,d,e) a=_mm256_macc_pd(d,e,_mm256_macc_pd(b,c,a))
#define vfmaseq(a,b,c,d,e) a=_mm256_nmacc_pd(d,e,_mm256_macc_pd(b,c,a))
#else
#define vfmaeq(a,b,c) a=_mm256_add_pd(a,_mm256_mul_pd(b,c))
#define vfmseq(a,b,c) a=_mm256_sub_pd(a,_mm256_mul_pd(b,c))
#define vfmaaeq(a,b,c,d,e) \
  a=_mm256_add_pd(a,_mm256_add_pd(_mm256_mul_pd(b,c),_mm256_mul_pd(d,e)))
#define vfmaseq(a,b,c,d,e) \
  a=_mm256_add_pd(a,_mm256_sub_pd(_mm256_mul_pd(b,c),_mm256_mul_pd(d,e)))
#endif
#define vneg(a) _mm256_xor_pd(_mm256_set1_pd(-0.),a)
#define vload(a) _mm256_set1_pd(a)
#define vabs(a) _mm256_andnot_pd(_mm256_set1_pd(-0.),a)
#define vsqrt(a) _mm256_sqrt_pd(a)
#define vlt(a,b) _mm256_cmp_pd(a,b,_CMP_LT_OQ)
#define vgt(a,b) _mm256_cmp_pd(a,b,_CMP_GT_OQ)
#define vge(a,b) _mm256_cmp_pd(a,b,_CMP_GE_OQ)
#define vne(a,b) _mm256_cmp_pd(a,b,_CMP_NEQ_OQ)
#define vand(a,b) _mm256_and_pd(a,b)
#define vor(a,b) _mm256_or_pd(a,b)
#define vmin(a,b) _mm256_min_pd(a,b)
#define vmax(a,b) _mm256_max_pd(a,b)
#define vanyTrue(a) (_mm256_movemask_pd(a)!=0)
#define vallTrue(a) (_mm256_movemask_pd(a)==15)
#define vblend(m,a,b) _mm256_blendv_pd(b,a,m)
#define vzero _mm256_setzero_pd()
#define vone _mm256_set1_pd(1.)

#endif

#endif

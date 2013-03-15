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

/*  \file sharp_complex_hacks.h
 *  support for converting vector types and complex numbers
 *
 *  Copyright (C) 2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef SHARP_COMPLEX_HACKS_H
#define SHARP_COMPLEX_HACKS_H

#ifdef __cplusplus
#error This header file cannot be included from C++, only from C
#endif

#include <math.h>
#include <complex.h>
#include "sharp_vecsupport.h"

#define UNSAFE_CODE

#if (VLEN==1)

static inline complex double vhsum_cmplx(Tv a, Tv b)
  { return a+_Complex_I*b; }

static inline void vhsum_cmplx2 (Tv a, Tv b, Tv c, Tv d,
  complex double * restrict c1, complex double * restrict c2)
  { *c1 += a+_Complex_I*b; *c2 += c+_Complex_I*d; }

#endif

#if (VLEN==2)

static inline complex double vhsum_cmplx (Tv a, Tv b)
  {
#if defined(__SSE3__)
  Tv tmp = _mm_hadd_pd(a,b);
#else
  Tv tmp = vadd(_mm_shuffle_pd(a,b,_MM_SHUFFLE2(0,1)),
                _mm_shuffle_pd(a,b,_MM_SHUFFLE2(1,0)));
#endif
  union {Tv v; complex double c; } u;
  u.v=tmp; return u.c;
  }

static inline void vhsum_cmplx2 (Tv a, Tv b, Tv c,
  Tv d, complex double * restrict c1, complex double * restrict c2)
  {
#ifdef UNSAFE_CODE
#if defined(__SSE3__)
  vaddeq(*((__m128d *)c1),_mm_hadd_pd(a,b));
  vaddeq(*((__m128d *)c2),_mm_hadd_pd(c,d));
#else
  vaddeq(*((__m128d *)c1),vadd(_mm_shuffle_pd(a,b,_MM_SHUFFLE2(0,1)),
                               _mm_shuffle_pd(a,b,_MM_SHUFFLE2(1,0))));
  vaddeq(*((__m128d *)c2),vadd(_mm_shuffle_pd(c,d,_MM_SHUFFLE2(0,1)),
                               _mm_shuffle_pd(c,d,_MM_SHUFFLE2(1,0))));
#endif
#else
  union {Tv v; complex double c; } u1, u2;
#if defined(__SSE3__)
  u1.v = _mm_hadd_pd(a,b); u2.v=_mm_hadd_pd(c,d);
#else
  u1.v = vadd(_mm_shuffle_pd(a,b,_MM_SHUFFLE2(0,1)),
              _mm_shuffle_pd(a,b,_MM_SHUFFLE2(1,0)));
  u2.v = vadd(_mm_shuffle_pd(c,d,_MM_SHUFFLE2(0,1)),
              _mm_shuffle_pd(c,d,_MM_SHUFFLE2(1,0)));
#endif
  *c1+=u1.c; *c2+=u2.c;
#endif
  }

#endif

#if (VLEN==4)

static inline complex double vhsum_cmplx (Tv a, Tv b)
  {
  Tv tmp=_mm256_hadd_pd(a,b);
  Tv tmp2=_mm256_permute2f128_pd(tmp,tmp,1);
  tmp=_mm256_add_pd(tmp,tmp2);
#ifdef UNSAFE_CODE
  complex double ret;
  *((__m128d *)&ret)=_mm256_extractf128_pd(tmp, 0);
  return ret;
#else
  union {Tv v; complex double c[2]; } u;
  u.v=tmp; return u.c[0];
#endif
  }

static inline void vhsum_cmplx2 (Tv a, Tv b, Tv c, Tv d,
  complex double * restrict c1, complex double * restrict c2)
  {
  Tv tmp1=_mm256_hadd_pd(a,b), tmp2=_mm256_hadd_pd(c,d);
  Tv tmp3=_mm256_permute2f128_pd(tmp1,tmp2,49),
     tmp4=_mm256_permute2f128_pd(tmp1,tmp2,32);
  tmp1=vadd(tmp3,tmp4);
#ifdef UNSAFE_CODE
  *((__m128d *)c1)=_mm_add_pd(*((__m128d *)c1),_mm256_extractf128_pd(tmp1, 0));
  *((__m128d *)c2)=_mm_add_pd(*((__m128d *)c2),_mm256_extractf128_pd(tmp1, 1));
#else
  union {Tv v; complex double c[2]; } u;
  u.v=tmp1;
  *c1+=u.c[0]; *c2+=u.c[1];
#endif
  }

#endif

#endif

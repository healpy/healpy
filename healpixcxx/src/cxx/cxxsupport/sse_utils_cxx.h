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

/*! \file sse_utils_cxx.h
 *  SSE/SSE2/SSE3-related functionality for C++
 *
 *  Copyright (C) 2011, 2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SSE_UTILS_CXX_H
#define PLANCK_SSE_UTILS_CXX_H

template<typename T, int sz> class svec;

#if (defined(__SSE2__))

#include <xmmintrin.h>
#include <emmintrin.h>

template<> class svec<int, 4>
  {
  public:
    typedef int Ts;
    typedef __m128i Tv;
    typedef union { Tv v; Ts d[4]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const Ts &val) : v(_mm_set1_epi32(val)) {}
    svec (const Ts &val1, const Ts &val2, const Ts &val3, const Ts &val4)
      : v(_mm_set_epi32(val4,val3,val2,val1)) {}

    const svec &operator= (const Ts &val)
      { v=_mm_set1_epi32(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    Ts operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, Ts val)
      { Tu u; u.v=v; u.d[p]=val; v=u.v; }

    const svec &operator+= (const svec &b)
      { v=_mm_add_epi32(v,b.v); return *this; }
    const svec &operator-= (const svec &b)
      { v=_mm_sub_epi32(v,b.v); return *this; }
    svec operator+ (const svec &b) const
      { return svec(_mm_add_epi32(v,b.v)); }
    svec operator- (const svec &b) const
      { return svec(_mm_sub_epi32(v,b.v)); }

    const svec &operator&= (const svec &b)
      { v=_mm_and_si128(v,b.v); return *this; }
    const svec &operator|= (const svec &b)
      { v=_mm_or_si128(v,b.v); return *this; }
    const svec &operator^= (const svec &b)
      { v=_mm_xor_si128(v,b.v); return *this; }
    svec operator& (const svec &b) const
      { return svec(_mm_and_si128(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_si128(v,b.v)); }
    svec operator^ (const svec &b) const
      { return svec(_mm_xor_si128(v,b.v)); }
    svec andnot (const svec &b) const
      { return svec(_mm_andnot_si128(v,b.v)); }

    const svec &operator<<= (int b)
      { v=_mm_slli_epi32(v,b); return *this; }
    svec operator<< (int b) const
      { return svec(_mm_slli_epi32(v,b)); }
    const svec &operator>>= (int b)
      { v=_mm_srai_epi32(v,b); return *this; }
    svec operator>> (int b) const
      { return svec(_mm_srai_epi32(v,b)); }

    svec eq (const svec &b) const
      { return svec(_mm_cmpeq_epi32(v,b.v)); }
    svec gt (const svec &b) const
      { return svec(_mm_cmpgt_epi32(v,b.v)); }
    svec lt (const svec &b) const
      { return svec(_mm_cmplt_epi32(v,b.v)); }
  };

typedef svec<int,4> V4si;

#if 0
template<> class svec<long long , 2>
  {
  public:
    typedef long long Ts;
    typedef __m128i Tv;
    typedef union { Tv v; Ts d[2]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const Ts &val) : v(_mm_set1_epi64x(val)) {}
    svec (const Ts &val1, const Ts &val2)
      : v(_mm_set_epi64x(val2,val1)) {}

    const svec &operator= (const Ts &val)
      { v=_mm_set1_epi64x(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    int operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, int val)
      { Tu u; u.v=v; u.d[p]=val; v=u.v; }

    const svec &operator+= (const svec &b)
      { v=_mm_add_epi64(v,b.v); return *this; }
    const svec &operator-= (const svec &b)
      { v=_mm_sub_epi64(v,b.v); return *this; }
    svec operator+ (const svec &b) const
      { return svec(_mm_add_epi64(v,b.v)); }
    svec operator- (const svec &b) const
      { return svec(_mm_sub_epi64(v,b.v)); }

    const svec &operator&= (const svec &b)
      { v=_mm_and_si128(v,b.v); return *this; }
    const svec &operator|= (const svec &b)
      { v=_mm_or_si128(v,b.v); return *this; }
    const svec &operator^= (const svec &b)
      { v=_mm_xor_si128(v,b.v); return *this; }
    svec operator& (const svec &b) const
      { return svec(_mm_and_si128(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_si128(v,b.v)); }
    svec operator^ (const svec &b) const
      { return svec(_mm_xor_si128(v,b.v)); }
    svec andnot (const svec &b) const
      { return svec(_mm_andnot_si128(v,b.v)); }

    const svec &operator<<= (int b)
      { v=_mm_slli_epi64(v,b); return *this; }
    svec operator<< (int b) const
      { return svec(_mm_slli_epi64(v,b)); }
  };

typedef svec<long long,2> V2di;
#endif

template<> class svec<float, 4>
  {
  public:
    typedef float Ts;
    typedef __m128 Tv;
    typedef union { Tv v; Ts d[4]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const Ts &val) : v(_mm_set1_ps(val)) {}
    svec (Ts val1, Ts val2, Ts val3, Ts val4)
      : v(_mm_set_ps(val4,val3,val2,val1)) {}
    explicit svec (const svec<int,4> &b) : v(_mm_cvtepi32_ps(b.v)) {}

    operator svec<int,4>() const
      { return svec<int,4> (_mm_cvtps_epi32(v)); }
    const svec &operator= (const Ts &val)
      { v=_mm_set1_ps(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    Ts operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, Ts val)
      { Tu u; u.v=v; u.d[p]=val; v=u.v; }

    const svec &operator+= (const svec &b)
      { v=_mm_add_ps(v,b.v); return *this; }
    const svec &operator-= (const svec &b)
      { v=_mm_sub_ps(v,b.v); return *this; }
    const svec &operator*= (const svec &b)
      { v=_mm_mul_ps(v,b.v); return *this; }
    const svec &operator/= (const svec &b)
      { v=_mm_div_ps(v,b.v); return *this; }

    svec operator+ (const svec &b) const
      { return svec(_mm_add_ps(v,b.v)); }
    svec operator- (const svec &b) const
      { return svec(_mm_sub_ps(v,b.v)); }
    svec operator* (const svec &b) const
      { return svec(_mm_mul_ps(v,b.v)); }
    svec operator/ (const svec &b) const
      { return svec(_mm_div_ps(v,b.v)); }

    const svec &operator&= (const svec &b)
      { v=_mm_and_ps(v,b.v); return *this; }
    const svec &operator|= (const svec &b)
      { v=_mm_or_ps(v,b.v); return *this; }
    const svec &operator^= (const svec &b)
      { v=_mm_xor_ps(v,b.v); return *this; }
    svec operator& (const svec &b) const
      { return svec(_mm_and_ps(v,b.v)); }
    svec andnot (const svec &b) const
      { return svec(_mm_andnot_ps(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_ps(v,b.v)); }
    svec operator^ (const svec &b) const
      { return svec(_mm_xor_ps(v,b.v)); }

    svec operator- () const
      { return svec(_mm_xor_ps(_mm_set1_ps(-0.),v)); }

    svec eq (const svec &b) const
      { return svec(_mm_cmpeq_ps(v,b.v)); }
    svec neq (const svec &b) const
      { return svec(_mm_cmpneq_ps(v,b.v)); }
    svec lt (const svec &b) const
      { return svec(_mm_cmplt_ps(v,b.v)); }
    svec le (const svec &b) const
      { return svec(_mm_cmple_ps(v,b.v)); }
    svec gt (const svec &b) const
      { return svec(_mm_cmpgt_ps(v,b.v)); }
    svec ge (const svec &b) const
      { return svec(_mm_cmpge_ps(v,b.v)); }

    void writeTo (Ts *val) const
      { _mm_storeu_ps (val, v); }
    void writeTo (Ts &a, Ts &b, Ts &c, Ts &d) const
      { Tu u; u.v=v; a=u.d[0]; b=u.d[1]; c=u.d[2]; d=u.d[3]; }
    void readFrom (const Ts *val)
      { v=_mm_loadu_ps(val); }
    void readFrom (Ts a, Ts b, Ts c, Ts d)
      { v=_mm_set_ps(d,c,b,a); }
  };

typedef svec<float,4> V4sf;

inline V4sf sqrt(const V4sf &v)
  { return V4sf(_mm_sqrt_ps(v.v)); }
inline V4sf abs(const V4sf &v)
  { return V4sf(_mm_andnot_ps(_mm_set1_ps(-0.),v.v)); }
inline V4sf blend(const V4sf &mask, const V4sf &a, const V4sf &b)
  { return (mask&a)|(mask.andnot(b)); }
inline bool any (const V4sf &a)
  { return _mm_movemask_ps(a.v)!=0; }
inline bool all (const V4sf &a)
  { return _mm_movemask_ps(a.v)==15; }
inline bool none (const V4sf &a)
  { return _mm_movemask_ps(a.v)==0; }
inline V4sf min (const V4sf &a, const V4sf &b)
  { return _mm_min_ps(a.v,b.v); }
inline V4sf max (const V4sf &a, const V4sf &b)
  { return _mm_max_ps(a.v,b.v); }

template<> class svec<double, 2>
  {
  public:
    typedef double Ts;
    typedef __m128d Tv;
    typedef union { Tv v; Ts d[2]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const Ts &val) : v(_mm_set1_pd(val)) {}
    svec (const Ts &val1, const Ts &val2)
      : v(_mm_set_pd(val2,val1)) {}
    explicit svec (const svec<int,4> &b) : v(_mm_cvtepi32_pd(b.v)) {}

    operator svec<int,4>() const
      { return svec<int,4> (_mm_cvtpd_epi32(v)); }
    const svec &operator= (const Ts &val)
      { v=_mm_set1_pd(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    Ts operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, Ts val)
      { Tu u; u.v=v; u.d[p]=val; v=u.v; }

    const svec &operator+= (const svec &b)
      { v=_mm_add_pd(v,b.v); return *this; }
    const svec &operator-= (const svec &b)
      { v=_mm_sub_pd(v,b.v); return *this; }
    const svec &operator*= (const svec &b)
      { v=_mm_mul_pd(v,b.v); return *this; }
    const svec &operator/= (const svec &b)
      { v=_mm_div_pd(v,b.v); return *this; }

    svec operator+ (const svec &b) const
      { return svec(_mm_add_pd(v,b.v)); }
    svec operator- (const svec &b) const
      { return svec(_mm_sub_pd(v,b.v)); }
    svec operator* (const svec &b) const
      { return svec(_mm_mul_pd(v,b.v)); }
    svec operator/ (const svec &b) const
      { return svec(_mm_div_pd(v,b.v)); }

    const svec &operator&= (const svec &b)
      { v=_mm_and_pd(v,b.v); return *this; }
    const svec &operator|= (const svec &b)
      { v=_mm_or_pd(v,b.v); return *this; }
    const svec &operator^= (const svec &b)
      { v=_mm_xor_pd(v,b.v); return *this; }
    svec operator& (const svec &b) const
      { return svec(_mm_and_pd(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_pd(v,b.v)); }
    svec operator^ (const svec &b) const
      { return svec(_mm_xor_pd(v,b.v)); }

    svec operator- () const
      { return svec(_mm_xor_pd(_mm_set1_pd(-0.),v)); }

    svec eq (const svec &b) const
      { return svec(_mm_cmpeq_pd(v,b.v)); }
    svec neq (const svec &b) const
      { return svec(_mm_cmpneq_pd(v,b.v)); }
    svec lt (const svec &b) const
      { return svec(_mm_cmplt_pd(v,b.v)); }
    svec le (const svec &b) const
      { return svec(_mm_cmple_pd(v,b.v)); }
    svec gt (const svec &b) const
      { return svec(_mm_cmpgt_pd(v,b.v)); }
    svec ge (const svec &b) const
      { return svec(_mm_cmpge_pd(v,b.v)); }

    void writeTo (Ts *val) const
      { _mm_storeu_pd (val, v); }
    void writeTo (Ts &a, Ts &b) const
      { _mm_store_sd(&a,v); _mm_storeh_pd(&b,v); }
    void readFrom (const Ts *val)
      { v=_mm_loadu_pd(val); }
    void readFrom (const Ts &a, const Ts &b)
      { v=_mm_set_pd(b,a); }
  };

typedef svec<double,2> V2df;

inline V2df sqrt(const V2df &v)
  { return V2df(_mm_sqrt_pd(v.v)); }
inline V2df abs(const V2df &v)
  { return V2df(_mm_andnot_pd(_mm_set1_pd(-0.),v.v)); }
inline V2df blend(const V2df &mask, const V2df &a, const V2df &b)
  { return V2df(_mm_or_pd(_mm_and_pd(a.v,mask.v),_mm_andnot_pd(mask.v,b.v))); }
inline bool any (const V2df &a)
  { return _mm_movemask_pd(a.v)!=0; }
inline bool all (const V2df &a)
  { return _mm_movemask_pd(a.v)==3; }
inline bool none (const V2df &a)
  { return _mm_movemask_pd(a.v)==0; }

template<typename T> inline T vcast(const V4si &a);
template<typename T> inline T vcast(const V4sf &a);
template<typename T> inline T vcast(const V2df &a);

template<> inline V4si vcast (const V4sf &a)
  { return V4si (_mm_castps_si128(a.v)); }
template<> inline V4sf vcast (const V4si &a)
  { return V4sf (_mm_castsi128_ps(a.v)); }
template<> inline V2df vcast (const V4si &a)
  { return V2df (_mm_castsi128_pd(a.v)); }

#endif

#endif

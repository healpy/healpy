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
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SSE_UTILS_CXX_H
#define PLANCK_SSE_UTILS_CXX_H

#ifndef PLANCK_DISABLE_SSE

template<typename T, int sz> class svec;

#if (defined(__SSE__))

#include <xmmintrin.h>

#define PLANCK_HAVE_SSE

template<> class svec<int, 4>
  {
  public:
    typedef __m128i Tv;
    typedef union { Tv v; int d[4]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const int &val) : v(_mm_set1_epi32(val)) {}
    svec (const int &val1, const int &val2, const int &val3, const int &val4)
      : v(_mm_set_epi32(val4,val3,val2,val1)) {}

    const svec &operator= (const int &val)
      { v=_mm_set1_epi32(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    int operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, int val)
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
    svec operator& (const svec &b) const
      { return svec(_mm_and_si128(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_si128(v,b.v)); }

    svec operator<< (int b) const
      { return svec(_mm_slli_epi32(v,b)); }

    svec eq (const svec &b) const
      { return svec(_mm_cmpeq_epi32(v,b.v)); }
  };

typedef svec<int,4> V4si;

inline V4si shuffle(const V4si &a, const V4si &b, int sh)
  { return V4si(reinterpret_cast<__m128i>(_mm_shuffle_ps
    (reinterpret_cast<__m128>(a.v), reinterpret_cast<__m128>(b.v), sh))); }

template<typename T> inline T vcast(const V4si &a)
  { return reinterpret_cast<const T &>(a); }


template<> class svec<float, 4>
  {
  public:
    typedef __m128 Tv;
    typedef union { Tv v; float d[4]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const float &val) : v(_mm_set1_ps(val)) {}
    svec (float val1, float val2, float val3, float val4)
      : v(_mm_set_ps(val4,val3,val2,val1)) {}
    const svec &operator= (const float &val)
      { v=_mm_set1_ps(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    float operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, float val)
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

    void writeTo (float *val) const
      { _mm_storeu_ps (val, v); }
    void writeTo (float &a, float &b, float &c, float &d) const
      { Tu u; u.v=v; a=u.d[0]; b=u.d[1]; c=u.d[2]; d=u.d[3]; }
    void readFrom (const float *val)
      { v=_mm_loadu_ps(val); }
    void readFrom (float a, float b, float c, float d)
      { v=_mm_set_ps(d,c,b,a); }
  };

typedef svec<float,4> V4sf;

inline V4sf sqrt(const V4sf &v)
  { return V4sf(_mm_sqrt_ps(v.v)); }
inline V4sf abs(const V4sf &v)
  { return V4sf(_mm_andnot_ps(_mm_set1_ps(-0.),v.v)); }
inline V4sf blend(const V4sf &mask, const V4sf &a, const V4sf &b)
  { return V4sf(_mm_or_ps(_mm_and_ps(a.v,mask.v),_mm_andnot_ps(mask.v,b.v))); }
inline bool any (const V4sf &a)
  { return _mm_movemask_ps(a.v)!=0; }
inline bool all (const V4sf &a)
  { return _mm_movemask_ps(a.v)==15; }
inline bool none (const V4sf &a)
  { return _mm_movemask_ps(a.v)==0; }

template<typename T> inline T vcast(const V4sf &a)
  { return reinterpret_cast<const T &>(a); }

#if (defined(__SSE2__))

#include <emmintrin.h>

#define PLANCK_HAVE_SSE2

template<> class svec<double, 2>
  {
  public:
    typedef __m128d Tv;
    typedef union { Tv v; double d[2]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const double &val) : v(_mm_set1_pd(val)) {}
    svec (const double &val1, const double &val2)
      : v(_mm_set_pd(val2,val1)) {}
    explicit svec (const svec<int,4> &b) : v(_mm_cvtepi32_pd(b.v)) {}

    operator svec<int,4>() const
      { return svec<int,4> (_mm_cvtpd_epi32(v)); }
    const svec &operator= (const double &val)
      { v=_mm_set1_pd(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    double operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, double val)
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

    void writeTo (double *val) const
      { _mm_storeu_pd (val, v); }
    void writeTo (double &a, double &b) const
      { _mm_store_sd(&a,v); _mm_storeh_pd(&b,v); }
    void readFrom (const double *val)
      { v=_mm_loadu_pd(val); }
    void readFrom (const double &a, const double &b)
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

template<typename T> inline T vcast(const V2df &a)
  { return reinterpret_cast<const T &>(a); }

#endif

#endif

#endif

#endif

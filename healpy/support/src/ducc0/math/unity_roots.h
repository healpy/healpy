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

#ifndef DUCC0_UNITY_ROOTS_H
#define DUCC0_UNITY_ROOTS_H

#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace ducc0 {

namespace detail_unity_roots {

using namespace std;

template<typename T, typename Tc> class UnityRoots
  {
  private:
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    struct cmplx_ { Thigh r, i; };
    size_t N, mask, shift;
    vector<cmplx_> v1, v2;

    static cmplx_ calc(size_t x, size_t n, Thigh ang)
      {
      x<<=3;
      if (x<4*n) // first half
        {
        if (x<2*n) // first quadrant
          {
          if (x<n) return {cos(Thigh(x)*ang), sin(Thigh(x)*ang)};
          return {sin(Thigh(2*n-x)*ang), cos(Thigh(2*n-x)*ang)};
          }
        else // second quadrant
          {
          x-=2*n;
          if (x<n) return {-sin(Thigh(x)*ang), cos(Thigh(x)*ang)};
          return {-cos(Thigh(2*n-x)*ang), sin(Thigh(2*n-x)*ang)};
          }
        }
      else
        {
        x=8*n-x;
        if (x<2*n) // third quadrant
          {
          if (x<n) return {cos(Thigh(x)*ang), -sin(Thigh(x)*ang)};
          return {sin(Thigh(2*n-x)*ang), -cos(Thigh(2*n-x)*ang)};
          }
        else // fourth quadrant
          {
          x-=2*n;
          if (x<n) return {-sin(Thigh(x)*ang), -cos(Thigh(x)*ang)};
          return {-cos(Thigh(2*n-x)*ang), -sin(Thigh(2*n-x)*ang)};
          }
        }
      }
#if 0  // alternative version, similar speed, but maybe a bit more accurate
    static cmplx_ calc2(size_t x, size_t n)
      {
      static constexpr Thigh pi = Thigh(3.141592653589793238462643383279502884197L);
      Thigh n4 = Thigh(n<<2);

      x<<=3;
      if (x<4*n) // first half
        {
        if (x<2*n) // first quadrant
          {
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {cos(ang), sin(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {sin(ang), cos(ang)};
          }
        else // second quadrant
          {
          x-=2*n;
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {-sin(ang), cos(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {-cos(ang), sin(ang)};
          }
        }
      else
        {
        x=8*n-x;
        if (x<2*n) // third quadrant
          {
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {cos(ang), -sin(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {sin(ang), -cos(ang)};
          }
        else // fourth quadrant
          {
          x-=2*n;
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {-sin(ang), -cos(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {-cos(ang), -sin(ang)};
          }
        }
      }
#endif

  public:
    UnityRoots(size_t n)
      : N(n)
      {
      constexpr auto pi = 3.141592653589793238462643383279502884197L;
      Thigh ang = Thigh(0.25L*pi/n);
      size_t nval = (n+2)/2;
      shift = 1;
      while((size_t(1)<<shift)*(size_t(1)<<shift) < nval) ++shift;
      mask = (size_t(1)<<shift)-1;
      v1.resize(mask+1);
      v1[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v1.size(); ++i)
        v1[i]=calc(i,n,ang);
      v2.resize((nval+mask)/(mask+1));
      v2[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v2.size(); ++i)
        v2[i]=calc(i*(mask+1),n,ang);
      }

    size_t size() const { return N; }

    Tc operator[](size_t idx) const
      {
      if (2*idx<=N)
        {
        auto x1=v1[idx&mask], x2=v2[idx>>shift];
        return Tc(T(x1.r*x2.r-x1.i*x2.i), T(x1.r*x2.i+x1.i*x2.r));
        }
      idx = N-idx;
      auto x1=v1[idx&mask], x2=v2[idx>>shift];
      return Tc(T(x1.r*x2.r-x1.i*x2.i), -T(x1.r*x2.i+x1.i*x2.r));
      }
  };

template<typename T, typename Tc> class MultiExp
  {
  private:
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    struct cmplx_ { Thigh r, i; };
    size_t N, mask, shift;
    vector<cmplx_> v1, v2;

  public:
    MultiExp(T ang0, size_t n)
      : N(n)
      {
      Thigh ang = ang0;
      size_t nval = n+2;
      shift = 1;
      while((size_t(1)<<shift)*(size_t(1)<<shift) < nval) ++shift;
      mask = (size_t(1)<<shift)-1;
      v1.resize(mask+1);
      v1[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v1.size(); ++i)
        v1[i] = {cos(i*ang), sin(i*ang)};
      v2.resize((nval+mask)/(mask+1));
      v2[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v2.size(); ++i)
        v2[i] = {cos((i*(mask+1))*ang), sin((i*(mask+1))*ang)};
      }

    size_t size() const { return N; }

    Tc operator[](size_t idx) const
      {
      auto x1=v1[idx&mask], x2=v2[idx>>shift];
      return Tc(T(x1.r*x2.r-x1.i*x2.i), T(x1.r*x2.i+x1.i*x2.r));
      }
  };

}

using detail_unity_roots::UnityRoots;
using detail_unity_roots::MultiExp;

}

#endif

/*
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

/*! \file alm.h
 *  Class for storing spherical harmonic coefficients.
 *
 *  Copyright (C) 2003-2021 Max-Planck-Society
 *  \author Martin Reinecke
 *
 *  For the a_lm rotation code:
 *  Copyright (c) 2018 Richard Mikael Slevinsky and Everett Dawes
 */

#ifndef DUCC0_ALM_H
#define DUCC0_ALM_H

#include <complex>
#include <cmath>
#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>
#include "ducc0/infra/simd.h"
#include "ducc0/infra/threading.h"
#include "ducc0/infra/mav.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_alm {

using namespace std;

/*! Base class for calculating the storage layout of spherical harmonic
    coefficients. */
class Alm_Base
  {
  protected:
    size_t lmax, arrsize;
    vector<size_t> mval;
    vector<ptrdiff_t> mstart;

  public:
    /*! Returns the total number of coefficients for maximum quantum numbers
        \a l and \a m. */
    static size_t Num_Alms (size_t l, size_t m)
      {
      MR_assert(m<=l,"mmax must not be larger than lmax");
      return ((m+1)*(m+2))/2 + (m+1)*(l-m);
      }
    size_t Num_Alms() const
      { return arrsize; }

    Alm_Base (size_t lmax_, const vector<size_t> &mval_,
              const vector<ptrdiff_t> &mstart_)
      : lmax(lmax_), mval(mval_)
      {
      MR_assert(mval.size()>0, "no m indices supplied");
      MR_assert(mstart_.size()==mval.size(), "mval and mstart have different sizes");
      for (size_t i=0; i<mval.size(); ++i)
        {
        MR_assert(mval[i]<=lmax, "m >= lmax");
        if (i>0)
          MR_assert(mval[i]>mval[i-1], "m not strictly ascending");
        }
      mstart.resize(mval.back()+1, -2*lmax);
      arrsize=0;
      for (size_t i=0; i<mval.size(); ++i)
        {
        mstart[mval[i]] = mstart_[i];
        arrsize = size_t(max(ptrdiff_t(arrsize), mstart_[i]+ptrdiff_t(lmax+1)));
        }
      }

    Alm_Base (size_t lmax_, const vector<size_t> &mval_)
      : lmax(lmax_), mval(mval_)
      {
      MR_assert(mval.size()>0, "no m indices supplied");
      for (size_t i=0; i<mval.size(); ++i)
        {
        MR_assert(mval[i]<=lmax, "m >= lmax");
        if (i>0)
          MR_assert(mval[i]>mval[i-1], "m not strictly ascending");
        }
      mstart.resize(mval.back()+1, -2*lmax);
      for (size_t i=0, cnt=0; i<mval.size(); ++i, cnt+=lmax-mval[i]+1)
        mstart[mval[i]] = ptrdiff_t(cnt)-ptrdiff_t(mval[i]);
      arrsize = size_t(mstart.back()+ptrdiff_t(lmax+1));
      }

    /*! Constructs an Alm_Base object with given \a lmax and \a mmax. */
    Alm_Base (size_t lmax_, size_t mmax_)
      : lmax(lmax_), mval(mmax_+1), mstart(mmax_+1)
      {
      ptrdiff_t idx = 0;
      for (size_t m=0; m<=mmax_; ++m)
        {
        mval[m] = m;
        mstart[m] = idx-m;
        idx += lmax-m+1;
        }
      arrsize = Num_Alms(lmax_, mmax_);
      }

    /*! Returns the maximum \a l */
    size_t Lmax() const { return lmax; }
    /*! Returns the maximum \a m */
    size_t Mmax() const { return mval.back(); }

    size_t n_entries() const { return arrsize; }

    /*! Returns an array index for a given m, from which the index of a_lm
        can be obtained by adding l. */
    size_t index_l0 (size_t m) const
      { return mstart[m]; }

    /*! Returns the array index of the specified coefficient. */
    size_t index (size_t l, size_t m) const
      { return index_l0(m) + l; }

    bool conformable(const Alm_Base &other) const
      {
      return (lmax==other.lmax) && (mval==other.mval) && (mstart==other.mstart);
      }
    bool complete() const
      { return mval.size() == lmax+1; }
  };


// the following struct is an adaptation of the algorithms found in
// https://github.com/MikaelSlevinsky/FastTransforms


struct ft_partial_sph_isometry_plan
  {
  static double Gy_index3_a1(int l, int i, int j) // precondition: i+j == 2*l
    {
    if (l+2 <= i && i <= 2*l)
      return (j+1)*(j+2);
    else if (0 <= i && i <= l-1)
      return -(j+1)*(j+2);
    else if (i==j+2) // j+1==l==i-1
      return 2*(j+1)*(j+2);
    else if (i==j)  // i==j==l
      return -2*(j+1)*(j+2);
    return 0;
    }
  static double Gy_index3_a2(int l, int i) // precondition: i+j == 2*l+2
    {
    if (2 <= i && i <= l)
      return -(i-1)*i;
    else if (l+3 <= i && i <= 2*l+2)
      return (i-1)*i;
    return 0;
    }

  static double Y_index_j_eq_i(int l, int i) // l>=0, i>=0, j>=i
    {
    auto r1 = abs(Gy_index3_a1(l, 2*l-i, i));
    auto r3 = abs(Gy_index3_a2(l, 2*l-i+2));
    return (r1+r3)*0.25;
    }
  static double Y_index_j_eq_i_plus_2(int l, int i) // l>=0, i>=0, j>=i
    {
    auto r1 = Gy_index3_a1(l, 2*l-i, i)*Gy_index3_a2(l, 2*l-i);
    return copysign(sqrt(abs(r1)),r1)*0.25;
    }
  static double Z_index(int l, int i, int j)
    {
    return (i==j) ? (j+1)*(2*l+1-j) : 0.0;
    }

  struct ft_symmetric_tridiagonal
    {
    vector<double> a, b;
    int n;

    ft_symmetric_tridiagonal(int N)
      : a(N), b(N-1), n(N) {}
    void resize(int N)
      {
      a.resize(N);
      b.resize(N-1);
      n=N;
      }
    };

  template<bool high_accuracy> class ft_symmetric_tridiagonal_symmetric_eigen
    {
    private:
      vector<double> A, B, C;
      int sign;

    public:
      vector<double> lambda;
      int n;

    private:
      template<typename Tv, size_t N> DUCC0_NOINLINE int eval_helper
        (int jmin, const vector<double> &c, vector<double> &f) const
        {
        constexpr double eps = 0x1p-52;
        constexpr double floatmin = 0x1p-300;

        if (n<1)
          {
          for (int j=jmin; j<n; ++j)
            f[j] = 0.0;
          return n;
          }
        constexpr size_t vlen=Tv::size();
        constexpr size_t step=vlen*N;
// FIXME: All of the juggling with Tvl and Tv can probably go away once
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=99728 is fixed.
        using Tvl = typename Tv::Tv;
        int j=jmin;
        for (; j+int(step)<=n; j+=int(step))
          {
          Tvl vk[N], vkp1[N], nrm[N], X[N], fj[N];
          for (size_t i=0; i<N; ++i)
            {
            vk[i] = Tv(1);
            vkp1[i] = Tv(0);
            nrm[i] = Tv(1);
            X[i] = Tv(&lambda[j+i*Tv::size()], element_aligned_tag());
            fj[i] = Tv(c[n-1]);
            }
          {
          int k=n-1;
          for (; k>2; k-=3)
            {
            Tvl maxnrm = Tv(0);
            for (size_t i=0; i<N; ++i)
              {
              Tvl vkm1, vkm2, vkm3;
              if constexpr(high_accuracy)
                {
                vkm1 = Tvl(Tv(A[k  ]))*((X[i]+B[k  ])*vk[i] - Tvl(Tv(C[k  ]))*vkp1[i]);
                vkm2 = Tvl(Tv(A[k-1]))*((X[i]+B[k-1])*vkm1  - Tvl(Tv(C[k-1]))*vk[i]);
                vkm3 = Tvl(Tv(A[k-2]))*((X[i]+B[k-2])*vkm2  - Tvl(Tv(C[k-2]))*vkm1);
                }
              else
                {
                vkm1 = (Tvl(Tv(A[k  ]))*X[i]+B[k  ])*vk[i] - Tvl(Tv(C[k  ]))*vkp1[i];
                vkm2 = (Tvl(Tv(A[k-1]))*X[i]+B[k-1])*vkm1  - Tvl(Tv(C[k-1]))*vk[i];
                vkm3 = (Tvl(Tv(A[k-2]))*X[i]+B[k-2])*vkm2  - Tvl(Tv(C[k-2]))*vkm1;
                }
              vkp1[i] = vkm2;
              vk[i] = vkm3;
              nrm[i] += vkm1*vkm1 + vkm2*vkm2 + vkm3*vkm3;
              maxnrm = max(Tv(maxnrm), Tv(nrm[i]));
              fj[i] += vkm1*c[k-1] + vkm2*c[k-2] + vkm3*c[k-3];
              }
            if (any_of(Tv(maxnrm) > eps/floatmin))
              for (size_t i=0; i<N; ++i)
                {
                nrm[i] = Tv(1.0)/sqrt(Tv(nrm[i]));
                vkp1[i] *= nrm[i];
                vk[i] *= nrm[i];
                fj[i] *= nrm[i];
                nrm[i] = Tv(1.0);
                }
            }
          for (; k>0; --k)
            {
            Tvl maxnrm = Tv(0);
            for (size_t i=0; i<N; ++i)
              {
              Tvl vkm1;
              if constexpr(high_accuracy)
                vkm1 = Tvl(Tv(A[k]))*((X[i]+B[k])*vk[i] - Tvl(Tv(C[k]))*vkp1[i]);
              else
                vkm1 = (Tvl(Tv(A[k]))*X[i]+B[k])*vk[i] - Tvl(Tv(C[k]))*vkp1[i];
              vkp1[i] = vk[i];
              vk[i] = vkm1;
              nrm[i] += vkm1*vkm1;
              maxnrm = max(Tv(maxnrm), Tv(nrm[i]));
              fj[i] += vkm1*c[k-1];
              }
            if (any_of(Tv(maxnrm) > eps/floatmin))
              for (size_t i=0; i<N; ++i)
                {
                nrm[i] = Tv(1.0)/sqrt(Tv(nrm[i]));
                vkp1[i] *= nrm[i];
                vk[i] *= nrm[i];
                fj[i] *= nrm[i];
                nrm[i] = Tv(1.0);
                }
            }
          }
          for (size_t i=0; i<N; ++i)
            for (size_t q=0; q<vlen; ++q)
              f[j+vlen*i+q] = fj[i][q]*copysign(1.0/sqrt(nrm[i][q]),sign*vk[i][q]);
          }
        return j;
        }

    public:
      ft_symmetric_tridiagonal_symmetric_eigen() {}
      ft_symmetric_tridiagonal_symmetric_eigen(size_t nmax)
        {
        A.reserve(nmax);
        B.reserve(nmax);
        C.reserve(nmax);
        lambda.reserve(nmax); // FIXME: maybe too much
        }

      void Set (const ft_symmetric_tridiagonal &T, const int sign_)
        {
        A.resize(T.n);
        B.resize(T.n);
        C.resize(T.n);
        sign = sign_;
        n = T.n;

        if (n>1)
          {
          A[n-1] = 1/T.b[n-2];
          if constexpr(high_accuracy)
            B[n-1] = -T.a[n-1];
          else
            B[n-1] = -T.a[n-1]/T.b[n-2];
          }
        for (int i=n-2; i>0; i--)
          {
          A[i] = 1/T.b[i-1];
          if constexpr(high_accuracy)
            {
            B[i] = -T.a[i];
            C[i] = T.b[i];
            }
          else
            {
            B[i] = -T.a[i]/T.b[i-1];
            C[i] = T.b[i]/T.b[i-1];
            }
          }
        }

      void eval (const vector<double> &x, vector<double> &y) const
        {
        int j=0;
        if constexpr (vectorizable<double>)
          {
          j = eval_helper<native_simd<double>,4>(j, x, y);
          j = eval_helper<native_simd<double>,2>(j, x, y);
          j = eval_helper<native_simd<double>,1>(j, x, y);
          }
        eval_helper<typename simd_select<double,1>::type,1>(j, x, y);
        }
    };

  ft_symmetric_tridiagonal T;
  ft_symmetric_tridiagonal_symmetric_eigen<true> F11, F21, F12, F22;
  int l;

  ft_partial_sph_isometry_plan(const int lmax)
    : T((lmax+2)/2), F11(lmax/2), F21((lmax+1)/2), F12((lmax+1)/2), F22((lmax+2)/2), l(-1) {}

  void Set(const int l_)
    {
    l = l_;
    int n11 = l/2;
    T.resize(n11);
    for (int i = 0; i < n11; i++)
      T.a[n11-1-i] = Y_index_j_eq_i(l, 2*i+1);
    for (int i = 0; i < n11-1; i++)
      T.b[n11-2-i] = Y_index_j_eq_i_plus_2(l, 2*i+1);
    int sign = (l%4)/2 == 1 ? 1 : -1;
    F11.Set(T, sign);
    F11.lambda.resize(n11);
    for (int i = 0; i < n11; i++)
      F11.lambda[n11-1-i] = Z_index(l, 2*i+1, 2*i+1);

    int n21 = (l+1)/2;
    T.resize(n21);
    for (int i = 0; i < n21; i++)
      T.a[n21-1-i] = Y_index_j_eq_i(l, 2*i);
    for (int i = 0; i < n21-1; i++)
      T.b[n21-2-i] = Y_index_j_eq_i_plus_2(l, 2*i);
    sign = ((l+1)%4)/2 == 1 ? -1 : 1;
    F21.Set(T, sign);
    F21.lambda.resize(n21);
    for (int i = 0; i < n21; i++)
      F21.lambda[i] = Z_index(l, l+1-l%2+2*i, l+1-l%2+2*i);

    int n12 = (l+1)/2;
    T.resize(n12);
    for (int i = 0; i < n12; i++)
      T.a[i] = Y_index_j_eq_i(l, 2*i+l-l%2+1);
    for (int i = 0; i < n12-1; i++)
      T.b[i] = Y_index_j_eq_i_plus_2(l, 2*i+l-l%2+1);
    F12.Set(T, sign);
    F12.lambda.resize(n12);
    for (int i = 0; i < n12; i++)
      F12.lambda[n12-1-i] = Z_index(l, 2*i, 2*i);

    int n22 = (l+2)/2;
    T.resize(n22);
    for (int i = 0; i < n22; i++)
      T.a[i] = Y_index_j_eq_i(l, 2*i+l+l%2);
    for (int i = 0; i < n22-1; i++)
      T.b[i] = Y_index_j_eq_i_plus_2(l, 2*i+l+l%2);
    sign = (l%4)/2 == 1 ? -1 : 1;
    F22.Set(T, sign);
    F22.lambda.resize(n22);
    for (int i = 0; i < n22; i++)
      F22.lambda[i] = Z_index(l, l+l%2+2*i, l+l%2+2*i);
    }
  };


template<typename T> void xchg_yz(const Alm_Base &base, mav<complex<T>,1> &alm,
  size_t nthreads)
  {
  auto lmax = base.Lmax();
  MR_assert(lmax==base.Mmax(), "lmax and mmax must be equal");

  if (lmax>0) // deal with l==1
    {
    auto t = T(-alm(base.index(1,0)).real()/sqrt(2.));
    alm.v(base.index(1,0)).real(T(-alm(base.index(1,1)).imag()*sqrt(2.)));
    alm.v(base.index(1,1)).imag(t);
    }
  if (lmax<=1) return;
  execDynamic(lmax-1,nthreads,1,[&](ducc0::Scheduler &sched)
    {
    vector<double> tin(2*lmax+3), tout(2*lmax+3), tin2(2*lmax+3);
    ft_partial_sph_isometry_plan F(lmax);
    // iterate downwards in l to get the smaller work packages at the end
    while (auto rng=sched.getNext()) for(auto l=lmax-rng.lo; l+rng.hi>lmax; --l)
      {
      F.Set(l);

      int mstart = 1+(l%2);
      for (int i=0; i<F.F11.n; ++i)
        tin[i] = alm(base.index(l,mstart+2*i)).imag();
      F.F11.eval(tin, tout);
      for (int i=0; i<F.F11.n; ++i)
        alm.v(base.index(l,mstart+2*i)).imag(T(tout[i]));

      mstart = l%2;
      for (int i=0; i<F.F22.n; ++i)
        tin[i] = alm(base.index(l,mstart+2*i)).real();
      if (mstart==0)
        tin[0]/=sqrt(2.);
      F.F22.eval(tin, tout);
      if (mstart==0)
        tout[0]*=sqrt(2.);
      for (int i=0; i<F.F22.n; ++i)
        alm.v(base.index(l,mstart+2*i)).real(T(tout[i]));

      mstart = 2-(l%2);
      for (int i=0; i<F.F21.n; ++i)
        tin[i] = alm(base.index(l,mstart+2*i)).imag();

      mstart = 1-(l%2);
      for (int i=0; i<F.F12.n; ++i)
        tin2[i] = alm(base.index(l,mstart+2*i)).real();
      if (mstart==0)
        tin2[0]/=sqrt(2.);
      F.F21.eval(tin,tout);
      if (mstart==0)
        tout[0]*=sqrt(2.);
      for (int i=0; i<F.F12.n; ++i)
        alm.v(base.index(l,mstart+2*i)).real(T(tout[i]));

      F.F12.eval(tin2,tout);
      mstart = 2-(l%2);
      for (int i=0; i<F.F21.n; ++i)
        alm.v(base.index(l,mstart+2*i)).imag(T(tout[i]));
      }
    });
  }

template<typename T> void rotate_alm (const Alm_Base &base, mav<complex<T>,1> &alm,
  double psi, double theta, double phi, size_t nthreads)
  {
  auto lmax=base.Lmax();
  MR_assert (base.complete(), "rotate_alm: need complete A_lm set");
  MR_assert (alm.shape(0)==base.Num_Alms(), "bad size of a_lm array");

  if (theta!=0)
    {
    if (psi!=0)
      for (size_t m=0; m<=lmax; ++m)
        {
        auto exppsi = complex<T>(polar(1.,-psi*m));
        for (size_t l=m; l<=lmax; ++l)
          alm.v(base.index(l,m))*=exppsi;
        }
    xchg_yz(base, alm, nthreads);
    for (size_t m=0; m<=lmax; ++m)
      {
      auto exptheta = complex<T>(polar(1.,-theta*m));
      for (size_t l=m; l<=lmax; ++l)
        alm.v(base.index(l,m))*=exptheta;
      }
    xchg_yz(base, alm, nthreads);
    if (phi!=0)
      for (size_t m=0; m<=lmax; ++m)
        {
        auto expphi = complex<T>(polar(1.,-phi*m));
        for (size_t l=m; l<=lmax; ++l)
          alm.v(base.index(l,m))*=expphi;
        }
    }
  else
    if (phi+psi!=0)
      for (size_t m=0; m<=lmax; ++m)
        {
        auto expang = complex<T>(polar(1.,-(psi+phi)*m));
        for (size_t l=m; l<=lmax; ++l)
          alm.v(base.index(l,m)) *= expang;
        }
  }
}

using detail_alm::Alm_Base;
using detail_alm::rotate_alm;
}

#endif

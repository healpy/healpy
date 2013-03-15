/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm_powspec_tools.h"
#include "string_utils.h"
#include "alm.h"
#include "planck_rng.h"
#include "powspec.h"
#include "xcomplex.h"
#include "rotmatrix.h"
#include "openmp_support.h"
#include "wigner.h"
#include "lsconstants.h"

using namespace std;
template<typename T> void create_alm
  (const PowSpec &powspec, Alm<xcomplex<T> > &alm, planck_rng &rng)
  {
  int lmax = alm.Lmax();
  int mmax = alm.Mmax();
  const double hsqrt2 = 1/sqrt(2.);

  for (int l=0; l<=lmax; ++l)
    {
    double rms_tt = sqrt(powspec.tt(l));
    double zeta1_r = rng.rand_gauss();
    alm(l,0) = T(zeta1_r * rms_tt);
    for (int m=1; m<=min(l,mmax); ++m)
      {
      zeta1_r = rng.rand_gauss()*hsqrt2;
      double zeta1_i = rng.rand_gauss()*hsqrt2;
      alm(l,m).Set (T(zeta1_r*rms_tt), T(zeta1_i*rms_tt));
      }
    }
  }

template void create_alm (const PowSpec &powspec,
  Alm<xcomplex<float> > &alm, planck_rng &rng);
template void create_alm (const PowSpec &powspec,
  Alm<xcomplex<double> > &alm, planck_rng &rng);


template<typename T> void create_alm_pol
  (const PowSpec &powspec,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   planck_rng &rng)
  {
  int lmax = almT.Lmax();
  int mmax = almT.Mmax();
  const double hsqrt2 = 1/sqrt(2.);

  for (int l=0; l<=lmax; ++l)
    {
    double rms_tt=0, rms_g1=0;
    if (powspec.tt(l) != 0)
      {
      rms_tt = sqrt(powspec.tt(l));
      rms_g1 = powspec.tg(l)/rms_tt;
      }

    double zeta1_r = rng.rand_gauss();
    almT(l,0) = T(zeta1_r * rms_tt);
    almG(l,0) = T(zeta1_r * rms_g1);
    for (int m=1; m<=min(l,mmax); ++m)
      {
      zeta1_r = rng.rand_gauss()*hsqrt2;
      double zeta1_i = rng.rand_gauss()*hsqrt2;
      almT(l,m).Set (T(zeta1_r*rms_tt), T(zeta1_i*rms_tt));
      almG(l,m).Set (T(zeta1_r*rms_g1), T(zeta1_i*rms_g1));
      }
    }

  for (int l=0; l<=lmax; ++l)
    {
    double rms_g2 = 0;
    double rms_cc = 0;
    if (powspec.tt(l) != 0)
      {
      rms_g2 = powspec.gg(l) - (powspec.tg(l)/powspec.tt(l))*powspec.tg(l);
      if (rms_g2 <= 0)
        {
        planck_assert (abs(rms_g2) <= 1e-8*abs(powspec.gg(l)),
          "Inconsistent TT, GG and TG spectra at l="+dataToString(l));
        rms_g2 = 0;
        }
      rms_g2 = sqrt(rms_g2);
      rms_cc = sqrt(powspec.cc(l));
      }
    almG(l,0) += T(rng.rand_gauss()*rms_g2);
    almC(l,0)  = T(rng.rand_gauss()*rms_cc);

    for (int m=1; m<=min(l,mmax); ++m)
      {
      double zeta2_r = rng.rand_gauss()*hsqrt2;
      double zeta2_i = rng.rand_gauss()*hsqrt2;
      double zeta3_r = rng.rand_gauss()*hsqrt2;
      double zeta3_i = rng.rand_gauss()*hsqrt2;

      almG(l,m) += xcomplex<T> (T(zeta2_r*rms_g2),T(zeta2_i*rms_g2));
      almC(l,m).Set (T(zeta3_r*rms_cc),T(zeta3_i*rms_cc));
      }
    }
  }

template void create_alm_pol
  (const PowSpec &powspec,
   Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG,
   Alm<xcomplex<float> > &almC,
   planck_rng &rng);
template void create_alm_pol
  (const PowSpec &powspec,
   Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   planck_rng &rng);


template<typename T> void extract_crosspowspec
  (const Alm<xcomplex<T> > &alm1,
   const Alm<xcomplex<T> > &alm2,PowSpec &powspec)
  {
  planck_assert (alm1.conformable(alm2), "a_lm are not conformable");
  arr<double> tt(alm1.Lmax()+1);
  for (int l=0; l<=alm1.Lmax(); ++l)
    {
    tt[l] = alm1(l,0).re*alm2(l,0).re;
    int limit = min(l,alm1.Mmax());
    for (int m=1; m<=limit; ++m)
      tt[l] += 2 * (alm1(l,m).re*alm2(l,m).re + alm1(l,m).im*alm2(l,m).im);
    tt[l] /= (2*l+1);
    }
  powspec.Set(tt);
  }

template void extract_crosspowspec
  (const Alm<xcomplex<float> > &alm1,
   const Alm<xcomplex<float> > &alm2, PowSpec &powspec);
template void extract_crosspowspec
  (const Alm<xcomplex<double> > &alm1,
   const Alm<xcomplex<double> > &alm2, PowSpec &powspec);


template<typename T> void extract_powspec
  (const Alm<xcomplex<T> > &alm, PowSpec &powspec)
  { extract_crosspowspec (alm,alm,powspec); }

template void extract_powspec
  (const Alm<xcomplex<float> > &alm, PowSpec &powspec);
template void extract_powspec
  (const Alm<xcomplex<double> > &alm, PowSpec &powspec);

namespace {

template<typename T> void extract_crosspowspec
  (const Alm<xcomplex<T> > &almT1,
   const Alm<xcomplex<T> > &almG1,
   const Alm<xcomplex<T> > &almC1,
   const Alm<xcomplex<T> > &almT2,
   const Alm<xcomplex<T> > &almG2,
   const Alm<xcomplex<T> > &almC2,
   PowSpec &powspec)
  {
  planck_assert (almT1.conformable(almG1) && almT1.conformable(almC1) &&
                 almT1.conformable(almT2) && almT1.conformable(almG2) &&
                 almT1.conformable(almC2), "a_lm are not conformable");

  int lmax = almT1.Lmax();
  arr<double> tt(lmax+1), gg(lmax+1), cc(lmax+1), tg(lmax+1),
              tc(lmax+1), gc(lmax+1);
  for (int l=0; l<=lmax; ++l)
    {
    tt[l] = almT1(l,0).re*almT2(l,0).re;
    gg[l] = almG1(l,0).re*almG2(l,0).re;
    cc[l] = almC1(l,0).re*almC2(l,0).re;
    tg[l] = almT1(l,0).re*almG2(l,0).re;
    tc[l] = almT1(l,0).re*almC2(l,0).re;
    gc[l] = almG1(l,0).re*almC2(l,0).re;
    int limit = min(l,almT1.Mmax());
    for (int m=1; m<=limit; ++m)
      {
      tt[l] += 2 * (almT1(l,m).re*almT2(l,m).re + almT1(l,m).im*almT2(l,m).im);
      gg[l] += 2 * (almG1(l,m).re*almG2(l,m).re + almG1(l,m).im*almG2(l,m).im);
      cc[l] += 2 * (almC1(l,m).re*almC2(l,m).re + almC1(l,m).im*almC2(l,m).im);
      tg[l] += 2 * (almT1(l,m).re*almG2(l,m).re + almT1(l,m).im*almG2(l,m).im);
      tc[l] += 2 * (almT1(l,m).re*almC2(l,m).re + almT1(l,m).im*almC2(l,m).im);
      gc[l] += 2 * (almG1(l,m).re*almC2(l,m).re + almG1(l,m).im*almC2(l,m).im);
      }
    tt[l] /= (2*l+1);
    gg[l] /= (2*l+1);
    cc[l] /= (2*l+1);
    tg[l] /= (2*l+1);
    tc[l] /= (2*l+1);
    gc[l] /= (2*l+1);
    }
  powspec.Set(tt,gg,cc,tg,tc,gc);
  }

} // unnamed namespace

template<typename T> void extract_powspec
  (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC,
   PowSpec &powspec)
  { extract_crosspowspec(almT,almG,almC,almT,almG,almC,powspec); }

template void extract_powspec
  (const Alm<xcomplex<float> > &almT,
   const Alm<xcomplex<float> > &almG,
   const Alm<xcomplex<float> > &almC,
   PowSpec &powspec);
template void extract_powspec
  (const Alm<xcomplex<double> > &almT,
   const Alm<xcomplex<double> > &almG,
   const Alm<xcomplex<double> > &almC,
   PowSpec &powspec);


template<typename T> void smoothWithGauss
  (Alm<xcomplex<T> > &alm, double fwhm)
  {
  int fct = (fwhm>=0) ? 1 : -1;
  double sigma = fwhm*fwhm2sigma;
  arr<double> gb(alm.Lmax()+1);
  for (int l=0; l<=alm.Lmax(); ++l)
    gb[l] = exp(-.5*fct*l*(l+1)*sigma*sigma);
  alm.ScaleL(gb);
  }

template void smoothWithGauss
  (Alm<xcomplex<float> > &alm, double fwhm);
template void smoothWithGauss
  (Alm<xcomplex<double> > &alm, double fwhm);


template<typename T> void smoothWithGauss
  (Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   double fwhm)
  {
  int fct = (fwhm>=0) ? 1 : -1;
  double sigma = fwhm*fwhm2sigma;
  double fact_pol = exp(2*fct*sigma*sigma);
  arr<double> gb(almT.Lmax()+1);
  for (int l=0; l<=almT.Lmax(); ++l)
    gb[l] = exp(-.5*fct*l*(l+1)*sigma*sigma);
  almT.ScaleL(gb);
  for (int l=0; l<=almT.Lmax(); ++l)
    gb[l] *= fact_pol;
  almG.ScaleL(gb); almC.ScaleL(gb);
  }

template void smoothWithGauss
  (Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG,
   Alm<xcomplex<float> > &almC,
   double fwhm);
template void smoothWithGauss
  (Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   double fwhm);

template<typename T> void rotate_alm (Alm<xcomplex<T> > &alm,
  double psi, double theta, double phi)
  {
  planck_assert (alm.Lmax()==alm.Mmax(),
    "rotate_alm: lmax must be equal to mmax");
  int lmax=alm.Lmax();
  arr<xcomplex<double> > exppsi(lmax+1), expphi(lmax+1);
  for (int m=0; m<=lmax; ++m)
    {
    exppsi[m].Set (cos(psi*m),-sin(psi*m));
    expphi[m].Set (cos(phi*m),-sin(phi*m));
    }

  wigner_d_risbo_openmp rec(lmax,theta);

  arr<xcomplex<double> > almtmp(lmax+1);

  for (int l=0; l<=lmax; ++l)
    {
    const arr2<double> &d(rec.recurse());

    for (int m=0; m<=l; ++m)
      almtmp[m] = xcomplex<double>(alm(l,0))*d[l][l+m];

#pragma omp parallel
{
    int64 lo,hi;
    openmp_calc_share(0,l+1,lo,hi);

    bool flip = true;
    for (int mm=1; mm<=l; ++mm)
      {
      xcomplex<double> t1 = xcomplex<double>(alm(l,mm))*exppsi[mm];
      bool flip2 = ((mm+lo)&1) ? true : false;
      for (int m=lo; m<hi; ++m)
        {
        double d1 = flip2 ? -d[l-mm][l-m] : d[l-mm][l-m];
        double d2 = flip  ? -d[l-mm][l+m] : d[l-mm][l+m];
        double f1 = d1+d2, f2 = d1-d2;
        almtmp[m].re += t1.re*f1; almtmp[m].im += t1.im*f2;
        flip2 = !flip2;
        }
      flip = !flip;
      }
}

    for (int m=0; m<=l; ++m)
      alm(l,m) = xcomplex<T>(almtmp[m]*expphi[m]);
    }
  }

template void rotate_alm (Alm<xcomplex<float> > &alm,
  double psi, double theta, double phi);
template void rotate_alm (Alm<xcomplex<double> > &alm,
  double psi, double theta, double phi);

template<typename T> void rotate_alm (Alm<xcomplex<T> > &almT,
  Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC,
  double psi, double theta, double phi)
  {
  planck_assert (almT.Lmax()==almT.Mmax(),
    "rotate_alm: lmax must be equal to mmax");
  planck_assert (almG.conformable(almT) && almC.conformable(almT),
    "rotate_alm: a_lm are not conformable");
  int lmax=almT.Lmax();
  arr<xcomplex<double> > exppsi(lmax+1), expphi(lmax+1);
  for (int m=0; m<=lmax; ++m)
    {
    exppsi[m].Set (cos(psi*m),-sin(psi*m));
    expphi[m].Set (cos(phi*m),-sin(phi*m));
    }

  wigner_d_risbo_openmp rec(lmax,theta);

  arr<xcomplex<double> > almtmpT(lmax+1), almtmpG(lmax+1), almtmpC(lmax+1);

  for (int l=0; l<=lmax; ++l)
    {
    const arr2<double> &d(rec.recurse());

    for (int m=0; m<=l; ++m)
      {
      almtmpT[m] = xcomplex<double>(almT(l,0))*d[l][m+l];
      almtmpG[m] = xcomplex<double>(almG(l,0))*d[l][m+l];
      almtmpC[m] = xcomplex<double>(almC(l,0))*d[l][m+l];
      }

#pragma omp parallel
{
    int64 lo,hi;
    openmp_calc_share(0,l+1,lo,hi);

    bool flip = true;
    for (int mm=1; mm<=l; ++mm)
      {
      xcomplex<double> t1T = xcomplex<double>(almT(l,mm))*exppsi[mm];
      xcomplex<double> t1G = xcomplex<double>(almG(l,mm))*exppsi[mm];
      xcomplex<double> t1C = xcomplex<double>(almC(l,mm))*exppsi[mm];
      bool flip2 = ((mm+lo)&1) ? true : false;
      for (int m=lo; m<hi; ++m)
        {
        double d1 = flip2 ? -d[l-mm][l-m] : d[l-mm][l-m];
        double d2 = flip  ? -d[l-mm][l+m] : d[l-mm][l+m];
        double f1 = d1+d2, f2 = d1-d2;
        almtmpT[m].re += t1T.re*f1; almtmpT[m].im += t1T.im*f2;
        almtmpG[m].re += t1G.re*f1; almtmpG[m].im += t1G.im*f2;
        almtmpC[m].re += t1C.re*f1; almtmpC[m].im += t1C.im*f2;
        flip2 = !flip2;
        }
      flip = !flip;
      }
}

    for (int m=0; m<=l; ++m)
      {
      almT(l,m) = xcomplex<T>(almtmpT[m]*expphi[m]);
      almG(l,m) = xcomplex<T>(almtmpG[m]*expphi[m]);
      almC(l,m) = xcomplex<T>(almtmpC[m]*expphi[m]);
      }
    }
  }

template void rotate_alm (Alm<xcomplex<float> > &almT,
  Alm<xcomplex<float> > &almG, Alm<xcomplex<float> > &almC,
  double psi, double theta, double phi);
template void rotate_alm (Alm<xcomplex<double> > &almT,
  Alm<xcomplex<double> > &almG, Alm<xcomplex<double> > &almC,
  double psi, double theta, double phi);


template<typename T> void rotate_alm (Alm<xcomplex<T> > &alm,
  const rotmatrix &mat)
  {
  double a1, a2, a3;
  mat.Extract_CPAC_Euler_Angles (a1, a2, a3);
  rotate_alm (alm, a3, a2, a1);
  }

template void rotate_alm (Alm<xcomplex<float> > &alm, const rotmatrix &mat);
template void rotate_alm (Alm<xcomplex<double> > &alm, const rotmatrix &mat);

template<typename T> void rotate_alm (Alm<xcomplex<T> > &almT,
  Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC,
  const rotmatrix &mat)
  {
  double a1, a2, a3;
  mat.Extract_CPAC_Euler_Angles (a1, a2, a3);
  rotate_alm (almT, almG, almC, a3, a2, a1);
  }

template void rotate_alm (Alm<xcomplex<float> > &almT,
  Alm<xcomplex<float> > &almG, Alm<xcomplex<float> > &almC,
  const rotmatrix &mat);
template void rotate_alm (Alm<xcomplex<double> > &almT,
  Alm<xcomplex<double> > &almG, Alm<xcomplex<double> > &almC,
  const rotmatrix &mat);

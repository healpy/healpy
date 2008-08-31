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
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm_powspec_tools.h"
#include "alm.h"
#include "planck_rng.h"
#include "powspec.h"
#include "xcomplex.h"
#include "rotmatrix.h"
#include "openmp_support.h"

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
    alm(l,0) = zeta1_r * rms_tt;
    for (int m=1; m<=min(l,mmax); ++m)
      {
      zeta1_r = rng.rand_gauss()*hsqrt2;
      double zeta1_i = rng.rand_gauss()*hsqrt2;
      alm(l,m).Set (zeta1_r*rms_tt, zeta1_i*rms_tt);
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
    almT(l,0) = zeta1_r * rms_tt;
    almG(l,0) = zeta1_r * rms_g1;
    for (int m=1; m<=min(l,mmax); ++m)
      {
      zeta1_r = rng.rand_gauss()*hsqrt2;
      double zeta1_i = rng.rand_gauss()*hsqrt2;
      almT(l,m).Set (zeta1_r*rms_tt, zeta1_i*rms_tt);
      almG(l,m).Set (zeta1_r*rms_g1, zeta1_i*rms_g1);
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
    almG(l,0) += rng.rand_gauss()*rms_g2;
    almC(l,0)  = rng.rand_gauss()*rms_cc;

    for (int m=1; m<=min(l,mmax); ++m)
      {
      double zeta2_r = rng.rand_gauss()*hsqrt2;
      double zeta2_i = rng.rand_gauss()*hsqrt2;
      double zeta3_r = rng.rand_gauss()*hsqrt2;
      double zeta3_i = rng.rand_gauss()*hsqrt2;

      almG(l,m) += xcomplex<T> (zeta2_r*rms_g2,zeta2_i*rms_g2);
      almC(l,m).Set (zeta3_r*rms_cc,zeta3_i*rms_cc);
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


template<typename T> void extract_powspec
  (const Alm<xcomplex<T> > &alm, PowSpec &powspec)
  {
  arr<double> tt(alm.Lmax()+1);
  for (int l=0; l<=alm.Lmax(); ++l)
    {
    tt[l] = norm(alm(l,0));
    int limit = min(l,alm.Mmax());
    for (int m=1; m<=limit; ++m)
      tt[l] += 2*norm(alm(l,m));
    tt[l] /= (2*l+1);
    }
  powspec.Set(tt);
  }

template void extract_powspec
  (const Alm<xcomplex<float> > &alm, PowSpec &powspec);
template void extract_powspec
  (const Alm<xcomplex<double> > &alm, PowSpec &powspec);


template<typename T> void extract_crosspowspec
  (const Alm<xcomplex<T> > &alm1,
   const Alm<xcomplex<T> > &alm2,PowSpec &powspec)
  {
  planck_assert (alm1.conformable(alm2),
    "extract_crosspowspec: a_lms are not conformable");
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
  (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC,
   PowSpec &powspec)
  {
  planck_assert (almT.conformable(almG) && almT.conformable(almC),
    "extract_powspec: a_lms are not conformable");
  int lmax = almT.Lmax();
  arr<double> tt(lmax+1), gg(lmax+1), cc(lmax+1), tg(lmax+1);
  for (int l=0; l<=lmax; ++l)
    {
    tt[l] = norm(almT(l,0));
    gg[l] = norm(almG(l,0));
    cc[l] = norm(almC(l,0));
    tg[l] = (almT(l,0)*conj(almG(l,0))).re;
    int limit = min(l,almT.Mmax());
    for (int m=1; m<=limit; ++m)
      {
      tt[l] += 2*norm(almT(l,m));
      gg[l] += 2*norm(almG(l,m));
      cc[l] += 2*norm(almC(l,m));
      tg[l] += 2*(almT(l,m)*conj(almG(l,m))).re;
      }
    tt[l] /= (2*l+1);
    gg[l] /= (2*l+1);
    cc[l] /= (2*l+1);
    tg[l] /= (2*l+1);
    }
  powspec.Set(tt,gg,cc,tg);
  }

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


template<typename T> void smooth_with_Gauss
  (Alm<xcomplex<T> > &alm, double fwhm_arcmin)
  {
  int fct = (fwhm_arcmin>=0) ? 1 : -1;
  double sigma = fwhm_arcmin/60*degr2rad*fwhm2sigma;
  arr<double> gb(alm.Lmax()+1);
  for (int l=0; l<=alm.Lmax(); ++l)
    gb[l] = exp(-.5*fct*l*(l+1)*sigma*sigma);
  alm.ScaleL(gb);
  }

template void smooth_with_Gauss
  (Alm<xcomplex<float> > &alm, double fwhm_arcmin);
template void smooth_with_Gauss
  (Alm<xcomplex<double> > &alm, double fwhm_arcmin);


template<typename T> void smooth_with_Gauss
  (Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   double fwhm_arcmin)
  {
  int fct = (fwhm_arcmin>=0) ? 1 : -1;
  double sigma = fwhm_arcmin/60*degr2rad*fwhm2sigma;
  double fact_pol = exp(2*fct*sigma*sigma);
  arr<double> gb(almT.Lmax()+1);
  for (int l=0; l<=almT.Lmax(); ++l)
    gb[l] = exp(-.5*fct*l*(l+1)*sigma*sigma);
  almT.ScaleL(gb);
  for (int l=0; l<=almT.Lmax(); ++l)
    gb[l] *= fact_pol;
  almG.ScaleL(gb); almC.ScaleL(gb);
  }

template void smooth_with_Gauss
  (Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG,
   Alm<xcomplex<float> > &almC,
   double fwhm_arcmin);
template void smooth_with_Gauss
  (Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   double fwhm_arcmin);

namespace {

#if 0

class wigner_d
  {
  private:
    double p,q;
    arr<double> sqt;
    arr2<double> d;
    int n;

    void do_line0 (double *l1, int j)
      {
      double xj = 1./j;
      l1[j] = -p*l1[j-1];
      for (int i=j-1; i>=1; --i)
        l1[i] = xj*sqt[j]*(q*sqt[j-i]*l1[i] - p*sqt[i]*l1[i-1]);
      l1[0] = q*l1[0];
      }
    void do_line (const double *l1, double *l2, int j, int k)
      {
      double xj = 1./j;
      double t1 = xj*sqt[j-k]*q, t2 = xj*sqt[j-k]*p;
      double t3 = xj*sqt[k]*p,   t4 = xj*sqt[k]*q;
      l2[j] = sqt[j] * (t4*l1[j-1]-t2*l2[j-1]);
      for (int i=j-1; i>=1; --i)
        l2[i] = t1*sqt[j-i]*l2[i] - t2*sqt[i]*l2[i-1]
               +t3*sqt[j-i]*l1[i] + t4*sqt[i]*l1[i-1];
      l2[0] = sqt[j] * (t3*l1[0]+t1*l2[0]);
      }

  public:
    wigner_d(int lmax, double ang)
      : p(sin(ang/2)), q(cos(ang/2)), sqt(2*lmax+1), d(lmax+1,2*lmax+1), n(-1)
      { for (int m=0; m<sqt.size(); ++m) sqt[m] = sqrt(double(m)); }

    const arr2<double> &recurse ()
      {
      ++n;
      if (n==0)
        d[0][0] = 1;
      else if (n==1)
        {
        d[0][0] = q*q; d[0][1] = -p*q*sqt[2]; d[0][2] = p*p;
        d[1][0] = -d[0][1]; d[1][1] = q*q-p*p; d[1][2] = d[0][1];
        }
      else
        {
        // padding
        int sign = (n&1)? -1 : 1;
        for (int i=0; i<=2*n-2; ++i)
          {
          d[n][i] = sign*d[n-2][2*n-2-i];
          sign=-sign;
          }
        do_line (d[n-1],d[n],2*n-1,n);
        for (int k=n; k>=2; --k)
          {
          do_line (d[k-2],d[k-1],2*n-1,k-1);
          do_line (d[k-1],d[k],2*n,k);
          }
        do_line0 (d[0],2*n-1);
        do_line (d[0],d[1],2*n,1);
        do_line0 (d[0],2*n);
        }
      return d;
      }
  };

#else

class wigner_d
  {
  private:
    double p,q;
    arr<double> sqt;
    arr2<double> d, dd;
    int n;

  public:
    wigner_d(int lmax, double ang)
      : p(sin(ang/2)), q(cos(ang/2)), sqt(2*lmax+1), d(lmax+1,2*lmax+1),
        dd(lmax+1,2*lmax+1), n(-1)
      { for (int m=0; m<sqt.size(); ++m) sqt[m] = sqrt(double(m)); }

    const arr2<double> &recurse ()
      {
      ++n;
      if (n==0)
        d[0][0] = 1;
      else if (n==1)
        {
        d[0][0] = q*q; d[0][1] = -p*q*sqt[2]; d[0][2] = p*p;
        d[1][0] = -d[0][1]; d[1][1] = q*q-p*p; d[1][2] = d[0][1];
        }
      else
        {
        // padding
        int sign = (n&1)? -1 : 1;
        for (int i=0; i<=2*n-2; ++i)
          {
          d[n][i] = sign*d[n-2][2*n-2-i];
          sign=-sign;
          }
        for (int j=2*n-1; j<=2*n; ++j)
          {
          double xj = 1./j;
          dd[0][0] = q*d[0][0];
          for (int i=1;i<j; ++i)
            dd[0][i] = xj*sqt[j]*(q*sqt[j-i]*d[0][i] - p*sqt[i]*d[0][i-1]);
          dd[0][j] = -p*d[0][j-1];
#pragma omp parallel
{
          int k;
#pragma omp for schedule(static)
          for (k=1; k<=n; ++k)
            {
            double t1 = xj*sqt[j-k]*q, t2 = xj*sqt[j-k]*p;
            double t3 = xj*sqt[k  ]*p, t4 = xj*sqt[k  ]*q;
            dd[k][0] = xj*sqt[j]*(q*sqt[j-k]*d[k][0] + p*sqt[k]*d[k-1][0]);
            for (int i=1; i<j; ++i)
              dd[k][i] = t1*sqt[j-i]*d[k  ][i] - t2*sqt[i]*d[k  ][i-1]
                       + t3*sqt[j-i]*d[k-1][i] + t4*sqt[i]*d[k-1][i-1];
            dd[k][j] = -t2*sqt[j]*d[k][j-1] + t4*sqt[j]*d[k-1][j-1];
            }
}
          dd.swap(d);
          }
        }
      return d;
      }
  };

#endif

} // namespace

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

  wigner_d rec(lmax,theta);

  arr<xcomplex<double> > almtmp(lmax+1);

  for (int l=0; l<=lmax; ++l)
    {
    announce_progress (pow(double(l),3),pow(l-1.,3),pow(double(lmax),3));
    const arr2<double> &d(rec.recurse());

    for (int m=0; m<=l; ++m)
      almtmp[m] = alm(l,0)*d[l][l+m];

#pragma omp parallel
{
    int lo,hi;
    openmp_calc_share(0,l+1,lo,hi);

    bool flip = true;
    for (int mm=1; mm<=l; ++mm)
      {
      xcomplex<double> t1 = alm(l,mm)*exppsi[mm];
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
      alm(l,m) = almtmp[m]*expphi[m];
    }
  end_announce_progress();
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

  wigner_d rec(lmax,theta);

  arr<xcomplex<double> > almtmpT(lmax+1), almtmpG(lmax+1), almtmpC(lmax+1);

  for (int l=0; l<=lmax; ++l)
    {
    announce_progress (pow(double(l),3),pow(l-1.,3),pow(double(lmax),3));
    const arr2<double> &d(rec.recurse());

    for (int m=0; m<=l; ++m)
      {
      almtmpT[m] = almT(l,0)*d[l][m+l];
      almtmpG[m] = almG(l,0)*d[l][m+l];
      almtmpC[m] = almC(l,0)*d[l][m+l];
      }

#pragma omp parallel
{
    int lo,hi;
    openmp_calc_share(0,l+1,lo,hi);

    bool flip = true;
    for (int mm=1; mm<=l; ++mm)
      {
      xcomplex<double> t1T = almT(l,mm)*exppsi[mm];
      xcomplex<double> t1G = almG(l,mm)*exppsi[mm];
      xcomplex<double> t1C = almC(l,mm)*exppsi[mm];
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
      almT(l,m) = almtmpT[m]*expphi[m];
      almG(l,m) = almtmpG[m]*expphi[m];
      almC(l,m) = almtmpC[m]*expphi[m];
      }
    }
  end_announce_progress();
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

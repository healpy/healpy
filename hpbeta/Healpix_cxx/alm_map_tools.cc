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
 *  Copyright (C) 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm_map_tools.h"
#include "alm.h"
#include "fftpack_support.h"
#include "ylmgen.h"
#include "xcomplex.h"

using namespace std;

namespace {

struct info_comparator
  {
  inline bool operator()(const ringinfo &a, const ringinfo &b)
    { return a.sth<b.sth; }
  };

struct pair_comparator
  {
  inline bool operator()(const ringpair &a, const ringpair &b)
    {
    if (a.r1.nph==b.r1.nph)
      return a.r1.phi0<b.r1.phi0;
    return a.r1.nph<b.r1.nph;
    }
  };

void init_lam_fact_1d (int m, arr<double> &lam_fact)
  {
  for (int l=m; l<lam_fact.size(); ++l)
    lam_fact[l] = (l<2) ? 0. : 2*sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
  }

void init_lam_fact_deriv_1d (int m, arr<double> &lam_fact)
  {
  lam_fact[m]=0;
  for (int l=m+1; l<lam_fact.size(); ++l)
    lam_fact[l] = sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
  }

void init_normal_l (arr<double> &normal_l)
  {
  for (int l=0; l<normal_l.size(); ++l)
    normal_l[l] = (l<2) ? 0. : sqrt(1./((l+2.)*(l+1.)*l*(l-1.)));
  }

void get_chunk_info (int nrings, int &nchunks, int &chunksize)
  {
  nchunks = nrings/max(100,nrings/10) + 1;
  chunksize = (nrings+nchunks-1)/nchunks;
  }

class ringhelper
  {
  private:
    double phi0_;
    arr<xcomplex<double> > shiftarr, work;
    rfft plan;
    bool norot;

    void update(int nph, int mmax, double phi0)
      {
      norot = (abs(phi0)<1e-14);
      if (!norot)
        {
        if ((mmax!=shiftarr.size()-1) || (!approx(phi0,phi0_,1e-12)))
          {
          shiftarr.alloc(mmax+1);
          phi0_=phi0;
          for (int m=0; m<=mmax; ++m)
            shiftarr[m].Set (cos(m*phi0),sin(m*phi0));
          }
        }
      if (nph!=plan.size())
        plan.Set(nph);
      if (nph>work.size())
        work.alloc(2*nph);
      }

  public:
    ringhelper() : phi0_(0), norot(true) {}

    template<typename T> void phase2ring (int nph, int mmax, double phi0,
      const xcomplex<double> *phase, T *ring)
      {
      update (nph, mmax, phi0);

      for (int m=1; m<nph; ++m) work[m]=0;
      work[0]=phase[0];

      if (norot)
        for (int m=1; m<=mmax; ++m)
          {
          work[m%nph] += phase[m];
          work[nph-1-((m-1)%nph)] += conj(phase[m]);
          }
      else
        for (int m=1; m<=mmax; ++m)
          {
          xcomplex<double> tmp = phase[m]*shiftarr[m];
          work[m%nph] += tmp;
          work[nph-1-((m-1)%nph)] += conj(tmp);
          }

      plan.backward_c(work);
      for (int m=0; m<nph; ++m) ring[m] = work[m].re;
      }
    template<typename T> void phase2ring (int mmax,
      const xcomplex<double> *phase, const ringinfo &info, T *data)
      {
      if (info.nph>0)
        phase2ring (info.nph, mmax, info.phi0, phase, data+info.ofs);
      }
    template<typename T> void phase2pair (int mmax,
      const xcomplex<double> *phase1, const xcomplex<double> *phase2,
      const ringpair &pair, T *data)
      {
      phase2ring (mmax, phase1, pair.r1, data);
      phase2ring (mmax, phase2, pair.r2, data);
      }

    template<typename T> void ring2phase (int nph, int mmax, double phi0,
      double weight, const T *ring, xcomplex<double> *phase)
      {
      update (nph, mmax, -phi0);
      for (int m=0; m<nph; ++m) work[m] = ring[m]*weight;
      plan.forward_c(work);

      if (norot)
        for (int m=0; m<=mmax; ++m)
          phase[m] = work[m%nph];
      else
        for (int m=0; m<=mmax; ++m)
          phase[m] = work[m%nph]*shiftarr[m];
      }
    template<typename T> void ring2phase (int mmax, const ringinfo &info,
      const T *data, xcomplex<double> *phase)
      {
      if (info.nph>0)
        ring2phase (info.nph, mmax, info.phi0, info.weight, data+info.ofs,
          phase);
      }
    template<typename T> void pair2phase (int mmax, const ringpair &pair,
      const T *data, xcomplex<double> *phase1, xcomplex<double> *phase2)
      {
      ring2phase (mmax, pair.r1, data, phase1);
      ring2phase (mmax, pair.r2, data, phase2);
      }
  };

} // namespace


void info2pair(const vector<ringinfo> &info, vector<ringpair> &pair)
  {
  pair.clear();
  vector<ringinfo> info2=info;
  sort(info2.begin(),info2.end(),info_comparator());
  unsigned int pos=0;
  while (pos<info2.size()-1)
    {
    if (approx(info2[pos].cth,-info2[pos+1].cth,1e-12))
      {
      pair.push_back(ringpair(info2[pos],info2[pos+1]));
      pos+=2;
      }
    else
      {
      pair.push_back(ringpair(info2[pos]));
      ++pos;
      }
    }
  if (pos<info2.size())
    pair.push_back(info2[pos]);

  sort(pair.begin(),pair.end(),pair_comparator());
  }


#define MAP2ALM_MACRO(px) \
  { \
  alm_tmp[l].re += px.re*Ylm[l]; \
  alm_tmp[l].im += px.im*Ylm[l]; \
  ++l; \
  }

template<typename T> void map2alm (const vector<ringpair> &pair,
  const T *map, Alm<xcomplex<T> > &alm, bool add_alm)
  {
  int lmax = alm.Lmax(), mmax = alm.Mmax();

  int nchunks, chunksize;
  get_chunk_info(pair.size(),nchunks,chunksize);
  arr2<xcomplex<double> > phas1(chunksize,mmax+1), phas2(chunksize,mmax+1);

  if (!add_alm) alm.SetToZero();

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,int(pair.size()));

#pragma omp parallel
{
    ringhelper helper;

    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      helper.pair2phase(mmax,pair[ith],map,phas1[ith-llim],phas2[ith-llim]);
} // end of parallel region

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<xcomplex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l) alm_tmp[l].Set(0.,0.);
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(pair[ith+llim].r1.cth,pair[ith+llim].r1.sth,m,Ylm,l);
        if (l<=lmax)
          {
          if (pair[ith+llim].r2.nph>0)
            {
            xcomplex<double> p1 = phas1[ith][m]+phas2[ith][m],
                             p2 = phas1[ith][m]-phas2[ith][m];

            if ((l-m)&1)
              MAP2ALM_MACRO(p2)
            for (;l<lmax;)
              {
              MAP2ALM_MACRO(p1)
              MAP2ALM_MACRO(p2)
              }
            if (l==lmax)
              MAP2ALM_MACRO(p1)
            }
          else
            {
            xcomplex<double> p1 = phas1[ith][m];
            for (;l<=lmax;)
              MAP2ALM_MACRO(p1)
            }
          }
        }
      xcomplex<T> *palm = alm.mstart(m);
      for (int l=m; l<=lmax; ++l)
        { palm[l].re += alm_tmp[l].re; palm[l].im += alm_tmp[l].im; }
      }
} // end of parallel region
    }
  }

template void map2alm (const vector<ringpair> &pair,
  const float *map, Alm<xcomplex<float> > &alm, bool add_alm);
template void map2alm (const vector<ringpair> &pair,
  const double *map, Alm<xcomplex<double> > &alm, bool add_alm);

#define SETUP_LAMBDA \
  const double t1  = lam_lm1m*lam_fact[l]; \
  const double a_w = (l-m2)*two_on_s2 + l*(l-1); \
  const double a_x = twocth*(l-1)*lam_lm; \
  const double lambda_w = a_w*lam_lm - t1*c_on_s2; \
  const double lambda_x = m_on_s2 * (a_x-t1);

#define MAP2ALM_POL_MACRO(Tx,Qx,Qy,Ux,Uy) \
  { \
  double lam_lm1m=lam_lm; \
  lam_lm=Ylm[l]; \
  alm_tmp[l][0].re += Tx.re*lam_lm; alm_tmp[l][0].im += Tx.im*lam_lm; \
  SETUP_LAMBDA \
  alm_tmp[l][1].re += Qx.re*lambda_w - Uy.im*lambda_x; \
  alm_tmp[l][1].im += Qx.im*lambda_w + Uy.re*lambda_x; \
  alm_tmp[l][2].re += Ux.re*lambda_w + Qy.im*lambda_x; \
  alm_tmp[l][2].im += Ux.im*lambda_w - Qy.re*lambda_x; \
  ++l; \
  }

template<typename T> void map2alm_pol
  (const vector<ringpair> &pair, const T *mapT, const T *mapQ, const T *mapU,
   Alm<xcomplex<T> > &almT, Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC,
   bool add_alm)
  {
  planck_assert (almT.conformable(almG) && almT.conformable(almC),
    "map2alm_pol: a_lm are not conformable");

  int lmax = almT.Lmax(), mmax = almT.Mmax();

  arr<double> normal_l (lmax+1);
  init_normal_l (normal_l);

  int nchunks, chunksize;
  get_chunk_info(pair.size(),nchunks,chunksize);

  arr2<xcomplex<double> > phas1T(chunksize,mmax+1),phas2T(chunksize,mmax+1),
                          phas1Q(chunksize,mmax+1),phas2Q(chunksize,mmax+1),
                          phas1U(chunksize,mmax+1),phas2U(chunksize,mmax+1);

  if (!add_alm)
    { almT.SetToZero(); almG.SetToZero(); almC.SetToZero(); }

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,int(pair.size()));

#pragma omp parallel
{
    ringhelper helper;

    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      helper.pair2phase (mmax, pair[ith], mapT,
        phas1T[ith-llim], phas2T[ith-llim]);
      helper.pair2phase (mmax, pair[ith], mapQ,
        phas1Q[ith-llim], phas2Q[ith-llim]);
      helper.pair2phase (mmax, pair[ith], mapU,
        phas1U[ith-llim], phas2U[ith-llim]);
      }
} // end of parallel region

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact(lmax+1);
    arr<xcomplex<double>[3] > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      init_lam_fact_1d (m,lam_fact);

      for (int l=m; l<alm_tmp.size(); ++l)
        alm_tmp[l][0]=alm_tmp[l][1]=alm_tmp[l][2] = 0;

      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        double cth=pair[ith+llim].r1.cth, sth=pair[ith+llim].r1.sth;
        generator.get_Ylm(cth,sth,m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s2 = 1/(sth*sth);
          double c_on_s2 = cth * one_on_s2;
          double two_on_s2 = 2*one_on_s2;
          double twocth = 2*cth;
          int m2 = m*m;
          double m_on_s2 = m*one_on_s2;

          if (pair[ith+llim].r2.nph>0)
            {
            xcomplex<double> T1 = phas1T[ith][m]+phas2T[ith][m],
                             T2 = phas1T[ith][m]-phas2T[ith][m],
                             Q1 = phas1Q[ith][m]+phas2Q[ith][m],
                             Q2 = phas1Q[ith][m]-phas2Q[ith][m],
                             U1 = phas1U[ith][m]+phas2U[ith][m],
                             U2 = phas1U[ith][m]-phas2U[ith][m];

            double lam_lm = 0;
            if ((l-m)&1)
              MAP2ALM_POL_MACRO(T2,Q2,Q1,U2,U1)
            for (;l<lmax;)
              {
              MAP2ALM_POL_MACRO(T1,Q1,Q2,U1,U2)
              MAP2ALM_POL_MACRO(T2,Q2,Q1,U2,U1)
              }
            if (l==lmax)
              MAP2ALM_POL_MACRO(T1,Q1,Q2,U1,U2)
            }
          else
            {
            xcomplex<double> T1 = phas1T[ith][m],
                             Q1 = phas1Q[ith][m],
                             U1 = phas1U[ith][m];
            double lam_lm = 0;
            for (;l<=lmax;)
              MAP2ALM_POL_MACRO(T1,Q1,Q1,U1,U1)
            }
          }
        }
      xcomplex<T> *palmT=almT.mstart(m), *palmG=almG.mstart(m),
                  *palmC=almC.mstart(m);
      for (int l=m;l<=lmax;++l)
        {
        palmT[l].re += alm_tmp[l][0].re;
        palmT[l].im += alm_tmp[l][0].im;
        palmG[l].re += alm_tmp[l][1].re*normal_l[l];
        palmG[l].im += alm_tmp[l][1].im*normal_l[l];
        palmC[l].re += alm_tmp[l][2].re*normal_l[l];
        palmC[l].im += alm_tmp[l][2].im*normal_l[l];
        }
      }
} // end of parallel region
    }
  }

template void map2alm_pol (const vector<ringpair> &pair, const float *mapT,
   const float *mapQ, const float *mapU, Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG, Alm<xcomplex<float> > &almC, bool add_alm);
template void map2alm_pol (const vector<ringpair> &pair, const double *mapT,
   const double *mapQ, const double *mapU, Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG, Alm<xcomplex<double> > &almC, bool add_alm);


#define ALM2MAP_MACRO(px) \
  { \
  px.re += alm_tmp[l].re*Ylm[l]; \
  px.im += alm_tmp[l].im*Ylm[l]; \
  ++l; \
  }

template<typename T> void alm2map (const Alm<xcomplex<T> > &alm,
  const vector<ringpair> &pair, T *map)
  {
  int lmax = alm.Lmax(), mmax = alm.Mmax();

  int nchunks, chunksize;
  get_chunk_info(pair.size(),nchunks,chunksize);

  arr2<xcomplex<double> > phas1(chunksize,mmax+1), phas2(chunksize,mmax+1);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,int(pair.size()));

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<xcomplex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l)
        alm_tmp[l]=alm(l,m);

      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(pair[ith+llim].r1.cth,pair[ith+llim].r1.sth,m,Ylm,l);
        if (l>lmax)
          phas1[ith][m] = phas2[ith][m] = 0;
        else
          {
          if (pair[ith+llim].r2.nph>0)
            {
            xcomplex<double> p1=0, p2=0;

            if ((l-m)&1)
              ALM2MAP_MACRO(p2)
            for (;l<lmax;)
              {
              ALM2MAP_MACRO(p1)
              ALM2MAP_MACRO(p2)
              }
            if (l==lmax)
              ALM2MAP_MACRO(p1)
            phas1[ith][m] = p1+p2; phas2[ith][m] = p1-p2;
            }
          else
            {
            xcomplex<double> p1=0;
            for (;l<=lmax;)
              ALM2MAP_MACRO(p1)
            phas1[ith][m] = p1;
            }
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    ringhelper helper;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      helper.phase2pair (mmax,phas1[ith-llim],phas2[ith-llim],pair[ith],map);
} // end of parallel region
    }
  }

template void alm2map (const Alm<xcomplex<float> > &alm,
  const vector<ringpair> &pair, float *map);
template void alm2map (const Alm<xcomplex<double> > &alm,
  const vector<ringpair> &pair, double *map);

#define ALM2MAP_POL_MACRO(Tx,Qx,Qy,Ux,Uy) \
  { \
  double lam_lm1m = lam_lm; \
  lam_lm = Ylm[l]; \
  Tx.re+=alm_tmp[l][0].re*lam_lm;   Tx.im+=alm_tmp[l][0].im*lam_lm; \
  SETUP_LAMBDA \
  Qx.re+=alm_tmp[l][1].re*lambda_w; Qx.im+=alm_tmp[l][1].im*lambda_w; \
  Ux.re-=alm_tmp[l][2].re*lambda_w; Ux.im-=alm_tmp[l][2].im*lambda_w; \
  Qy.re-=alm_tmp[l][2].im*lambda_x; Qy.im+=alm_tmp[l][2].re*lambda_x; \
  Uy.re-=alm_tmp[l][1].im*lambda_x; Uy.im+=alm_tmp[l][1].re*lambda_x; \
  ++l; \
  }

template<typename T> void alm2map_pol (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG, const Alm<xcomplex<T> > &almC,
   const vector<ringpair> &pair, T *mapT, T *mapQ, T *mapU)
  {
  int lmax = almT.Lmax(), mmax = almT.Mmax();

  planck_assert (almT.conformable(almG) && almT.conformable(almC),
    "alm2map_pol: a_lm are not conformable");

  arr<double> normal_l (lmax+1);
  init_normal_l (normal_l);

  int nchunks, chunksize;
  get_chunk_info(pair.size(),nchunks,chunksize);

  arr2<xcomplex<double> >
    phas1T(chunksize,mmax+1), phas2T(chunksize,mmax+1),
    phas1Q(chunksize,mmax+1), phas2Q(chunksize,mmax+1),
    phas1U(chunksize,mmax+1), phas2U(chunksize,mmax+1);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,int(pair.size()));

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact (lmax+1);
    arr<xcomplex<double>[3]> alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      int m2 = m*m;
      init_lam_fact_1d (m,lam_fact);
      for (int l=m; l<=lmax; ++l)
        {
        alm_tmp[l][0] = almT(l,m);
        alm_tmp[l][1] = almG(l,m)*(-normal_l[l]);
        alm_tmp[l][2] = almC(l,m)*(-normal_l[l]);
        }
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        double cth=pair[ith+llim].r1.cth, sth=pair[ith+llim].r1.sth;
        int l;
        generator.get_Ylm(cth,sth,m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s2 = 1/(sth*sth);
          double c_on_s2 = cth * one_on_s2;
          double two_on_s2 = 2*one_on_s2;
          double m_on_s2 = m*one_on_s2;
          double twocth = 2*cth;

          if (pair[ith+llim].r2.nph>0)
            {
            xcomplex<double> T1=0, T2=0, Q1=0, Q2=0, U1=0, U2=0;
            double lam_lm = 0;
            if ((l-m)&1)
              ALM2MAP_POL_MACRO(T2,Q2,Q1,U2,U1)
            for (;l<lmax;)
              {
              ALM2MAP_POL_MACRO(T1,Q1,Q2,U1,U2)
              ALM2MAP_POL_MACRO(T2,Q2,Q1,U2,U1)
              }
            if (l==lmax)
              ALM2MAP_POL_MACRO(T1,Q1,Q2,U1,U2)

            phas1T[ith][m] = T1+T2; phas2T[ith][m] = T1-T2;
            phas1Q[ith][m] =-Q1-Q2; phas2Q[ith][m] =-Q1+Q2;
            phas1U[ith][m] = U1+U2; phas2U[ith][m] = U1-U2;
            }
          else
            {
            xcomplex<double> T1=0, Q1=0, U1=0;
            double lam_lm = 0;
            for (;l<=lmax;)
              { ALM2MAP_POL_MACRO(T1,Q1,Q1,U1,U1) }
            phas1T[ith][m] = T1;
            phas1Q[ith][m] =-Q1;
            phas1U[ith][m] = U1;
            }
          }
        else
          {
          phas1T[ith][m] = phas2T[ith][m] = 0;
          phas1Q[ith][m] = phas2Q[ith][m] = 0;
          phas1U[ith][m] = phas2U[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    ringhelper helper;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      helper.phase2pair(mmax,phas1T[ith-llim],phas2T[ith-llim],pair[ith],mapT);
      helper.phase2pair(mmax,phas1Q[ith-llim],phas2Q[ith-llim],pair[ith],mapQ);
      helper.phase2pair(mmax,phas1U[ith-llim],phas2U[ith-llim],pair[ith],mapU);
      }
} // end of parallel region
    }
  }

template void alm2map_pol (const Alm<xcomplex<float> > &almT,
   const Alm<xcomplex<float> > &almG, const Alm<xcomplex<float> > &almC,
   const vector<ringpair> &pair, float *mapT, float *mapQ, float *mapU);
template void alm2map_pol (const Alm<xcomplex<double> > &almT,
   const Alm<xcomplex<double> > &almG, const Alm<xcomplex<double> > &almC,
   const vector<ringpair> &pair, double *mapT, double *mapQ, double *mapU);

#define ALM2MAP_DER1_MACRO(px,pdthx,pdphx) \
  { \
  double lam_lm1m = lam_lm; \
  lam_lm = Ylm[l]; \
  const double t1 = alm_tmp[l].re*lam_lm; \
  const double t2 = alm_tmp[l].im*lam_lm; \
  const double t3 = l*cotanth; \
  const double t4 = one_on_s*lam_lm1m*lam_fact[l]; \
  px.re+=t1; px.im+=t2; \
  pdthx.re+=t3*t1-t4*alm_tmp[l].re; pdthx.im+=t3*t2-t4*alm_tmp[l].im; \
  pdphx.re-=m*t2; pdphx.im+=m*t1; \
  ++l; \
  }

template<typename T> void alm2map_der1 (const Alm<xcomplex<T> > &alm,
   const vector<ringpair> &pair, T *map, T *dth, T *dph)
  {
  int lmax = alm.Lmax(), mmax = alm.Mmax();

  int nchunks, chunksize;
  get_chunk_info(pair.size(),nchunks,chunksize);

  arr2<xcomplex<double> >
    phas1(chunksize,mmax+1), phas2(chunksize,mmax+1),
    phas1dth(chunksize,mmax+1), phas2dth(chunksize,mmax+1),
    phas1dph(chunksize,mmax+1), phas2dph(chunksize,mmax+1);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,int(pair.size()));

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact (lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      const xcomplex<T> *alm_tmp=alm.mstart(m);
      init_lam_fact_deriv_1d (m,lam_fact);
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        double cth=pair[ith+llim].r1.cth, sth=pair[ith+llim].r1.sth;
        int l;
        generator.get_Ylm(cth,sth,m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s = 1/sth;
          double cotanth = cth*one_on_s;

          if (pair[ith+llim].r2.nph>0)
            {
            xcomplex<double> p1=0, p2=0, pth1=0, pth2=0, pph1=0, pph2=0;

            double lam_lm = 0;
            if ((l-m)&1)
              ALM2MAP_DER1_MACRO(p2,pth2,pph2)
            for(;l<lmax;)
              {
              ALM2MAP_DER1_MACRO(p1,pth1,pph1)
              ALM2MAP_DER1_MACRO(p2,pth2,pph2)
              }
            if (l==lmax)
              ALM2MAP_DER1_MACRO(p1,pth1,pph1)

            phas1[ith][m] = p1+p2; phas2[ith][m] = p1-p2;
            phas1dth[ith][m] = pth1+pth2; phas2dth[ith][m] = pth2-pth1;
            phas1dph[ith][m] = (pph1+pph2)*one_on_s;
            phas2dph[ith][m] = (pph1-pph2)*one_on_s;
            }
          else
            {
            xcomplex<double> p1=0, pth1=0, pph1=0;

            double lam_lm = 0;
            for(;l<=lmax;)
              ALM2MAP_DER1_MACRO(p1,pth1,pph1)

            phas1[ith][m] = p1;
            phas1dth[ith][m] = pth1;
            phas1dph[ith][m] = pph1*one_on_s;
            }
          }
        else
          {
          phas1[ith][m] = phas2[ith][m] = 0;
          phas1dth[ith][m] = phas2dth[ith][m] = 0;
          phas1dph[ith][m] = phas2dph[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    ringhelper helper;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      helper.phase2pair(mmax,phas1[ith-llim],phas2[ith-llim],pair[ith],map);
      helper.phase2pair(mmax,phas1dth[ith-llim],phas2dth[ith-llim],
        pair[ith],dth);
      helper.phase2pair(mmax,phas1dph[ith-llim],phas2dph[ith-llim],
        pair[ith],dph);
      }
} // end of parallel region
    }
  }

template void alm2map_der1 (const Alm<xcomplex<float> > &alm,
   const vector<ringpair> &pair, float *map, float *dth, float *dph);
template void alm2map_der1 (const Alm<xcomplex<double> > &alm,
   const vector<ringpair> &pair, double *map, double *dth, double *dph);

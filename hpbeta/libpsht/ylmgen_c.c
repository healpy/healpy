/*
 *  This file is part of libpsht.
 *
 *  libpsht is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libpsht is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libpsht; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libpsht is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Code for efficient calculation of Y_lm(theta,phi=0)
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "ylmgen_c.h"
#include "c_utils.h"

enum { large_exponent2=90, minscale=-4, maxscale=11, max_spin=100 };

static void sylmgen_init (sylmgen_d *gen, const Ylmgen_C *ygen, int spin)
  {
  int i;
  UTIL_ASSERT(spin>=1,"incorrect spin");
  gen->s=spin;
  gen->m=gen->mlo=gen->mhi=-1234567890;
  ALLOC(gen->fx,ylmgen_dbl3,ygen->lmax+2);

  for (i=0; i<ygen->lmax+2; ++i)
    gen->fx[i][0]=gen->fx[i][1]=gen->fx[i][2]=0.;

  gen->cth_crit = 2.;
  gen->mdist_crit = ygen->lmax+1;
  }

static void sylmgen_destroy (sylmgen_d *gen)
  { DEALLOC(gen->fx); }

static void sylmgen_prepare (sylmgen_d *gen, const Ylmgen_C *ygen, int m_)
  {
  int mlo_, mhi_, ms_similar, l;

  if (m_==gen->m) return;
  UTIL_ASSERT(m_>=0,"incorrect m");

  mlo_=m_; mhi_=gen->s;
  if (mhi_<mlo_) SWAP(mhi_,mlo_,int);
  ms_similar = ((gen->mhi==mhi_) && (gen->mlo==mlo_));

  gen->m=m_;
  gen->mlo = mlo_; gen->mhi = mhi_;

  if (!ms_similar)
    {
    for (l=gen->mhi; l<ygen->lmax; ++l)
      {
      double t = ygen->flm1[l+gen->m]*ygen->flm1[l-gen->m]
                *ygen->flm1[l+gen->s]*ygen->flm1[l-gen->s];
      double lt = 2*l+1;
      double l1 = l+1;
      gen->fx[l+1][0]=l1*lt*t;
      gen->fx[l+1][1]=gen->m*gen->s*ygen->xl[l]*ygen->xl[l+1];
      t = ygen->flm2[l+gen->m]*ygen->flm2[l-gen->m]
         *ygen->flm2[l+gen->s]*ygen->flm2[l-gen->s];
      gen->fx[l+1][2]=t*l1*ygen->xl[l];
      }
    gen->prefactor = 0.5L*(ygen->logsum[2*gen->mhi]
      -ygen->logsum[gen->mhi+gen->mlo]-ygen->logsum[gen->mhi-gen->mlo]);
    }

  gen->preMinus_p = gen->preMinus_m = 0;
  if (gen->mhi==gen->m)
    {
    gen->cosPow = gen->mhi+gen->s; gen->sinPow = gen->mhi-gen->s;
    gen->preMinus_p = gen->preMinus_m = ((gen->mhi-gen->s)&1);
    }
  else
    {
    gen->cosPow = gen->mhi+gen->m; gen->sinPow = gen->mhi-gen->m;
    gen->preMinus_m = ((gen->mhi+gen->m)&1);
    }
  }

static void sylmgen_recalc (sylmgen_d *gen, const Ylmgen_C *ygen, int ith,
  ylmgen_dbl2 *res, int *firstl)
  {
  const double ln2     = 0.6931471805599453094172321214581766;
  const double inv_ln2 = 1.4426950408889634073599246810018921;
  int l=gen->mhi;
  int lmax = ygen->lmax;
  ylmgen_dbl3 *fy = gen->fx;
  const double fsmall = ygen->fsmall, fbig = ygen->fbig, eps = ygen->eps;
  const double cth = ygen->cth[ith];
  long double logvalp = inv_ln2*(gen->prefactor
    + ygen->lc05[ith]*gen->cosPow + ygen->ls05[ith]*gen->sinPow);
  long double logvalm = inv_ln2*(gen->prefactor
    + ygen->lc05[ith]*gen->sinPow + ygen->ls05[ith]*gen->cosPow);
  int scalep = (int)(logvalp/large_exponent2)-minscale;
  int scalem = (int)(logvalm/large_exponent2)-minscale;
  double rec1p=0., rec1m=0.;
  double rec2p = exp(ln2*(double)(logvalp-(scalep+minscale)*large_exponent2));
  double rec2m = exp(ln2*(double)(logvalm-(scalem+minscale)*large_exponent2));
  double corfacp,corfacm;
  double tp,tm;

  if ((abs(gen->m-gen->s)>=gen->mdist_crit)&&(fabs(cth)>=gen->cth_crit))
    { *firstl=ygen->lmax+1; return; }

  if (gen->preMinus_p)
    rec2p=-rec2p;
  if (gen->preMinus_m)
    rec2m=-rec2m;
  if (gen->s&1)
    rec2p=-rec2p;

  /* iterate until we reach the realm of IEEE numbers */
  while((scalem<0)&&(scalep<0))
    {
    if (++l>lmax) break;
    rec1p = (cth - fy[l][1])*fy[l][0]*rec2p - fy[l][2]*rec1p;
    rec1m = (cth + fy[l][1])*fy[l][0]*rec2m - fy[l][2]*rec1m;
    if (++l>lmax) break;
    rec2p = (cth - fy[l][1])*fy[l][0]*rec1p - fy[l][2]*rec2p;
    rec2m = (cth + fy[l][1])*fy[l][0]*rec1m - fy[l][2]*rec2m;

    while (fabs(rec2p)>fbig)
      { rec1p *= fsmall; rec2p *= fsmall; ++scalep; }
    while (fabs(rec2m)>fbig)
      { rec1m *= fsmall; rec2m *= fsmall; ++scalem; }
    }

  corfacp = (scalep<0) ? 0. : ygen->cf[scalep];
  corfacm = (scalem<0) ? 0. : ygen->cf[scalem];

  if (l<=lmax)
    {
    while (1)
      {
      if ((fabs(rec2p*corfacp)>eps) || (fabs(rec2m*corfacm)>eps))
        break;
      if (++l>lmax) break;
      rec1p = (cth - fy[l][1])*fy[l][0]*rec2p - fy[l][2]*rec1p;
      rec1m = (cth + fy[l][1])*fy[l][0]*rec2m - fy[l][2]*rec1m;
      if ((fabs(rec1p*corfacp)>eps) || (fabs(rec1m*corfacm)>eps))
        { SWAP(rec1p,rec2p,double); SWAP(rec1m,rec2m,double); break; }
      if (++l>lmax) break;
      rec2p = (cth - fy[l][1])*fy[l][0]*rec1p - fy[l][2]*rec2p;
      rec2m = (cth + fy[l][1])*fy[l][0]*rec1m - fy[l][2]*rec2m;

      if ((fabs(rec2p)>fbig) || (fabs(rec2m)>fbig))
        {
        while (fabs(rec2p)>fbig)
          { rec1p *= fsmall; rec2p *= fsmall; ++scalep; }
        while (fabs(rec2m)>fbig)
          { rec1m *= fsmall; rec2m *= fsmall; ++scalem; }
        corfacp = (scalep<0) ? 0. : ygen->cf[scalep];
        corfacm = (scalem<0) ? 0. : ygen->cf[scalem];
        }
      }
    }

  *firstl=l;
  if (l>lmax)
    {
    gen->mdist_crit=abs(gen->m-gen->s);
    gen->cth_crit= fabs(cth);
    return;
    }

  tp = rec2p*corfacp; tm = rec2m*corfacm;
  res[l][0]=tp+tm;
  res[l][1]=tm-tp;

  while (1)
    {
    if ((fabs(tp)>eps) && (fabs(tm)>eps))
      break;
    if (++l>lmax) break;
    rec1p = (cth - fy[l][1])*fy[l][0]*rec2p - fy[l][2]*rec1p;
    rec1m = (cth + fy[l][1])*fy[l][0]*rec2m - fy[l][2]*rec1m;
    tp=rec1p*corfacp; tm=rec1m*corfacm;
    res[l][0]=tp+tm; res[l][1]=tm-tp;
    if ((fabs(tp)>eps) && (fabs(tm)>eps))
      { SWAP(rec1p,rec2p,double); SWAP(rec1m,rec2m,double); break; }
    if (++l>lmax) break;
    rec2p = (cth - fy[l][1])*fy[l][0]*rec1p - fy[l][2]*rec2p;
    rec2m = (cth + fy[l][1])*fy[l][0]*rec1m - fy[l][2]*rec2m;
    tp=rec2p*corfacp; tm=rec2m*corfacm;
    res[l][0]=tp+tm; res[l][1]=tm-tp;

    if ((fabs(rec2p)>fbig) || (fabs(rec2m)>fbig))
      {
      while (fabs(rec2p)>fbig)
        { rec1p *= fsmall; rec2p *= fsmall; ++scalep; }
      while (fabs(rec2m)>fbig)
        { rec1m *= fsmall; rec2m *= fsmall; ++scalem; }
      corfacp = (scalep<0) ? 0. : ygen->cf[scalep];
      corfacm = (scalem<0) ? 0. : ygen->cf[scalem];
      }
    }

  rec1p *= corfacp; rec2p *= corfacp;
  rec1m *= corfacm; rec2m *= corfacm;

  for (;l<lmax-1;l+=2)
    {
    rec1p = (cth - fy[l+1][1])*fy[l+1][0]*rec2p - fy[l+1][2]*rec1p;
    rec1m = (cth + fy[l+1][1])*fy[l+1][0]*rec2m - fy[l+1][2]*rec1m;
    res[l+1][0] = rec1p+rec1m; res[l+1][1] = rec1m-rec1p;
    rec2p = (cth - fy[l+2][1])*fy[l+2][0]*rec1p - fy[l+2][2]*rec2p;
    rec2m = (cth + fy[l+2][1])*fy[l+2][0]*rec1m - fy[l+2][2]*rec2m;
    res[l+2][0] = rec2p+rec2m; res[l+2][1] = rec2m-rec2p;
    }
  while (1)
    {
    if (++l>lmax) break;
    rec1p = (cth - fy[l][1])*fy[l][0]*rec2p - fy[l][2]*rec1p;
    rec1m = (cth + fy[l][1])*fy[l][0]*rec2m - fy[l][2]*rec1m;
    res[l][0] = rec1p+rec1m; res[l][1] = rec1m-rec1p;
    if (++l>lmax) break;
    rec2p = (cth - fy[l][1])*fy[l][0]*rec1p - fy[l][2]*rec2p;
    rec2m = (cth + fy[l][1])*fy[l][0]*rec1m - fy[l][2]*rec2m;
    res[l][0] = rec2p+rec2m; res[l][1] = rec2m-rec2p;
    }
  }

#ifdef PLANCK_HAVE_SSE2

#define ADVANCE(L,ap,am,bp,bm) \
  { \
  v2df f0=_mm_set1_pd(fy[L][0]), \
       f1=_mm_set1_pd(fy[L][1]), \
       f2=_mm_set1_pd(fy[L][2]); \
  ap = _mm_sub_pd( \
       _mm_mul_pd(_mm_sub_pd(cth,f1),_mm_mul_pd(f0,bp)), \
       _mm_mul_pd(f2,ap)); \
  am = _mm_sub_pd( \
       _mm_mul_pd(_mm_add_pd(cth,f1),_mm_mul_pd(f0,bm)), \
       _mm_mul_pd(f2,am)); \
  }

#define RENORMX(r1,r2,corf,sca,scb) \
  do \
    { \
    double rec1a, rec1b, rec2a, rec2b, corfaca, corfacb; \
    read_v2df (r1, &rec1a, &rec1b); read_v2df (r2, &rec2a, &rec2b); \
    read_v2df (corf, &corfaca, &corfacb); \
    while (fabs(rec2a)>fbig) \
      { \
      rec1a*=fsmall; rec2a*=fsmall; ++sca; \
      corfaca = (sca<0) ? 0. : ygen->cf[sca]; \
      } \
    while (fabs(rec2b)>fbig) \
      { \
      rec1b*=fsmall; rec2b*=fsmall; ++scb; \
      corfacb = (scb<0) ? 0. : ygen->cf[scb]; \
      } \
    r1=build_v2df(rec1a,rec1b); r2=build_v2df(rec2a,rec2b); \
    corf=build_v2df(corfaca,corfacb); \
    } \
  while(0)

#define RENORM \
  RENORMX(rec1p,rec2p,corfacp,scale1p,scale2p); \
  RENORMX(rec1m,rec2m,corfacm,scale1m,scale2m)

static void sylmgen_recalc_sse2 (sylmgen_d *gen, const Ylmgen_C *ygen,
  int ith1, int ith2, v2df2 *res, int *firstl)
  {
  const double ln2     = 0.6931471805599453094172321214581766;
  const double inv_ln2 = 1.4426950408889634073599246810018921;
  int l=gen->mhi;
  int lmax = ygen->lmax;
  ylmgen_dbl3 *fy = gen->fx;
  const double fbig=ygen->fbig, fsmall=ygen->fsmall;
  v2df eps2    = build_v2df(ygen->eps,ygen->eps);
  const double cth1=ygen->cth[ith1], cth2=ygen->cth[ith2];
  v2df cth = build_v2df(cth1,cth2);
  long double
    logval1p = inv_ln2*(gen->prefactor
             + ygen->lc05[ith1]*gen->cosPow + ygen->ls05[ith1]*gen->sinPow),
    logval2p = inv_ln2*(gen->prefactor
             + ygen->lc05[ith2]*gen->cosPow + ygen->ls05[ith2]*gen->sinPow),
    logval1m = inv_ln2*(gen->prefactor
             + ygen->lc05[ith1]*gen->sinPow + ygen->ls05[ith1]*gen->cosPow),
    logval2m = inv_ln2*(gen->prefactor
             + ygen->lc05[ith2]*gen->sinPow + ygen->ls05[ith2]*gen->cosPow);

  int scale1p = (int)(logval1p/large_exponent2)-minscale,
      scale2p = (int)(logval2p/large_exponent2)-minscale,
      scale1m = (int)(logval1m/large_exponent2)-minscale,
      scale2m = (int)(logval2m/large_exponent2)-minscale;

  v2df rec1p =_mm_setzero_pd(), rec1m=_mm_setzero_pd();
  v2df rec2p = build_v2df(
    exp(ln2*(double)(logval1p-(scale1p+minscale)*large_exponent2)),
    exp(ln2*(double)(logval2p-(scale2p+minscale)*large_exponent2)));
  v2df rec2m = build_v2df(
    exp(ln2*(double)(logval1m-(scale1m+minscale)*large_exponent2)),
    exp(ln2*(double)(logval2m-(scale2m+minscale)*large_exponent2)));
  v2df corfacp=build_v2df((scale1p<0) ? 0. : ygen->cf[scale1p],
                          (scale2p<0) ? 0. : ygen->cf[scale2p]),
       corfacm=build_v2df((scale1m<0) ? 0. : ygen->cf[scale1m],
                          (scale2m<0) ? 0. : ygen->cf[scale2m]);
  v2df tp,tm;

  if ((abs(gen->m-gen->s)>=gen->mdist_crit)
     &&(fabs(ygen->cth[ith1])>=gen->cth_crit)
     &&(fabs(ygen->cth[ith2])>=gen->cth_crit))
    { *firstl=ygen->lmax+1; return; }

  if (gen->preMinus_p)
    rec2p = _mm_xor_pd (rec2p,V2DF_SIGNMASK); /* negate */
  if (gen->preMinus_m)
    rec2m = _mm_xor_pd (rec2m,V2DF_SIGNMASK); /* negate */
  if (gen->s&1)
    rec2p = _mm_xor_pd (rec2p,V2DF_SIGNMASK); /* negate */

  /* iterate until we reach the realm of IEEE numbers */
  while((scale1m<0)&&(scale1p<0)&&(scale2m<0)&&(scale2p<0))
    {
    if (++l>lmax) break;
    ADVANCE (l,rec1p,rec1m,rec2p,rec2m)
    if (++l>lmax) break;
    ADVANCE (l,rec2p,rec2m,rec1p,rec1m)

    RENORM;
    }

  if (l<=lmax)
    {
    while (1)
      {
      if (v2df_any_gt(_mm_mul_pd(rec2p,corfacp),eps2) ||
          v2df_any_gt(_mm_mul_pd(rec2m,corfacm),eps2))
        break;
      if (++l>lmax) break;
      ADVANCE (l,rec1p,rec1m,rec2p,rec2m)
      if (v2df_any_gt(_mm_mul_pd(rec1p,corfacp),eps2) ||
          v2df_any_gt(_mm_mul_pd(rec1m,corfacm),eps2))
        { SWAP(rec1p,rec2p,v2df); SWAP(rec1m,rec2m,v2df); break; }
      if (++l>lmax) break;
      ADVANCE (l,rec2p,rec2m,rec1p,rec1m)

      RENORM;
      }
    }

  *firstl=l;
  if (l>lmax)
    {
    gen->mdist_crit=abs(gen->m-gen->s);
    gen->cth_crit= (fabs(cth1)<fabs(cth2)) ? fabs(cth1) : fabs(cth2);
    return;
    }
  tp = _mm_mul_pd(rec2p,corfacp); tm = _mm_mul_pd(rec2m,corfacm);
  res[l].a=_mm_add_pd(tp,tm);
  res[l].b=_mm_sub_pd(tm,tp);

  while (1)
    {
    if (v2df_all_ge(tp,eps2) && v2df_all_ge(tm,eps2))
      break;
    if (++l>lmax) break;
    ADVANCE(l,rec1p,rec1m,rec2p,rec2m)
    tp=_mm_mul_pd(rec1p,corfacp); tm=_mm_mul_pd(rec1m,corfacm);
    res[l].a=_mm_add_pd(tp,tm); res[l].b=_mm_sub_pd(tm,tp);
    if (v2df_all_ge(tp,eps2) && v2df_all_ge(tm,eps2))
      { SWAP(rec1p,rec2p,v2df); SWAP(rec1m,rec2m,v2df); break; }
    if (++l>lmax) break;
    ADVANCE (l,rec2p,rec2m,rec1p,rec1m)
    tp=_mm_mul_pd(rec2p,corfacp); tm=_mm_mul_pd(rec2m,corfacm);
    res[l].a=_mm_add_pd(tp,tm); res[l].b=_mm_sub_pd(tm,tp);

    RENORM;
    }

  rec1p = _mm_mul_pd(rec1p,corfacp); rec2p = _mm_mul_pd(rec2p,corfacp);
  rec1m = _mm_mul_pd(rec1m,corfacm); rec2m = _mm_mul_pd(rec2m,corfacm);

  for (;l<lmax-1;l+=2)
    {
    ADVANCE(l+1,rec1p,rec1m,rec2p,rec2m)
    res[l+1].a=_mm_add_pd(rec1p,rec1m); res[l+1].b=_mm_sub_pd(rec1m,rec1p);
    ADVANCE(l+2,rec2p,rec2m,rec1p,rec1m)
    res[l+2].a=_mm_add_pd(rec2p,rec2m); res[l+2].b=_mm_sub_pd(rec2m,rec2p);
    }
  while (1)
    {
    if (++l>lmax) break;
    ADVANCE(l,rec1p,rec1m,rec2p,rec2m)
    res[l].a=_mm_add_pd(rec1p,rec1m); res[l].b=_mm_sub_pd(rec1m,rec1p);
    if (++l>lmax) break;
    ADVANCE(l,rec2p,rec2m,rec1p,rec1m)
    res[l].a=_mm_add_pd(rec2p,rec2m); res[l].b=_mm_sub_pd(rec2m,rec2p);
    }
  }

#endif

void Ylmgen_init (Ylmgen_C *gen, int l_max, int m_max, int spinrec,
  double epsilon)
  {
  int m;
  const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;
  const double inv_ln2 = 1.4426950408889634073599246810018921;

  gen->fsmall = ldexp(1.,-large_exponent2);
  gen->fbig   = ldexp(1., large_exponent2);
  gen->eps = epsilon;
  gen->cth_crit = 2.;
  gen->ith = -1;
  gen->nth = 0;
  gen->lmax = l_max;
  gen->mmax = m_max;
  gen->spinrec = spinrec;
  gen->m_cur = -1;
  gen->m_crit = gen->mmax+1;
  gen->firstl = RALLOC(int,max_spin+1);
  for (m=0; m<=max_spin; ++m) gen->firstl[m]=-1;
  gen->cf = RALLOC(double,maxscale-minscale+1);
  for (m=0; m<(maxscale-minscale+1); ++m)
    gen->cf[m] = ldexp(1.,(m+minscale)*large_exponent2);
  gen->recfac = RALLOC(ylmgen_dbl2,gen->lmax+1);
  gen->mfac = RALLOC(double,gen->mmax+1);
  gen->mfac[0] = 1;
  for (m=1; m<=gen->mmax; ++m)
    gen->mfac[m] = gen->mfac[m-1]*sqrt((2*m+1.)/(2*m));
  for (m=0; m<=gen->mmax; ++m)
    gen->mfac[m] = inv_ln2*log(inv_sqrt4pi*gen->mfac[m]);

  gen->t1fac = RALLOC(double,gen->lmax+1);
  for (m=0; m<=gen->lmax; ++m)
    gen->t1fac[m] = sqrt(4.*(m+1)*(m+1)-1.);
  gen->t2fac = RALLOC(double,gen->lmax+gen->mmax+1);
  for (m=0; m<=gen->lmax+gen->mmax; ++m)
    gen->t2fac[m] = 1./sqrt(m+1.);

  gen->lamfact = RALLOC(double,gen->lmax+1);
  gen->ylm = RALLOC(double,gen->lmax+1);
  ALLOC(gen->lambda_wx,ylmgen_dbl2 *,max_spin+1);
  for (m=0; m<=max_spin; ++m)
    gen->lambda_wx[m]=NULL;

  gen->sylm = RALLOC(sylmgen_d *,max_spin+1);
  for (m=0; m<=max_spin; ++m)
    gen->sylm[m]=NULL;

  gen->ylm_uptodate = 0;
  gen->lwx_uptodate = RALLOC(int,max_spin+1);
  SET_ARRAY(gen->lwx_uptodate,0,max_spin+1,0);
  gen->recfac_uptodate = 0;
  gen->lamfact_uptodate = 0;

  gen->th = gen->cth = gen->sth = gen->logsth = NULL;

#ifdef PLANCK_HAVE_SSE2
  gen->ith1 = gen->ith2 = -1;
  gen->ylm_sse2 = RALLOC(v2df,gen->lmax+1);
  ALLOC(gen->lambda_wx_sse2,v2df2 *,max_spin+1);
  for (m=0; m<=max_spin; ++m)
    gen->lambda_wx_sse2[m]=NULL;
  gen->ylm_uptodate_sse2 = 0;
  gen->lwx_uptodate_sse2 = RALLOC(int,max_spin+1);
  SET_ARRAY(gen->lwx_uptodate_sse2,0,max_spin+1,0);
#endif

  ALLOC(gen->logsum,long double,2*gen->lmax+1);
  gen->lc05 = gen->ls05 = NULL;
  ALLOC(gen->flm1,double,2*gen->lmax+1);
  ALLOC(gen->flm2,double,2*gen->lmax+1);
  ALLOC(gen->xl,double,gen->lmax+1);

  gen->logsum[0] = 0.;
  for (m=1; m<2*gen->lmax+1; ++m)
    gen->logsum[m] = gen->logsum[m-1]+logl((long double)m);
  for (m=0; m<2*gen->lmax+1; ++m)
    {
    gen->flm1[m] = sqrt(1./(m+1.));
    gen->flm2[m] = sqrt(m/(m+1.));
    }

  gen->xl[0]=0;
  for (m=1; m<gen->lmax+1; ++m) gen->xl[m]=1./m;
  }

void Ylmgen_destroy (Ylmgen_C *gen)
  {
  int m;

  DEALLOC(gen->firstl);
  DEALLOC(gen->cf);
  DEALLOC(gen->recfac);
  DEALLOC(gen->mfac);
  DEALLOC(gen->t1fac);
  DEALLOC(gen->t2fac);
  DEALLOC(gen->lamfact);
  DEALLOC(gen->ylm);
  DEALLOC(gen->lwx_uptodate);
  for (m=0; m<=max_spin; ++m)
    DEALLOC(gen->lambda_wx[m]);
  DEALLOC(gen->lambda_wx);
  for (m=0; m<=max_spin; ++m)
    if (gen->sylm[m])
      {
      sylmgen_destroy (gen->sylm[m]);
      DEALLOC(gen->sylm[m]);
      }
  DEALLOC(gen->sylm);
  DEALLOC(gen->th);
  DEALLOC(gen->cth);
  DEALLOC(gen->sth);
  DEALLOC(gen->logsth);
  DEALLOC(gen->logsum);
  DEALLOC(gen->lc05);
  DEALLOC(gen->ls05);
  DEALLOC(gen->flm1);
  DEALLOC(gen->flm2);
  DEALLOC(gen->xl);
#ifdef PLANCK_HAVE_SSE2
  DEALLOC(gen->ylm_sse2);
  for (m=0; m<=max_spin; ++m)
    DEALLOC(gen->lambda_wx_sse2[m]);
  DEALLOC(gen->lambda_wx_sse2);
  DEALLOC(gen->lwx_uptodate_sse2);
#endif
  }

void Ylmgen_set_theta (Ylmgen_C *gen, const double *theta, int nth)
  {
  const double inv_ln2 = 1.4426950408889634073599246810018921;
  int m;
  DEALLOC(gen->th);
  DEALLOC(gen->cth);
  DEALLOC(gen->sth);
  DEALLOC(gen->logsth);
  DEALLOC(gen->lc05);
  DEALLOC(gen->ls05);
  gen->th = RALLOC(double,nth);
  gen->cth = RALLOC(double,nth);
  gen->sth = RALLOC(double,nth);
  gen->logsth = RALLOC(double,nth);
  gen->lc05 = RALLOC(long double,nth);
  gen->ls05 = RALLOC(long double,nth);
  for (m=0; m<nth; ++m)
    {
    const double pi = 3.141592653589793238462643383279502884197;
    double th=theta[m];
    UTIL_ASSERT ((th>=0.)&&(th<=pi),"bad theta angle specified");
    /* tiny adjustments to make sure cos and sin (theta/2) are positive */
    if (th==0.) th=1e-16;
    if (ABSAPPROX(th,pi,1e-15)) th=pi-1e-15;
    gen->th[m] = th;
    gen->cth[m] = cos(th);
    gen->sth[m] = sin(th);
    gen->logsth[m] = inv_ln2*log(gen->sth[m]);
    gen->lc05[m]=logl(cosl(0.5L*th));
    gen->ls05[m]=logl(sinl(0.5L*th));
    }

  gen->nth = nth;
  gen->ith = -1;
#ifdef PLANCK_HAVE_SSE2
  gen->ith1 = gen->ith2 = -1;
#endif
  }

void Ylmgen_prepare (Ylmgen_C *gen, int ith, int m)
  {
  if ((ith==gen->ith) && (m==gen->m_cur)) return;

  gen->ylm_uptodate = 0;
  SET_ARRAY(gen->lwx_uptodate,0,max_spin+1,0);

  gen->ith = ith;

  if (m!=gen->m_cur)
    {
    gen->recfac_uptodate = 0;
    gen->lamfact_uptodate = 0;
    gen->m_cur = m;
    }
  }

static void Ylmgen_recalc_recfac (Ylmgen_C *gen)
  {
  double f_old=1;
  int l, m;

  if (gen->recfac_uptodate) return;
  gen->recfac_uptodate = 1;

  m = gen->m_cur;
  for (l=m; l<=gen->lmax; ++l)
    {
    gen->recfac[l][0] = gen->t1fac[l]*gen->t2fac[l+m]*gen->t2fac[l-m];
    gen->recfac[l][1] = gen->recfac[l][0]/f_old;
    f_old = gen->recfac[l][0];
    }
  }

static void Ylmgen_recalc_lamfact (Ylmgen_C *gen)
  {
  int l, m;

  if (gen->lamfact_uptodate) return;
  gen->lamfact_uptodate = 1;

  m = gen->m_cur;
  gen->lamfact[m] = 0;
  for (l=m+1; l<=gen->lmax; ++l)
    gen->lamfact[l] = sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
  }

#define RENORMALIZE_SCALAR \
  do \
    { \
    while (fabs(lam_2)>fbig) \
      { lam_1*=fsmall; lam_2*=fsmall; ++scale; } \
    corfac = (scale<0) ? 0. : gen->cf[scale]; \
    } \
  while(0)

void Ylmgen_recalc_Ylm (Ylmgen_C *gen)
  {
  const double ln2 = 0.6931471805599453094172321214581766;

  double logval,lam_1,lam_2,corfac;
  double eps=gen->eps, fbig=gen->fbig, fsmall=gen->fsmall;
  ylmgen_dbl2 *recfac = gen->recfac;
  int lmax=gen->lmax;
  int scale,l;
  int m = gen->m_cur;
  double cth=gen->cth[gen->ith], sth=gen->sth[gen->ith];
  double *result = gen->ylm;

  if (gen->ylm_uptodate) return;
  gen->ylm_uptodate=1;

  if (((m>=gen->m_crit)&&(fabs(cth)>=gen->cth_crit)) || ((m>0)&&(sth==0)))
    { gen->firstl[0]=gen->lmax+1; return; }

  Ylmgen_recalc_recfac(gen);

  logval = gen->mfac[m];
  if (m>0) logval += m*gen->logsth[gen->ith];
  scale = (int) (logval/large_exponent2)-minscale;
  corfac = (scale<0) ? 0. : gen->cf[scale];

  lam_1 = 0;
  lam_2 = exp(ln2*(logval-(scale+minscale)*large_exponent2));
  if (m&1) lam_2 = -lam_2;

  l=m;
  if (scale<0)
    {
    while (1)
      {
      if (++l>lmax) break;
      lam_1 = cth*lam_2*recfac[l-1][0] - lam_1*recfac[l-1][1];
      if (++l>lmax) break;
      lam_2 = cth*lam_1*recfac[l-1][0] - lam_2*recfac[l-1][1];
      if (fabs(lam_2)>fbig)
        {
        RENORMALIZE_SCALAR;
        if (scale>=0) break;
        }
      }
    }

  lam_1*=corfac;
  lam_2*=corfac;

  if (l<=lmax)
    {
    while (1)
      {
      if (fabs(lam_2)>eps) break;
      if (++l>lmax) break;
      lam_1 = cth*lam_2*recfac[l-1][0] - lam_1*recfac[l-1][1];
      if (fabs(lam_1)>eps)
        { double x=lam_1; lam_1=lam_2; lam_2=x; break; }
      if (++l>lmax) break;
      lam_2 = cth*lam_1*recfac[l-1][0] - lam_2*recfac[l-1][1];
      }
    }

  gen->firstl[0]=l;
  if (l>lmax)
    { gen->m_crit=m; gen->cth_crit=fabs(cth); return; }

  for(;l<lmax-3;l+=4)
    {
    result[l]=lam_2;
    lam_1 = cth*lam_2*recfac[l][0] - lam_1*recfac[l][1];
    result[l+1] = lam_1;
    lam_2 = cth*lam_1*recfac[l+1][0] - lam_2*recfac[l+1][1];
    result[l+2]=lam_2;
    lam_1 = cth*lam_2*recfac[l+2][0] - lam_1*recfac[l+2][1];
    result[l+3] = lam_1;
    lam_2 = cth*lam_1*recfac[l+3][0] - lam_2*recfac[l+3][1];
    }

  while (1)
    {
    result[l]=lam_2;
    if (++l>lmax) break;
    lam_1 = cth*lam_2*recfac[l-1][0] - lam_1*recfac[l-1][1];
    result[l] = lam_1;
    if (++l>lmax) break;
    lam_2 = cth*lam_1*recfac[l-1][0] - lam_2*recfac[l-1][1];
    }
  }

static void Ylmgen_recalc_lambda_wx1 (Ylmgen_C *gen)
  {
  if (gen->lwx_uptodate[1]) return;
  Ylmgen_recalc_Ylm(gen);
  gen->firstl[1] = gen->firstl[0];
  if (gen->firstl[1]>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate[1] = 1;

  {
  double cth=gen->cth[gen->ith];
  double xsth=1./(gen->sth[gen->ith]);
  double m=gen->m_cur;
  double m_on_sth = m*xsth;
  double lam_lm=0;
  ylmgen_dbl2 *lambda_wx = gen->lambda_wx[1];
  int l;
  double ell;
  for (ell=l=gen->firstl[1]; l<=gen->lmax; ++l, ell+=1.)
    {
    double lam_lm1m=lam_lm;
    lam_lm=gen->ylm[l];
    lambda_wx[l][0] = xsth*(gen->lamfact[l]*lam_lm1m - ell*cth*lam_lm);
    lambda_wx[l][1] = m_on_sth*lam_lm;
    }
  }
  }

static void Ylmgen_recalc_lambda_wx2 (Ylmgen_C *gen)
  {
  if (gen->lwx_uptodate[2]) return;
  Ylmgen_recalc_Ylm(gen);
  gen->firstl[2] = gen->firstl[0];
  if (gen->firstl[2]>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate[2] = 1;

  {
  double cth=gen->cth[gen->ith];
  double sth=gen->sth[gen->ith];
  double m=gen->m_cur;
  double one_on_s2 = 1./(sth*sth);
  double two_on_s2 = 2*one_on_s2;
  double two_c_on_s2 = cth * two_on_s2;
  double m2 = m*m;
  double two_m_on_s2 = m*two_on_s2;
  double lam_lm=0;
  ylmgen_dbl2 *lambda_wx = gen->lambda_wx[2];
  int l;
  double ell;
  for (ell=l=gen->firstl[2]; l<=gen->lmax; ++l, ell+=1.)
    {
    double lam_lm1m=lam_lm;
    lam_lm=gen->ylm[l];
    {
    const double t1  = lam_lm1m*gen->lamfact[l];
    const double a_w = (m2-ell)*two_on_s2 - ell*(ell-1.);
    const double a_x = cth*(ell-1.)*lam_lm;
    lambda_wx[l][0] = a_w*lam_lm + t1*two_c_on_s2;
    lambda_wx[l][1] = two_m_on_s2 * (t1-a_x);
    }
    }
  }
  }

void Ylmgen_recalc_lambda_wx (Ylmgen_C *gen, int spin)
  {
  UTIL_ASSERT ((spin>0) && (spin<=max_spin),
    "invalid spin in Ylmgen_recalc_lambda_wx");

  if (!gen->lambda_wx[spin])
    gen->lambda_wx[spin]=RALLOC(ylmgen_dbl2,gen->lmax+1);

  if (gen->spinrec && spin==1) { Ylmgen_recalc_lambda_wx1(gen); return; }
  if (gen->spinrec && spin==2) { Ylmgen_recalc_lambda_wx2(gen); return; }

  if (!gen->sylm[spin])
    {
    gen->sylm[spin]=RALLOC(sylmgen_d,1);
    sylmgen_init(gen->sylm[spin],gen,spin);
    }
  if (gen->lwx_uptodate[spin]) return;
  sylmgen_prepare(gen->sylm[spin],gen,gen->m_cur);
  sylmgen_recalc(gen->sylm[spin],gen,gen->ith,gen->lambda_wx[spin],
    &gen->firstl[spin]);
  gen->lwx_uptodate[spin] = 1;
  }

#ifdef PLANCK_HAVE_SSE2

void Ylmgen_prepare_sse2 (Ylmgen_C *gen, int ith1, int ith2, int m)
  {
  if ((ith1==gen->ith1) && (ith2==gen->ith2) && (m==gen->m_cur)) return;

  gen->ylm_uptodate_sse2 = 0;
  SET_ARRAY(gen->lwx_uptodate_sse2,0,max_spin+1,0);

  gen->ith1 = ith1; gen->ith2 = ith2;

  if (m!=gen->m_cur)
    {
    gen->recfac_uptodate = gen->lamfact_uptodate = 0;
    gen->m_cur = m;
    }
  }


#define RENORMALIZE \
  do \
    { \
    double lam1a, lam1b, lam2a, lam2b, corfaca, corfacb; \
    read_v2df (lam_1, &lam1a, &lam1b); read_v2df (lam_2, &lam2a, &lam2b); \
    read_v2df (corfac, &corfaca, &corfacb); \
    while (fabs(lam2a)>fbig) \
      { \
      lam1a*=fsmall; lam2a*=fsmall; ++scale1; \
      corfaca = (scale1<0) ? 0. : gen->cf[scale1]; \
      } \
    while (fabs(lam2b)>fbig) \
      { \
      lam1b*=fsmall; lam2b*=fsmall; ++scale2; \
      corfacb = (scale2<0) ? 0. : gen->cf[scale2]; \
      } \
    lam_1=build_v2df(lam1a,lam1b); lam_2=build_v2df(lam2a,lam2b); \
    corfac=build_v2df(corfaca,corfacb); \
    } \
  while(0)
#define GETPRE(prea,preb,lv) \
  { \
  prea=_mm_mul_pd(_mm_set1_pd(recfac[lv][0]),cth); \
  preb=_mm_set1_pd(recfac[lv][1]); \
  }
#define NEXTSTEP(prea,preb,prec,pred,reca,recb,lv) \
  { \
  preb = _mm_mul_pd(preb,reca); \
  prea = _mm_mul_pd(prea,recb); \
  prec = _mm_set1_pd(recfac[lv][0]); \
  pred = _mm_set1_pd(recfac[lv][1]); \
  reca = _mm_sub_pd(prea,preb); \
  prec = _mm_mul_pd(cth,prec); \
  }

void Ylmgen_recalc_Ylm_sse2 (Ylmgen_C *gen)
  {
  const double ln2 = 0.6931471805599453094172321214581766;

  v2df lam_1,lam_2,corfac;
  double logval1,logval2;
  double eps=gen->eps, fbig=gen->fbig, fsmall=gen->fsmall;
  v2df eps2=build_v2df(eps,eps);
  v2df fbig2=build_v2df(fbig,fbig);
  ylmgen_dbl2 *recfac = gen->recfac;
  int lmax=gen->lmax;
  int scale1,scale2,l;
  int m = gen->m_cur;
  double cth1=gen->cth[gen->ith1], cth2=gen->cth[gen->ith2];
  v2df cth=build_v2df(cth1,cth2);
  v2df *result = gen->ylm_sse2;
  v2df pre0,pre1,pre2,pre3;

  if (gen->ylm_uptodate_sse2) return;
  gen->ylm_uptodate_sse2=1;

  if ((m>=gen->m_crit)&&(fabs(cth1)>=gen->cth_crit)&&(fabs(cth2)>=gen->cth_crit))
    { gen->firstl[0]=gen->lmax+1; return; }

  Ylmgen_recalc_recfac(gen);

  logval1 = logval2 = gen->mfac[m];
  if (m>0) logval1 += m*gen->logsth[gen->ith1];
  if (m>0) logval2 += m*gen->logsth[gen->ith2];
  scale1 = (int) (logval1/large_exponent2)-minscale;
  scale2 = (int) (logval2/large_exponent2)-minscale;
  corfac = build_v2df((scale1<0) ? 0. : gen->cf[scale1],
                      (scale2<0) ? 0. : gen->cf[scale2]);

  lam_1 = _mm_setzero_pd();
  lam_2 = build_v2df(exp(ln2*(logval1-(scale1+minscale)*large_exponent2)),
                     exp(ln2*(logval2-(scale2+minscale)*large_exponent2)));
  if (m&1) lam_2 = _mm_xor_pd (lam_2,V2DF_SIGNMASK); /* negate */

  l=m;
  if ((scale1<0) && (scale2<0))
    {
    GETPRE(pre0,pre1,l)
    while (1)
      {
      if (++l>lmax) break;
      NEXTSTEP(pre0,pre1,pre2,pre3,lam_1,lam_2,l)
      if (++l>lmax) break;
      NEXTSTEP(pre2,pre3,pre0,pre1,lam_2,lam_1,l)
      if (v2df_any_gt(lam_2,fbig2))
        {
        RENORMALIZE;
        if ((scale1>=0) || (scale2>=0)) break;
        }
      }
    }

  if (l<=lmax)
    {
    GETPRE(pre0,pre1,l)
    while (1)
      {
      v2df t1;
      result[l]=t1=_mm_mul_pd(lam_2,corfac);
      if (v2df_any_gt(t1,eps2))
        break;
      if (++l>lmax) break;
      NEXTSTEP(pre0,pre1,pre2,pre3,lam_1,lam_2,l)

      result[l]=t1=_mm_mul_pd(lam_1,corfac);
      if (v2df_any_gt(t1,eps2))
        { v2df tmp=lam_1;lam_1=lam_2;lam_2=tmp; break; }
      if (++l>lmax) break;
      NEXTSTEP(pre2,pre3,pre0,pre1,lam_2,lam_1,l)

      if (v2df_any_gt(lam_2,fbig2))
        RENORMALIZE;
      }
    }

  gen->firstl[0]=l;
  if (l>lmax)
    {
    gen->m_crit=m;
    gen->cth_crit= (fabs(cth1)<fabs(cth2)) ? fabs(cth1) : fabs(cth2);
    return;
    }

  GETPRE(pre0,pre1,l)
  while (1)
    {
    v2df t1;
    result[l]=t1=_mm_mul_pd(lam_2,corfac);
    if (v2df_all_ge(t1,eps2))
      break;
    if (++l>lmax) return;
    NEXTSTEP(pre0,pre1,pre2,pre3,lam_1,lam_2,l)

    result[l]=t1=_mm_mul_pd(lam_1,corfac);
    if (v2df_all_ge(t1,eps2))
      { v2df tmp=lam_1;lam_1=lam_2;lam_2=tmp; break; }
    if (++l>lmax) return;
    NEXTSTEP(pre2,pre3,pre0,pre1,lam_2,lam_1,l)

    if (v2df_any_gt(lam_2,fbig2))
      RENORMALIZE;
    }

  lam_1 = _mm_mul_pd (lam_1,corfac);
  lam_2 = _mm_mul_pd (lam_2,corfac);

  GETPRE(pre0,pre1,l)
  for(;l<lmax-2;l+=2)
    {
    result[l]=lam_2;
    NEXTSTEP(pre0,pre1,pre2,pre3,lam_1,lam_2,l+1)
    result[l+1]=lam_1;
    NEXTSTEP(pre2,pre3,pre0,pre1,lam_2,lam_1,l+2)
    }

  while (1)
    {
    result[l]=lam_2;
    if (++l>lmax) break;
    NEXTSTEP(pre0,pre1,pre2,pre3,lam_1,lam_2,l)
    result[l] = lam_1;
    if (++l>lmax) break;
    NEXTSTEP(pre2,pre3,pre0,pre1,lam_2,lam_1,l)
    }
  }

static void Ylmgen_recalc_lambda_wx1_sse2 (Ylmgen_C *gen)
  {
  if (gen->lwx_uptodate_sse2[1]) return;
  Ylmgen_recalc_Ylm_sse2(gen);
  gen->firstl[1] = gen->firstl[0];
  if (gen->firstl[1]>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate_sse2[1] = 1;

  {
  v2df cth=build_v2df(gen->cth[gen->ith1],gen->cth[gen->ith2]);
  v2df xsth=build_v2df(1./gen->sth[gen->ith1],1./gen->sth[gen->ith2]);
  v2df m=build_v2df(gen->m_cur,gen->m_cur);
  v2df m_on_sth = _mm_mul_pd(m,xsth);
  v2df lam_lm=_mm_setzero_pd();
  v2df2 *lambda_wx = gen->lambda_wx_sse2[1];
  int l;
  v2df ell=build_v2df(gen->firstl[1],gen->firstl[1]);
  v2df uno=_mm_set1_pd(1.);
  for (l=gen->firstl[1]; l<=gen->lmax; ++l, ell=_mm_add_pd(ell,uno))
    {
    v2df lamfact=_mm_load1_pd(&gen->lamfact[l]);
    v2df lam_lm1m=lam_lm;
    lam_lm=gen->ylm_sse2[l];
    lambda_wx[l].a = _mm_mul_pd(xsth,_mm_sub_pd(_mm_mul_pd(lamfact,lam_lm1m),
                     _mm_mul_pd(_mm_mul_pd(ell,cth),lam_lm)));
    lambda_wx[l].b = _mm_mul_pd(m_on_sth,lam_lm);
    }
  }
  }

static void Ylmgen_recalc_lambda_wx2_sse2 (Ylmgen_C *gen)
  {
  if (gen->lwx_uptodate_sse2[2]) return;
  Ylmgen_recalc_Ylm_sse2(gen);
  gen->firstl[2] = gen->firstl[0];
  if (gen->firstl[2]>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate_sse2[2] = 1;

  {
  v2df cth=build_v2df(gen->cth[gen->ith1],gen->cth[gen->ith2]);
  v2df sth=build_v2df(gen->sth[gen->ith1],gen->sth[gen->ith2]);
  v2df m=build_v2df(gen->m_cur,gen->m_cur);
  v2df uno=_mm_set1_pd(1.);
  v2df one_on_s2 = _mm_div_pd(uno,_mm_mul_pd(sth,sth));
  v2df two_on_s2 = _mm_mul_pd(_mm_set1_pd(2.),one_on_s2);
  v2df two_c_on_s2 = _mm_mul_pd(cth,two_on_s2);
  v2df m2 = _mm_mul_pd(m,m);
  v2df two_m_on_s2 = _mm_mul_pd(m,two_on_s2);
  v2df lam_lm=_mm_setzero_pd();
  v2df2 *lambda_wx = gen->lambda_wx_sse2[2];
  int l;
  v2df ell=build_v2df(gen->firstl[2],gen->firstl[2]);
  for (l=gen->firstl[2]; l<=gen->lmax; ++l, ell=_mm_add_pd(ell,uno))
    {
    v2df lamfact=_mm_load1_pd(&gen->lamfact[l]);
    v2df lam_lm1m=lam_lm;
    lam_lm=gen->ylm_sse2[l];
    {
    const v2df t1  = _mm_mul_pd(lam_lm1m,lamfact);
    const v2df ellm1 = _mm_sub_pd(ell,uno);
    const v2df a_w = _mm_sub_pd
      (_mm_mul_pd(_mm_sub_pd(m2,ell),two_on_s2),_mm_mul_pd(ell,ellm1));
    const v2df a_x = _mm_mul_pd(_mm_mul_pd(cth,ellm1),lam_lm);
    lambda_wx[l].a =
      _mm_add_pd(_mm_mul_pd(a_w,lam_lm),_mm_mul_pd(t1,two_c_on_s2));
    lambda_wx[l].b = _mm_mul_pd(two_m_on_s2,_mm_sub_pd(t1,a_x));
    }
    }
  }
  }

void Ylmgen_recalc_lambda_wx_sse2 (Ylmgen_C *gen, int spin)
  {
  UTIL_ASSERT ((spin>0) && (spin<=max_spin),
    "invalid spin in Ylmgen_recalc_lambda_wx_sse2");

  if (!gen->lambda_wx_sse2[spin])
    gen->lambda_wx_sse2[spin]=RALLOC(v2df2,gen->lmax+1);

  if (gen->spinrec && spin==1) { Ylmgen_recalc_lambda_wx1_sse2(gen); return; }
  if (gen->spinrec && spin==2) { Ylmgen_recalc_lambda_wx2_sse2(gen); return; }

  if (!gen->sylm[spin])
    {
    gen->sylm[spin]=RALLOC(sylmgen_d,1);
    sylmgen_init(gen->sylm[spin],gen,spin);
    }
  if (gen->lwx_uptodate_sse2[spin]) return;
  sylmgen_prepare(gen->sylm[spin],gen,gen->m_cur);
  sylmgen_recalc_sse2(gen->sylm[spin],gen,gen->ith1,gen->ith2,
    gen->lambda_wx_sse2[spin],&gen->firstl[spin]);
  gen->lwx_uptodate_sse2[spin] = 1;
  }

#endif /* PLANCK_HAVE_SSE2 */

double *Ylmgen_get_norm (int lmax, int spin, int spinrec)
  {
  const double pi = 3.141592653589793238462643383279502884197;
  double *res=RALLOC(double,lmax+1);
  int l;
  double spinsign;
  /* sign convention for H=1 (LensPix paper) */
#if 1
  spinsign = (spin>0) ? -1.0 : 1.0;
#else
  spinsign = 1.0;
#endif

  if (spin==0)
    {
    for (l=0; l<=lmax; ++l)
      res[l]=1.;
    return res;
    }

  if ((!spinrec) || (spin>=3))
    {
    spinsign = (spin&1) ? -spinsign : spinsign;
    for (l=0; l<=lmax; ++l)
      res[l] = (l<spin) ? 0. : spinsign*0.5*sqrt((2*l+1)/(4*pi));
    return res;
    }

  if (spin==1)
    {
    for (l=0; l<=lmax; ++l)
      res[l] = (l<spin) ? 0. : spinsign*sqrt(1./((l+1.)*l));
    return res;
    }

  if (spin==2)
    {
    for (l=0; l<=lmax; ++l)
      res[l] = (l<spin) ? 0. : spinsign*sqrt(1./((l+2.)*(l+1.)*l*(l-1.)));
    return res;
    }

  UTIL_FAIL ("error in Ylmgen_get_norm");
  return NULL;
  }

int Ylmgen_maxspin(void)
  { return max_spin; }

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
 *  Copyright (C) 2005-2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "ylmgen_c.h"
#include "c_utils.h"

enum { large_exponent2=90, minscale=-4, maxscale=11, max_spin=3 };

void Ylmgen_init (Ylmgen_C *gen, int l_max, int m_max, double epsilon)
  {
  int m;
  const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;
  const double inv_ln2 = 1.4426950408889634073599246810018921;

  gen->fsmall = ldexp(1.,-large_exponent2);
  gen->fbig   = ldexp(1., large_exponent2);
  gen->eps = epsilon;
  gen->cth_crit = 2.;
  gen->ith = -1;
  gen->lmax = l_max;
  gen->mmax = m_max;
  gen->m_cur = -1;
  gen->m_crit = gen->mmax+1;
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
  ALLOC2D(gen->lambda_wx,ylmgen_dbl2,max_spin+1,gen->lmax+1);

  gen->ylm_uptodate = 0;
  gen->lwx_uptodate = RALLOC(int,max_spin+1);
  SET_ARRAY(gen->lwx_uptodate,0,max_spin+1,0);
  gen->recfac_uptodate = 0;
  gen->lamfact_uptodate = 0;

  gen->cth = gen->sth = gen->logsth = NULL;

#ifdef PLANCK_HAVE_SSE2
  gen->ith1 = gen->ith2 = -1;
  gen->ylm_sse2 = RALLOC(v2df,gen->lmax+1);
  ALLOC2D(gen->lambda_wx_sse2,v2df2,max_spin+1,gen->lmax+1);
  gen->ylm_uptodate_sse2 = 0;
  gen->lwx_uptodate_sse2 = RALLOC(int,max_spin+1);
  SET_ARRAY(gen->lwx_uptodate_sse2,0,max_spin+1,0);
#endif
  }

void Ylmgen_destroy (Ylmgen_C *gen)
  {
  DEALLOC(gen->cf);
  DEALLOC(gen->recfac);
  DEALLOC(gen->mfac);
  DEALLOC(gen->t1fac);
  DEALLOC(gen->t2fac);
  DEALLOC(gen->lamfact);
  DEALLOC(gen->ylm);
  DEALLOC(gen->lwx_uptodate);
  DEALLOC2D(gen->lambda_wx);
  DEALLOC(gen->cth);
  DEALLOC(gen->sth);
  DEALLOC(gen->logsth);
#ifdef PLANCK_HAVE_SSE2
  DEALLOC(gen->ylm_sse2);
  DEALLOC2D(gen->lambda_wx_sse2);
  DEALLOC(gen->lwx_uptodate_sse2);
#endif
  }

void Ylmgen_set_theta (Ylmgen_C *gen, const double *theta, int nth)
  {
  const double inv_ln2 = 1.4426950408889634073599246810018921;
  int m;
  DEALLOC(gen->cth);
  DEALLOC(gen->sth);
  DEALLOC(gen->logsth);
  gen->cth = RALLOC(double,nth);
  gen->sth = RALLOC(double,nth);
  gen->logsth = RALLOC(double,nth);
  for (m=0; m<nth; ++m)
    {
    gen->cth[m] = cos(theta[m]);
    gen->sth[m] = sin(theta[m]);
    /* if ring is exactly on a pole, move it a bit */
    if (fabs(gen->sth[m])<1e-15) gen->sth[m]=1e-15;
    gen->logsth[m] = inv_ln2*log(gen->sth[m]);
    UTIL_ASSERT(gen->sth[m]>0,"bad theta angle specified");
    }

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
    { gen->firstl=gen->lmax+1; return; }

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

  gen->firstl=l;
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
  if (gen->firstl>gen->lmax) return;
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
  for (ell=l=gen->firstl; l<=gen->lmax; ++l, ell+=1.)
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
  if (gen->firstl>gen->lmax) return;
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
  for (ell=l=gen->firstl; l<=gen->lmax; ++l, ell+=1.)
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

static void Ylmgen_recalc_lambda_wx_recursive (Ylmgen_C *gen, int spin)
  {
  if (gen->lwx_uptodate[spin]) return;
  Ylmgen_recalc_lambda_wx (gen, spin-1);
  if (gen->firstl>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate[spin] = 1;

  {
  ylmgen_dbl2 *lwxm1 = gen->lambda_wx[spin-1];
  ylmgen_dbl2 *lambda_wx = gen->lambda_wx[spin];
  double m=gen->m_cur;
  double cth=gen->cth[gen->ith];
  double sth=gen->sth[gen->ith];
  int ifirst = IMAX(1,gen->firstl);
  double spm1 = spin-1.;
  double last_w=0, last_x=0;

  int l;
  double ell;
  lambda_wx[0][0] = lambda_wx[0][1] = 0;
  for (ell=l=ifirst; l<=gen->lmax; ++l, ell+=1.)
    {
    double xlsth = 1./(ell*sth);
    double lcth = ell*cth;
    lambda_wx[l][0] = ( (l-spm1)*(m*lwxm1[l][1] - lcth*lwxm1[l][0])
                      + (l+spm1)*gen->lamfact[l]*last_w)
                    * xlsth;
    lambda_wx[l][1] = ( (l-spm1)*(m*lwxm1[l][0] - lcth*lwxm1[l][1])
                      + (l+spm1)*gen->lamfact[l]*last_x)
                    * xlsth;
    last_w = lwxm1[l][0];
    last_x = lwxm1[l][1];
    }
  }
  }

void Ylmgen_recalc_lambda_wx (Ylmgen_C *gen, int spin)
  {
  UTIL_ASSERT ((spin>0) && (spin<=max_spin),
    "invalid spin in Ylmgen_recalc_lambda_wx");
  switch(spin)
    {
    case 1 : Ylmgen_recalc_lambda_wx1(gen); break;
    case 2 : Ylmgen_recalc_lambda_wx2(gen); break;
    default:
      Ylmgen_recalc_lambda_wx_recursive(gen,spin);
      break;
    }
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
    { gen->firstl=gen->lmax+1; return; }

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

  gen->firstl=l;
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
  if (gen->firstl>gen->lmax) return;
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
  v2df ell=build_v2df(gen->firstl,gen->firstl);
  v2df uno=_mm_set1_pd(1.);
  for (l=gen->firstl; l<=gen->lmax; ++l, ell=_mm_add_pd(ell,uno))
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
  if (gen->firstl>gen->lmax) return;
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
  v2df ell=build_v2df(gen->firstl,gen->firstl);
  for (l=gen->firstl; l<=gen->lmax; ++l, ell=_mm_add_pd(ell,uno))
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

static void Ylmgen_recalc_lambda_wx_recursive_sse2 (Ylmgen_C *gen, int spin)
  {
  if (gen->lwx_uptodate_sse2[spin]) return;
  Ylmgen_recalc_lambda_wx_sse2 (gen, spin-1);
  if (gen->firstl>gen->lmax) return;
  Ylmgen_recalc_lamfact(gen);
  gen->lwx_uptodate_sse2[spin] = 1;

  {
  v2df2 *lwxm1 = gen->lambda_wx_sse2[spin-1];
  v2df2 *lambda_wx = gen->lambda_wx_sse2[spin];
  v2df m=build_v2df(gen->m_cur,gen->m_cur);
  v2df uno=_mm_set1_pd(1.);
  v2df cth=build_v2df(gen->cth[gen->ith1],gen->cth[gen->ith2]);
  v2df sth=build_v2df(gen->sth[gen->ith1],gen->sth[gen->ith2]);
  int ifirst = IMAX(1,gen->firstl);
  v2df spm1 = build_v2df(spin-1.,spin-1.);
  v2df last_w=_mm_setzero_pd(), last_x=_mm_setzero_pd();

  int l;
  v2df ell=build_v2df(ifirst,ifirst);
  lambda_wx[0].a = lambda_wx[0].b = _mm_setzero_pd();
  for (l=ifirst; l<=gen->lmax; ++l, ell=_mm_add_pd(ell,uno))
    {
    v2df xlsth = _mm_div_pd(uno,_mm_mul_pd(ell,sth));
    v2df lcth = _mm_mul_pd(ell,cth);
    v2df lmspm1 = _mm_sub_pd(ell,spm1);
    v2df lpspm1 = _mm_add_pd(ell,spm1);
    v2df t0 = _mm_mul_pd(lpspm1,_mm_load1_pd(&gen->lamfact[l]));
    v2df t1 = _mm_mul_pd(m,lwxm1[l].b);
    v2df t2 = _mm_mul_pd(lcth,lwxm1[l].a);
    v2df t3 = _mm_sub_pd(t1,t2);
    v2df t4 = _mm_mul_pd(lmspm1,t3);
    v2df t5 = _mm_mul_pd(t0,last_w);
    v2df t6 = _mm_add_pd(t4,t5);
    lambda_wx[l].a = _mm_mul_pd(t6, xlsth);
    t1 = _mm_mul_pd(m,lwxm1[l].a);
    t2 = _mm_mul_pd(lcth,lwxm1[l].b);
    t3 = _mm_sub_pd(t1,t2);
    t4 = _mm_mul_pd(lmspm1,t3);
    t5 = _mm_mul_pd(t0,last_x);
    t6 = _mm_add_pd(t4,t5);
    lambda_wx[l].b = _mm_mul_pd(t6, xlsth);
    last_w = lwxm1[l].a;
    last_x = lwxm1[l].b;
    }
  }
  }

void Ylmgen_recalc_lambda_wx_sse2 (Ylmgen_C *gen, int spin)
  {
  UTIL_ASSERT ((spin>0) && (spin<=max_spin),
    "invalid spin in Ylmgen_recalc_lambda_wx_sse2");
  switch(spin)
    {
    case 1 : Ylmgen_recalc_lambda_wx1_sse2(gen); break;
    case 2 : Ylmgen_recalc_lambda_wx2_sse2(gen); break;
    default:
      Ylmgen_recalc_lambda_wx_recursive_sse2(gen,spin);
      break;
    }
  }

#endif /* PLANCK_HAVE_SSE2 */

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

/*
 *  Helper code for efficient calculation of Y_lm(theta,phi=0)
 *
 *  Copyright (C) 2005-2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "sharp_ylmgen_c.h"
#include "c_utils.h"

static inline void normalize (double *val, int *scale, double xfmax)
  {
  while (fabs(*val)>xfmax) { *val*=sharp_fsmall; ++*scale; }
  if (*val!=0.)
    while (fabs(*val)<xfmax*sharp_fsmall) { *val*=sharp_fbig; --*scale; }
  }

void sharp_Ylmgen_init (sharp_Ylmgen_C *gen, int l_max, int m_max, int spin)
  {
  const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;

  gen->lmax = l_max;
  gen->mmax = m_max;
  UTIL_ASSERT(spin>=0,"incorrect spin: must be nonnegative");
  UTIL_ASSERT(m_max>=spin,"incorrect m_max: must be >= spin");
  UTIL_ASSERT(l_max>=m_max,"incorrect l_max: must be >= m_max");
  gen->s = spin;
  UTIL_ASSERT((sharp_minscale<=0)&&(sharp_maxscale>0),
    "bad value for min/maxscale");
  gen->cf=RALLOC(double,sharp_maxscale-sharp_minscale+1);
  gen->cf[-sharp_minscale]=1.;
  for (int m=-sharp_minscale-1; m>=0; --m)
    gen->cf[m]=gen->cf[m+1]*sharp_fsmall;
  for (int m=-sharp_minscale+1; m<(sharp_maxscale-sharp_minscale+1); ++m)
    gen->cf[m]=gen->cf[m-1]*sharp_fbig;

  gen->m = -1;
  if (spin==0)
    {
    gen->rf = RALLOC(sharp_ylmgen_dbl2,gen->lmax+1);
    gen->mfac = RALLOC(double,gen->mmax+1);
    gen->mfac[0] = inv_sqrt4pi;
    for (int m=1; m<=gen->mmax; ++m)
      gen->mfac[m] = gen->mfac[m-1]*sqrt((2*m+1.)/(2*m));
    gen->root = RALLOC(double,2*gen->lmax+5);
    gen->iroot = RALLOC(double,2*gen->lmax+5);
    for (int m=0; m<2*gen->lmax+5; ++m)
      {
      gen->root[m] = sqrt(m);
      gen->iroot[m] = (m==0) ? 0. : 1./gen->root[m];
      }
    }
  else
    {
    gen->m=gen->mlo=gen->mhi=-1234567890;
    ALLOC(gen->fx,sharp_ylmgen_dbl3,gen->lmax+2);
    for (int m=0; m<gen->lmax+2; ++m)
      gen->fx[m].f[0]=gen->fx[m].f[1]=gen->fx[m].f[2]=0.;
    ALLOC(gen->inv,double,gen->lmax+1);
    gen->inv[0]=0;
    for (int m=1; m<gen->lmax+1; ++m) gen->inv[m]=1./m;
    ALLOC(gen->flm1,double,2*gen->lmax+1);
    ALLOC(gen->flm2,double,2*gen->lmax+1);
    for (int m=0; m<2*gen->lmax+1; ++m)
      {
      gen->flm1[m] = sqrt(1./(m+1.));
      gen->flm2[m] = sqrt(m/(m+1.));
      }
    ALLOC(gen->prefac,double,gen->mmax+1);
    ALLOC(gen->fscale,int,gen->mmax+1);
    double *fac = RALLOC(double,2*gen->lmax+1);
    int *facscale = RALLOC(int,2*gen->lmax+1);
    fac[0]=1; facscale[0]=0;
    for (int m=1; m<2*gen->lmax+1; ++m)
      {
      fac[m]=fac[m-1]*sqrt(m);
      facscale[m]=facscale[m-1];
      normalize(&fac[m],&facscale[m],sharp_fbighalf);
      }
    for (int m=0; m<=gen->mmax; ++m)
      {
      int mlo=gen->s, mhi=m;
      if (mhi<mlo) SWAP(mhi,mlo,int);
      double tfac=fac[2*mhi]/fac[mhi+mlo];
      int tscale=facscale[2*mhi]-facscale[mhi+mlo];
      normalize(&tfac,&tscale,sharp_fbighalf);
      tfac/=fac[mhi-mlo];
      tscale-=facscale[mhi-mlo];
      normalize(&tfac,&tscale,sharp_fbighalf);
      gen->prefac[m]=tfac;
      gen->fscale[m]=tscale;
      }
    DEALLOC(fac);
    DEALLOC(facscale);
    }
  }

void sharp_Ylmgen_destroy (sharp_Ylmgen_C *gen)
  {
  DEALLOC(gen->cf);
  if (gen->s==0)
    {
    DEALLOC(gen->rf);
    DEALLOC(gen->mfac);
    DEALLOC(gen->root);
    DEALLOC(gen->iroot);
    }
  else
    {
    DEALLOC(gen->fx);
    DEALLOC(gen->prefac);
    DEALLOC(gen->fscale);
    DEALLOC(gen->flm1);
    DEALLOC(gen->flm2);
    DEALLOC(gen->inv);
    }
  }

void sharp_Ylmgen_prepare (sharp_Ylmgen_C *gen, int m)
  {
  if (m==gen->m) return;
  UTIL_ASSERT(m>=0,"incorrect m");
  gen->m = m;

  if (gen->s==0)
    {
    gen->rf[m].f[0] = gen->root[2*m+3];
    gen->rf[m].f[1] = 0.;
    for (int l=m+1; l<=gen->lmax; ++l)
      {
      double tmp=gen->root[2*l+3]*gen->iroot[l+1+m]*gen->iroot[l+1-m];
      gen->rf[l].f[0] = tmp*gen->root[2*l+1];
      gen->rf[l].f[1] = tmp*gen->root[l+m]*gen->root[l-m]*gen->iroot[2*l-1];
      }
    }
  else
    {
    int mlo_=m, mhi_=gen->s;
    if (mhi_<mlo_) SWAP(mhi_,mlo_,int);
    int ms_similar = ((gen->mhi==mhi_) && (gen->mlo==mlo_));

    gen->mlo = mlo_; gen->mhi = mhi_;

    if (!ms_similar)
      {
      for (int l=gen->mhi; l<gen->lmax; ++l)
        {
        double t = gen->flm1[l+gen->m]*gen->flm1[l-gen->m]
                  *gen->flm1[l+gen->s]*gen->flm1[l-gen->s];
        double lt = 2*l+1;
        double l1 = l+1;
        gen->fx[l+1].f[0]=l1*lt*t;
        gen->fx[l+1].f[1]=gen->m*gen->s*gen->inv[l]*gen->inv[l+1];
        t = gen->flm2[l+gen->m]*gen->flm2[l-gen->m]
           *gen->flm2[l+gen->s]*gen->flm2[l-gen->s];
        gen->fx[l+1].f[2]=t*l1*gen->inv[l];
        }
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
  }

double *sharp_Ylmgen_get_norm (int lmax, int spin)
  {
  const double pi = 3.141592653589793238462643383279502884197;
  double *res=RALLOC(double,lmax+1);
  /* sign convention for H=1 (LensPix paper) */
#if 1
   double spinsign = (spin>0) ? -1.0 : 1.0;
#else
   double spinsign = 1.0;
#endif

  if (spin==0)
    {
    for (int l=0; l<=lmax; ++l)
      res[l]=1.;
    return res;
    }

  spinsign = (spin&1) ? -spinsign : spinsign;
  for (int l=0; l<=lmax; ++l)
    res[l] = (l<spin) ? 0. : spinsign*0.5*sqrt((2*l+1)/(4*pi));
  return res;
  }

double *sharp_Ylmgen_get_d1norm (int lmax)
  {
  const double pi = 3.141592653589793238462643383279502884197;
  double *res=RALLOC(double,lmax+1);

  for (int l=0; l<=lmax; ++l)
    res[l] = (l<1) ? 0. : 0.5*sqrt(l*(l+1.)*(2*l+1.)/(4*pi));
  return res;
  }

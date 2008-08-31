/*
 *  This file is part of libfftpack.
 *
 *  libfftpack is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libfftpack is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libfftpack; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libfftpack is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include <stdlib.h>
#include "fftpack.h"
#include "bluestein.h"

/* returns the sum of all prime factors of n */
int prime_factor_sum (int n)
  {
  int result=0,x,limit,tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  limit=sqrt(n+0.01);
  for (x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+=x;
    n=tmp;
    limit=sqrt(n+0.01);
    }
  if (n>1) result+=n;

  return result;
  }

/* returns the smallest composite of 2, 3 and 5 which is >= n */
static int good_size(int n)
  {
  int maxfactors=1, i, j, k, f2=1, f3, f5, bestfac, guessfac;
  while ((n>>maxfactors)>0)
    ++maxfactors;
  bestfac=1<<maxfactors;

  for (i=0; i<maxfactors; ++i)
    {
    f3=1;
    for (j=0; j<maxfactors-i; ++j)
      {
      if (f2*f3>bestfac) break;
      f5=1;
      for (k=0; k<maxfactors-i-j; ++k)
        {
        guessfac = f2*f3*f5;
        if (guessfac>=bestfac) break;
        if ((guessfac>=n) && (guessfac<bestfac))
          bestfac=guessfac;
        f5*=5;
        }
      f3*=3;
      }
    f2*=2;
    }
  return bestfac;
  }

void bluestein_i (int n, double **tstorage)
  {
  static const double pi=3.14159265358979323846;
  int n2=good_size(n*2-1);
  int m, coeff;
  double angle, xn2;
  double *bk, *bkf, *work;
  double pibyn=pi/n;
  *tstorage = (double *)malloc (sizeof(double)*(1+2*n+8*n2+15));
  ((int *)(*tstorage))[0]=n2;
  bk  = *tstorage+1;
  bkf = *tstorage+1+2*n;
  work= *tstorage+1+2*(n+n2);

/* initialize b_k */
  bk[0] = 1;
  bk[1] = 0;

  coeff=0;
  for (m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    angle = pibyn*coeff;
    bk[2*m] = cos(angle);
    bk[2*m+1] = sin(angle);
    }

/* initialize the zero-padded, Fourier transformed b_k. Add normalisation.  */
  xn2 = 1./n2;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
    }
  for (m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;
  cffti (n2,work);
  cfftf (n2,bkf,work);
  }

void bluestein (int n, double *data, double *tstorage, int isign)
  {
  int n2=*((int *)tstorage);
  int m;
  double *bk, *bkf, *akf, *work;
  bk  = tstorage+1;
  bkf = tstorage+1+2*n;
  work= tstorage+1+2*(n+n2);
  akf = tstorage+1+2*n+6*n2+15;  

/* initialize a_k and FFT it */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   - data[m+1]*bk[m+1];
      akf[m+1] = data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   + data[m+1]*bk[m+1];
      akf[m+1] =-data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  for (m=2*n; m<2*n2; ++m)
    akf[m]=0;

  cfftf (n2,akf,work);

/* do the convolution */
  if (isign>0)
    for (m=0; m<2*n2; m+=2)
      {
      double im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (m=0; m<2*n2; m+=2)
      {
      double im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }


/* inverse FFT */
  cfftb (n2,akf,work);

/* multiply by b_k* */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] - bk[m+1]*akf[m+1];
      data[m+1] = bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] + bk[m+1]*akf[m+1];
      data[m+1] =-bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  }

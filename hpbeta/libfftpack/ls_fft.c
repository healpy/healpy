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
 *  Copyright (C) 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bluestein.h"
#include "fftpack.h"
#include "ls_fft.h"

complex_plan make_complex_plan (int length)
  {
  complex_plan plan = (complex_plan) malloc(sizeof(complex_plan_i));
  int pfsum = prime_factor_sum(length);
  double comp1 = length*pfsum;
  double comp2 = 2*3*length*log(3.*length);
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work));
  else
    {
    plan->work=(double *)malloc((4*length+15)*sizeof(double));
    cffti(length, plan->work);
    }
  return plan;
  }

void kill_complex_plan (complex_plan plan)
  {
  free(plan->work);
  free(plan);
  }

void complex_plan_forward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, -1);
  else
    cfftf (plan->length, data, plan->work);
  }

void complex_plan_backward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, 1);
  else
    cfftb (plan->length, data, plan->work);
  }


real_plan make_real_plan (int length)
  {
  real_plan plan = (real_plan) malloc(sizeof(real_plan_i));
  int pfsum = prime_factor_sum(length);
  double comp1 = .5*length*pfsum;
  double comp2 = 2*3*length*log(3.*length);
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work));
  else
    {
    plan->work=(double *)malloc((2*length+15)*sizeof(double));
    rffti(length, plan->work);
    }
  return plan;
  }

void kill_real_plan (real_plan plan)
  {
  free(plan->work);
  free(plan);
  }

void real_plan_forward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    int m;
    int n=plan->length;
    double *tmp = (double *)malloc(2*n*sizeof(double));
    for (m=0; m<n; ++m)
      {
      tmp[2*m] = data[m];
      tmp[2*m+1] = 0.;
      }
    bluestein(n,tmp,plan->work,-1);
    data[0] = tmp[0];
    memcpy (data+1, tmp+2, (n-1)*sizeof(double));
    free (tmp);
    }
  else
    rfftf (plan->length, data, plan->work);
  }

static void fftpack2halfcomplex (double *data, int n)
  {
  int m;
  double *tmp = (double *)malloc(n*sizeof(double));
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[m]=data[2*m-1];
    tmp[n-m]=data[2*m];
    }
  if (!(n&1))
    tmp[n/2]=data[n-1];
  memcpy (data,tmp,n*sizeof(double));
  free(tmp);
  }

static void halfcomplex2fftpack (double *data, int n)
  {
  int m;
  double *tmp = (double *)malloc(n*sizeof(double));
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[2*m-1]=data[m];
    tmp[2*m]=data[n-m];
    }
  if (!(n&1))
    tmp[n-1]=data[n/2];
  memcpy (data,tmp,n*sizeof(double));
  free(tmp);
  }

void real_plan_forward_fftw (real_plan plan, double *data)
  {
  real_plan_forward_fftpack (plan, data);
  fftpack2halfcomplex (data,plan->length);
  }

void real_plan_backward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    int m;
    int n=plan->length;
    double *tmp = (double *)malloc(2*n*sizeof(double));
    tmp[0]=data[0];
    tmp[1]=0.;
    memcpy (tmp+2,data+1, (n-1)*sizeof(double));
    if ((n&1)==0) tmp[n+1]=0.;
    for (m=2; m<n; m+=2)
      {
      tmp[2*n-m]=tmp[m];
      tmp[2*n-m+1]=-tmp[m+1];
      }
    bluestein (n, tmp, plan->work, 1);
    for (m=0; m<n; ++m)
      data[m] = tmp[2*m];
    free (tmp);
    }
  else
    rfftb (plan->length, data, plan->work);
  }

void real_plan_backward_fftw (real_plan plan, double *data)
  {
  halfcomplex2fftpack (data,plan->length);
  real_plan_backward_fftpack (plan, data);
  }

void real_plan_forward_c (real_plan plan, double *data)
  {
  int m;
  int n=plan->length;

  if (plan->bluestein)
    {
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    bluestein (plan->length, data, plan->work, -1);
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  else
    {
    for (m=0; m<n; ++m) data[m+1] = data[2*m];
    rfftf (n, data+1, plan->work);
    data[0] = data[1];
    data[1] = 0;
    for (m=2; m<n; m+=2)
      {
      data[2*n-m]   =  data[m];
      data[2*n-m+1] = -data[m+1];
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  }

void real_plan_backward_c (real_plan plan, double *data)
  {
  int m;
  int n=plan->length;

  if (plan->bluestein)
    {
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    bluestein (plan->length, data, plan->work, 1);
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    }
  else
    {
    data[1] = data[0];
    rfftb (n, data+1, plan->work);
    for (m=n-1; m>=0; --m)
      {
      data[2*m]   = data[m+1];
      data[2*m+1] = 0.;
      }
    }
  }

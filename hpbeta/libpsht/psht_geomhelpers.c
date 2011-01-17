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

/*! \file psht_geomhelpers.c
 *  Spherical transform library
 *
 *  Copyright (C) 2006-2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include "psht_geomhelpers.h"
#include "c_utils.h"

void psht_make_healpix_geom_info (int nside, int stride,
  psht_geom_info **geom_info)
  {
  double *weight=RALLOC(double,2*nside);
  SET_ARRAY(weight,0,2*nside,1);
  psht_make_weighted_healpix_geom_info (nside, stride, weight, geom_info);
  DEALLOC(weight);
  }

void psht_make_weighted_healpix_geom_info (int nside, int stride,
  const double *weight, psht_geom_info **geom_info)
  {
  const double pi=3.141592653589793238462643383279502884197;
  ptrdiff_t npix=(ptrdiff_t)nside*nside*12;
  ptrdiff_t ncap=2*(ptrdiff_t)nside*(nside-1);
  int nrings=4*nside-1;

  double *theta=RALLOC(double,nrings);
  double *weight_=RALLOC(double,nrings);
  int *nph=RALLOC(int,nrings);
  double *phi0=RALLOC(double,nrings);
  ptrdiff_t *ofs=RALLOC(ptrdiff_t,nrings);
  int *stride_=RALLOC(int,nrings);
  int m;
  for (m=0; m<nrings; ++m)
    {
    int ring=m+1;
    ptrdiff_t northring = (ring>2*nside) ? 4*nside-ring : ring;
    stride_[m] = stride;
    if (northring < nside)
      {
      theta[m] = 2*asin(northring/(sqrt(6.)*nside));
      nph[m] = 4*northring;
      phi0[m] = pi/nph[m];
      ofs[m] = 2*northring*(northring-1);
      }
    else
      {
      double fact1 = (8.*nside)/npix;
      double costheta = (2*nside-northring)*fact1;
      theta[m] = acos(costheta);
      nph[m] = 4*nside;
      if ((northring-nside) & 1)
        phi0[m] = 0;
      else
        phi0[m] = pi/nph[m];
      ofs[m] = ncap + (northring-nside)*nph[m];
      }
    if (northring != ring) /* southern hemisphere */
      {
      theta[m] = pi-theta[m];
      ofs[m] = npix - nph[m] - ofs[m];
      }
    weight_[m]=4.*pi/npix*weight[northring-1];
    }

  psht_make_geom_info (nrings, nph, ofs, stride_, phi0, theta, weight_,
    geom_info);

  DEALLOC(theta);
  DEALLOC(weight_);
  DEALLOC(nph);
  DEALLOC(phi0);
  DEALLOC(ofs);
  DEALLOC(stride_);
  }

static void gauleg (double x1, double x2, double *x, double *w, int n)
  {
  const double pi = 3.141592653589793238462643383279502884197;
  const double eps = 3.0E-14;

  int m = (n+1)/2;
  double xm = 0.5*(x2+x1);
  double xl = 0.5*(x2-x1);
  int i;
  for(i=1; i<=m; ++i)
    {
    double z = cos(pi*(i-0.25)/(n+0.5));
    double pp;
    int dobreak=0;
    while(1)
      {
      double p1 = 1.0, p2 = 0.0;
      double z1 = z;
      int j;
      for(j=1; j<=n; ++j)
        {
        double p3 = p2;
        p2 = p1;
        p1 = ((2*j-1)*z*p2-(j-1)*p3)/j;
        }
      pp = n*(z*p1-p2)/(z*z-1);
      z = z1 - p1/pp;
      if (dobreak) break;
      if (fabs(z-z1) <= eps) dobreak=1;
      }
    x[i-1] = xm - xl*z;
    x[n-i] = xm + xl*z;
    w[i-1] = w[n-i] = 2*xl/((1-z*z)*pp*pp);
    }
  }

static void makeweights (int bw, double *weights)
  {
  const double pi = 3.141592653589793238462643383279502884197;
  const double fudge = pi/(4*bw);
  int j;
  for (j=0; j<2*bw; ++j)
    {
    double tmpsum = 0;
    int k;
    for (k=0; k<bw; ++k)
      tmpsum += 1./(2*k+1) * sin((2*j+1)*(2*k+1)*fudge);
    tmpsum *= sin((2*j+1)*fudge);
    tmpsum *= 2./bw;
    weights[j] = tmpsum ;
    /* weights[j + 2*bw] = tmpsum * sin((2*j+1)*fudge); */
    }
  }

void psht_make_gauss_geom_info (int nrings, int nphi, int stride,
  psht_geom_info **geom_info)
  {
  const double pi=3.141592653589793238462643383279502884197;

  double *theta=RALLOC(double,nrings);
  double *weight=RALLOC(double,nrings);
  int *nph=RALLOC(int,nrings);
  double *phi0=RALLOC(double,nrings);
  ptrdiff_t *ofs=RALLOC(ptrdiff_t,nrings);
  int *stride_=RALLOC(int,nrings);
  int m;

  gauleg(-1,1,theta,weight,nrings);

  for (m=0; m<nrings; ++m)
    {
    theta[m] = acos(theta[m]);
    nph[m]=nphi;
    phi0[m]=0;
    ofs[m]=(ptrdiff_t)m*nphi;
    stride_[m]=stride;
    weight[m]*=2*pi/nphi;
    }

  psht_make_geom_info (nrings, nph, ofs, stride_, phi0, theta, weight,
    geom_info);

  DEALLOC(theta);
  DEALLOC(weight);
  DEALLOC(nph);
  DEALLOC(phi0);
  DEALLOC(ofs);
  DEALLOC(stride_);
  }

void psht_make_ecp_geom_info (int nrings, int nphi, double phi0, int stride,
  psht_geom_info **geom_info)
  {
  const double pi=3.141592653589793238462643383279502884197;

  double *theta=RALLOC(double,nrings);
  double *weight=RALLOC(double,nrings);
  int *nph=RALLOC(int,nrings);
  double *phi0_=RALLOC(double,nrings);
  ptrdiff_t *ofs=RALLOC(ptrdiff_t,nrings);
  int *stride_=RALLOC(int,nrings);

  int m;

  UTIL_ASSERT((nrings&1)==0,
    "Even number of rings needed for equidistant grid!");
  makeweights(nrings/2,weight);
  for (m=0; m<nrings; ++m)
    {
    theta[m] = (m+0.5)*pi/nrings;
    nph[m]=nphi;
    phi0_[m]=phi0;
    ofs[m]=(ptrdiff_t)m*nphi;
    stride_[m]=stride;
    weight[m]*=2*pi/nphi;
    }

  psht_make_geom_info (nrings, nph, ofs, stride_, phi0_, theta, weight,
    geom_info);

  DEALLOC(theta);
  DEALLOC(weight);
  DEALLOC(nph);
  DEALLOC(phi0_);
  DEALLOC(ofs);
  DEALLOC(stride_);
  }

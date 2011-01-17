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

/*! \file psht_test.c
    Accuracy test for libpsht's map analysis.

    This program first generates a_lm coefficients up to
    a user-specified lmax (with mmax=lmax); where applicable, the
    real and imaginary parts of the coefficients are uniform
    random numbers of the interval [-1;1[.
    Afterwards, the random a_lm are converted to a map.
    This map is analyzed (optionally using an iterative scheme
    with a user-supplied number of steps).
    After every iteration, the code then outputs the RMS of the residual a_lm
    (i.e. the difference between the current and original a_lm), divided by
    the RMS of the original a_lm, as well as the maximum absolute change of any
    real or imaginary part between the current and original a_lm.

    This operation can be performed for several different pixelisations:
      - a Gaussian with the minimal number of rings for exact analysis
        and a user-defined ring resolution
      - an ECP grid with the minimal number of rings for exact analysis
        and a user-defined ring resolution
      - a Healpix grid with a user-defined Nside parameter.

    The user can specify the spin of the desired transform.

    Copyright (C) 2006-2010 Max-Planck-Society
    \author Martin Reinecke
*/

#include <stdio.h>
#include <string.h>
#include "psht.h"
#include "psht_geomhelpers.h"
#include "psht_almhelpers.h"
#include "c_utils.h"
#include "walltime_c.h"

static double drand (double min, double max)
  {
  return min + (max-min)*rand()/(RAND_MAX+1.0);
  }

static void random_alm (pshtd_cmplx *alm, psht_alm_info *helper, int spin)
  {
  int l,m;
  for (m=0;m<=helper->mmax; ++m)
    for (l=m;l<=helper->lmax; ++l)
      {
      if ((l<spin)&&(m<spin))
        alm[psht_alm_index(helper,l,m)] = pshtd_cmplx_null;
      else
        {
        alm[psht_alm_index(helper,l,m)].re = drand(-1,1);
        alm[psht_alm_index(helper,l,m)].im = (m==0) ? 0 : drand(-1,1);
        }
      }
  }

static void measure_errors (pshtd_cmplx **alm, pshtd_cmplx **alm2,
  ptrdiff_t nalms, int ncomp)
  {
  int i;
  ptrdiff_t m;

  for (i=0; i<ncomp; ++i)
    {
    double sum=0, sum2=0, maxdiff=0;
    for (m=0; m<nalms; ++m)
      {
      double x=alm[i][m].re-alm2[i][m].re, y=alm[i][m].im-alm2[i][m].im;
      sum+=x*x+y*y;
      sum2+=alm[i][m].re*alm[i][m].re+alm[i][m].im*alm[i][m].im;
      if (fabs(x)>maxdiff) maxdiff=fabs(x);
      if (fabs(y)>maxdiff) maxdiff=fabs(y);
      }
    sum=sqrt(sum/nalms);
    sum2=sqrt(sum2/nalms);
    printf("component %i: rms %e, maxerr %e\n",i, sum/sum2, maxdiff);
    }
  }

static void map2alm_iter (psht_geom_info *tinfo, double **map,
  pshtd_cmplx **alm_orig, pshtd_cmplx **alm, int lmax, int mmax,
  ptrdiff_t npix, ptrdiff_t nalms, int spin, int niter)
  {
  psht_alm_info *alms;
  pshtd_joblist *joblist;
  int ncomp = (spin==0) ? 1 : 2;
  int iter,i;
  ptrdiff_t m;
  double timer;

  psht_make_triangular_alm_info(lmax,mmax,1,&alms);
  pshtd_make_joblist (&joblist);

  if (spin==0)
    pshtd_add_job_map2alm(joblist,map[0],alm[0],0);
  else
    pshtd_add_job_map2alm_spin(joblist,map[0],map[1],alm[0],alm[1],spin,0);
  timer=wallTime();
  pshtd_execute_jobs (joblist, tinfo, alms);
  printf("wall time for map2alm: %fs\n",wallTime()-timer);
  pshtd_clear_joblist (joblist);
  measure_errors(alm_orig,alm,nalms,ncomp);

  for (iter=0; iter<niter; ++iter)
    {
    double **map2;
    ALLOC2D(map2,double,ncomp,npix);
    printf ("\niteration %i:\n", iter+1);
    if (spin==0)
      pshtd_add_job_alm2map(joblist,alm[0],map2[0],0);
    else
      pshtd_add_job_alm2map_spin(joblist,alm[0],alm[1],map2[0],map2[1],spin,0);
    timer=wallTime();
    pshtd_execute_jobs (joblist, tinfo, alms);
    printf("wall time for alm2map: %fs\n",wallTime()-timer);
    pshtd_clear_joblist (joblist);
    for (i=0; i<ncomp; ++i)
      for (m=0; m<npix; ++m)
        map2[i][m] = map[i][m]-map2[i][m];

    if (spin==0)
      pshtd_add_job_map2alm(joblist,map2[0],alm[0],1);
    else
      pshtd_add_job_map2alm_spin(joblist,map2[0],map2[1],alm[0],alm[1],spin,1);
    timer=wallTime();
    pshtd_execute_jobs (joblist, tinfo, alms);
    printf("wall time for map2alm: %fs\n",wallTime()-timer);
    pshtd_clear_joblist (joblist);
    DEALLOC2D(map2);
    measure_errors(alm_orig,alm,nalms,ncomp);
    }

  psht_destroy_alm_info(alms);
  pshtd_destroy_joblist(joblist);
  }

static void check_accuracy (psht_geom_info *tinfo, ptrdiff_t lmax,
  ptrdiff_t mmax, ptrdiff_t npix, int spin, int niter)
  {
  psht_alm_info *alms;
  pshtd_joblist *joblist;
  double **map;
  pshtd_cmplx **alm, **alm2;
  ptrdiff_t nalms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);
  int ncomp = (spin==0) ? 1 : 2;
  double timer;

  ALLOC2D(map,double,ncomp,npix);

  psht_make_triangular_alm_info(lmax,mmax,1,&alms);
  pshtd_make_joblist (&joblist);

  srand(4);
  ALLOC2D(alm,pshtd_cmplx,ncomp,nalms);
  random_alm(alm[0],alms,spin);
  if (spin>0)
    random_alm(alm[1],alms,spin);

  ALLOC2D(alm2,pshtd_cmplx,ncomp,nalms);

  printf ("\niteration 0:\n");
  if (spin==0)
    pshtd_add_job_alm2map(joblist,alm[0],map[0],0);
  else
    pshtd_add_job_alm2map_spin(joblist,alm[0],alm[1],map[0],map[1],spin,0);
  timer=wallTime();
  pshtd_execute_jobs (joblist, tinfo, alms);
  printf("wall time for alm2map: %fs\n",wallTime()-timer);
  pshtd_clear_joblist (joblist);

  map2alm_iter(tinfo, map, alm, alm2, lmax, mmax, npix, nalms, spin, niter);

  DEALLOC2D(map);
  DEALLOC2D(alm);
  DEALLOC2D(alm2);

  psht_destroy_alm_info(alms);
  pshtd_destroy_joblist(joblist);
  }

int main(int argc, char **argv)
  {
  int lmax;
  int spin;
  int niter;
  psht_geom_info *tinfo;

  UTIL_ASSERT (argc==6,
    "usage: psht_test <healpix|ecp|gauss> <lmax> <nside|nphi> <niter> <spin>");
  lmax=atoi(argv[2]);
  niter=atoi(argv[4]);
  spin=atoi(argv[5]);

  printf("Testing map analysis accuracy.\n");
  printf("lmax=%d, %d iterations, spin=%d\n", lmax, niter, spin);

  if (strcmp(argv[1],"gauss")==0)
    {
    int nrings=lmax+1;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting Gaussian grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    psht_make_gauss_geom_info (nrings, ppring, 1, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,niter);
    psht_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"ecp")==0)
    {
    int nrings=2*lmax+2;
    int ppring=atoi(argv[3]);
    ptrdiff_t npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting ECP grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    psht_make_ecp_geom_info (nrings, ppring, 0., 1, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,niter);
    psht_destroy_geom_info(tinfo);
    }
  else if (strcmp(argv[1],"healpix")==0)
    {
    int nside=atoi(argv[3]);
    ptrdiff_t npix;
    if (nside<1) nside=1;
    npix=12*(ptrdiff_t)nside*nside;
    printf("\nTesting Healpix grid (nside=%d, %ld pixels)\n",
          nside,(long)npix);
    psht_make_healpix_geom_info (nside, 1, &tinfo);
    check_accuracy(tinfo,lmax,lmax,npix,spin,niter);
    psht_destroy_geom_info(tinfo);
    }
  else
    UTIL_FAIL("unknown grid geometry");

  return 0;
  }

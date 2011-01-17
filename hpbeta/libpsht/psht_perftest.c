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

/*! \file psht_perftest.c
    Performance test of libpsht's map synthesis and analysis algorithms.

    This program creates a list of transform jobs with a user-specified
    maximum multipole lmax that operate on a HEALPix, ECP or Gaussian grid.
    Depending on the geometry type, the user has to specify the Nside parameter
    (for HEALPix) or the number of pixels per ring.
    The individual job types are also given by the user
    and can be "alm2map", "map2alm", "alm2map_pol", "map2alm_pol",
    "alm2map_spin[1-3]", "map2alm_spin[1-3]", and "alm2map_deriv1".
    Any combination of job types is allowed.

    All requested jobs are executed simultaneously.

    Copyright (C) 2006-2010 Max-Planck-Society
    \author Martin Reinecke
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psht.h"
#include "psht_geomhelpers.h"
#include "psht_almhelpers.h"
#include "c_utils.h"
#include "walltime_c.h"

static void get_map(float **map, ptrdiff_t npix)
  {
  *map=RALLOC(float,npix);
  SET_ARRAY(*map,0,npix,1);
  }
static void get_alm(pshts_cmplx **alm, ptrdiff_t nalm)
  {
  static const pshts_cmplx pshts_cmplx_one={1,1};
  *alm=RALLOC(pshts_cmplx,nalm);
  SET_ARRAY(*alm,0,nalm,pshts_cmplx_one);
  }

static void prepare_job (const char *jobname, float **map,
  pshts_cmplx **alm, ptrdiff_t npix, ptrdiff_t nalm, int ofs_m, int ofs_a,
  int num_m, int num_a)
  {
  int m;
  printf("adding job: %s\n", jobname);
  for (m=0; m<num_m; ++m)
    get_map (&map[ofs_m+m],npix);
  for (m=0; m<num_a; ++m)
    get_alm (&alm[ofs_a+m],nalm);
  }

int main(int argc, char **argv)
  {
  ptrdiff_t npix=0,lmax,nalm;
  float *map[100];
  pshts_cmplx *alm[100];
  psht_alm_info *alms;
  psht_geom_info *tinfo;
  pshts_joblist *joblist;
  int ofs_m, ofs_a, m;
  double wtimer;

  UTIL_ASSERT (argc>=5,
    "usage: psht_perftest <healpix|ecp|gauss> <lmax> <nside|nphi> <type>+\n"
    "  where <type> can be 'alm2map', 'map2alm', 'alm2map_pol',\n"
    "  'map2alm_pol', 'alm2map_spin[1-3]', 'map2alm_spin[1-3]',\n"
    "  or 'alm2map_deriv1'");
  lmax=atoi(argv[2]);

  if (strcmp(argv[1],"gauss")==0)
    {
    int nrings=lmax+1;
    int ppring=atoi(argv[3]);
    npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting Gaussian grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    psht_make_gauss_geom_info (nrings, ppring, 1, &tinfo);
    }
  else if (strcmp(argv[1],"ecp")==0)
    {
    int nrings=2*lmax+2;
    int ppring=atoi(argv[3]);
    npix=(ptrdiff_t)nrings*ppring;
    printf("\nTesting ECP grid (%d rings, %d pixels/ring, %ld pixels)\n",
          nrings,ppring,(long)npix);
    psht_make_ecp_geom_info (nrings, ppring, 0., 1, &tinfo);
    }
  else if (strcmp(argv[1],"healpix")==0)
    {
    int nside=atoi(argv[3]);
    if (nside<1) nside=1;
    npix=12*(ptrdiff_t)nside*nside;
    printf("\nTesting Healpix grid (nside=%d, %ld pixels)\n",
          nside,(long)npix);
    psht_make_healpix_geom_info (nside, 1, &tinfo);
    }
  else
    UTIL_FAIL("unknown command");

  psht_make_triangular_alm_info(lmax,lmax,1,&alms);
  nalm = ((ptrdiff_t)(lmax+1)*(lmax+2))/2;
  pshts_make_joblist (&joblist);

  ofs_m=ofs_a=0;
  for (m=4; m<argc; ++m)
    {
    if (strcmp(argv[m],"alm2map")==0)
      {
      prepare_job ("alm2map",map,alm,npix,nalm,ofs_m,ofs_a,1,1);
      pshts_add_job_alm2map(joblist,alm[ofs_a],map[ofs_m],0);
      ++ofs_m; ++ofs_a;
      }
    else if (strcmp(argv[m],"map2alm")==0)
      {
      prepare_job ("map2alm",map,alm,npix,nalm,ofs_m,ofs_a,1,1);
      pshts_add_job_map2alm(joblist,map[ofs_m],alm[ofs_a],0);
      ++ofs_m; ++ofs_a;
      }
    else if (strcmp(argv[m],"alm2map_pol")==0)
      {
      prepare_job ("alm2map_pol",map,alm,npix,nalm,ofs_m,ofs_a,3,3);
      pshts_add_job_alm2map_pol(joblist,alm[ofs_a],alm[ofs_a+1],alm[ofs_a+2],
                                map[ofs_m],map[ofs_m+1],map[ofs_m+2],0);
      ofs_m+=3; ofs_a+=3;
      }
    else if (strcmp(argv[m],"map2alm_pol")==0)
      {
      prepare_job ("map2alm_pol",map,alm,npix,nalm,ofs_m,ofs_a,3,3);
      pshts_add_job_map2alm_pol(joblist,map[ofs_m],map[ofs_m+1],map[ofs_m+2],
                                alm[ofs_a],alm[ofs_a+1],alm[ofs_a+2],0);
      ofs_m+=3; ofs_a+=3;
      }
    else if (strcmp(argv[m],"alm2map_spin1")==0)
      {
      prepare_job ("alm2map_spin1",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_alm2map_spin(joblist,alm[ofs_a],alm[ofs_a+1],
                                 map[ofs_m],map[ofs_m+1],1,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"map2alm_spin1")==0)
      {
      prepare_job ("map2alm_spin1",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_map2alm_spin(joblist,map[ofs_m],map[ofs_m+1],
                                 alm[ofs_a],alm[ofs_a+1],1,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"alm2map_spin2")==0)
      {
      prepare_job ("alm2map_spin2",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_alm2map_spin(joblist,alm[ofs_a],alm[ofs_a+1],
                                 map[ofs_m],map[ofs_m+1],2,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"map2alm_spin2")==0)
      {
      prepare_job ("map2alm_spin2",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_map2alm_spin(joblist,map[ofs_m],map[ofs_m+1],
                                 alm[ofs_a],alm[ofs_a+1],2,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"alm2map_spin3")==0)
      {
      prepare_job ("alm2map_spin3",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_map2alm_spin(joblist,map[ofs_m],map[ofs_m+1],
                                 alm[ofs_a],alm[ofs_a+1],3,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"map2alm_spin3")==0)
      {
      prepare_job ("map2alm_spin3",map,alm,npix,nalm,ofs_m,ofs_a,2,2);
      pshts_add_job_map2alm_spin(joblist,map[ofs_m],map[ofs_m+1],
                                 alm[ofs_a],alm[ofs_a+1],3,0);
      ofs_m+=2; ofs_a+=2;
      }
    else if (strcmp(argv[m],"alm2map_deriv1")==0)
      {
      prepare_job ("alm2map_deriv1",map,alm,npix,nalm,ofs_m,ofs_a,2,1);
      pshts_add_job_alm2map_deriv1(joblist,alm[ofs_a], map[ofs_m],
                                   map[ofs_m+1],0);
      ofs_m+=2; ofs_a+=1;
      }
    else
      UTIL_FAIL("unknown transform type");
    }

  wtimer=wallTime();
  pshts_execute_jobs (joblist, tinfo, alms);
  printf("wall time for transform: %fs\n",wallTime()-wtimer);

  pshts_destroy_joblist(joblist);
  psht_destroy_geom_info(tinfo);
  psht_destroy_alm_info(alms);
  for (m=0; m<ofs_m; ++m)
    DEALLOC(map[m]);
  for (m=0; m<ofs_a; ++m)
    DEALLOC(alm[m]);

  return 0;
  }

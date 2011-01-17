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

/*! \file psht.c
 *  Spherical transform library
 *
 *  Copyright (C) 2006-2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include "ls_fft.h"
#include "sse_utils.h"
#include "ylmgen_c.h"
#include "psht.h"
#include "c_utils.h"

const pshts_cmplx pshts_cmplx_null={0,0};
const pshtd_cmplx pshtd_cmplx_null={0,0};

static void get_chunk_info (int ndata, int *nchunks, int *chunksize)
  {
  static const int chunksize_min=100;
  static const int nchunks_max=10;
  *chunksize = IMAX(chunksize_min,(ndata+nchunks_max-1)/nchunks_max);
  if ((*chunksize)&1) ++(*chunksize);
  *nchunks = (ndata+*chunksize-1) / *chunksize;
  }

typedef struct
  {
  double phi0_;
  pshtd_cmplx *shiftarr, *work;
  int s_shift, s_work;
  real_plan plan;
  int norot;
  } ringhelper;

static void ringhelper_init (ringhelper *self)
  {
  static ringhelper rh_null = { 0, NULL, NULL, 0, 0, NULL, 0 };
  *self = rh_null;
  }

static void ringhelper_destroy (ringhelper *self)
  {
  if (self->plan) kill_real_plan(self->plan);
  DEALLOC(self->shiftarr);
  DEALLOC(self->work);
  ringhelper_init(self);
  }

static void ringhelper_update (ringhelper *self, int nph, int mmax, double phi0)
  {
  int m;
  self->norot = (fabs(phi0)<1e-14);
  if (!(self->norot))
    if ((mmax!=self->s_shift-1) || (!FAPPROX(phi0,self->phi0_,1e-12)))
      {
      RESIZE (self->shiftarr,pshtd_cmplx,mmax+1);
      self->s_shift = mmax+1;
      self->phi0_ = phi0;
      for (m=0; m<=mmax; ++m)
        {
        self->shiftarr[m].re = cos(m*phi0);
        self->shiftarr[m].im = sin(m*phi0);
        }
      }
  if (!self->plan) self->plan=make_real_plan(nph);
  if (nph!=(int)self->plan->length)
    {
    kill_real_plan(self->plan);
    self->plan=make_real_plan(nph);
    }
  GROW(self->work,pshtd_cmplx,self->s_work,nph);
  }

static int ringinfo_compare (const void *xa, const void *xb)
  {
  const psht_ringinfo *a = xa, *b=xb;
  if (a->sth < b->sth) return -1;
  if (a->sth > b->sth) return 1;
  return 0;
  }
static int ringpair_compare (const void *xa, const void *xb)
  {
  const psht_ringpair *a = xa, *b=xb;
  if (a->r1.nph==b->r1.nph)
    {
    if (a->r1.phi0<b->r1.phi0) return -1;
    if (a->r1.phi0>b->r1.phi0) return 1;
    return 0;
    }
  if (a->r1.nph<b->r1.nph) return -1;
  if (a->r1.nph>b->r1.nph) return 1;
  return 0;
  }

static void assert_jobspace (int njobs_now)
  {
  UTIL_ASSERT (njobs_now<psht_maxjobs,
    "Too many jobs added to a libpsht joblist. Exiting ...");
  }

void psht_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, psht_alm_info **alm_info)
  {
  int m;
  psht_alm_info *info = RALLOC(psht_alm_info,1);
  info->lmax = lmax;
  info->mmax = mmax;
  info->mstart = RALLOC(ptrdiff_t,mmax+1);
  info->stride = stride;
  for (m=0; m<=mmax; ++m)
    info->mstart[m] = mstart[m];
  *alm_info = info;
  }

ptrdiff_t psht_alm_index (const psht_alm_info *self, int l, int m)
  { return self->mstart[m]+self->stride*l; }

void psht_destroy_alm_info (psht_alm_info *info)
  {
  DEALLOC (info->mstart);
  DEALLOC (info);
  }

void psht_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, psht_geom_info **geom_info)
  {
  psht_geom_info *info = RALLOC(psht_geom_info,1);
  psht_ringinfo *infos = RALLOC(psht_ringinfo,nrings);

  int pos=0;
  int m;
  info->pair=RALLOC(psht_ringpair,nrings);
  info->npairs=0;
  *geom_info = info;

  for (m=0; m<nrings; ++m)
    {
    infos[m].theta = theta[m];
    infos[m].cth = cos(theta[m]);
    infos[m].sth = sin(theta[m]);
    infos[m].weight = weight[m];
    infos[m].phi0 = phi0[m];
    infos[m].ofs = ofs[m];
    infos[m].stride = stride[m];
    infos[m].nph = nph[m];
    }
  qsort(infos,nrings,sizeof(psht_ringinfo),ringinfo_compare);
  while (pos<nrings)
    {
    if ((pos<nrings-1) && FAPPROX(infos[pos].cth,-infos[pos+1].cth,1e-12))
      {
      info->pair[info->npairs].r1=infos[pos];
      info->pair[info->npairs].r2=infos[pos+1];
      pos+=2;
      ++info->npairs;
      }
    else
      {
      info->pair[info->npairs].r1=infos[pos];
      info->pair[info->npairs].r2.nph=-1;
      ++pos;
      ++info->npairs;
      }
    }
  DEALLOC(infos);

  qsort(info->pair,info->npairs,sizeof(psht_ringpair),ringpair_compare);
  }

void psht_destroy_geom_info (psht_geom_info *geom_info)
  {
  DEALLOC (geom_info->pair);
  DEALLOC (geom_info);
  }

#define CONCAT(a,b) a ## b

#define FLT double
#define X(arg) CONCAT(pshtd_,arg)
#include "psht_inc.c"
#undef FLT
#undef X

#define FLT float
#define X(arg) CONCAT(pshts_,arg)
#include "psht_inc.c"
#undef FLT
#undef X

#undef CONCAT

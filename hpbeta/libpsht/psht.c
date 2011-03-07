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
 *  Copyright (C) 2006-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <math.h>
#include "ls_fft.h"
#include "sse_utils.h"
#include "ylmgen_c.h"
#ifdef USE_MPI
#include "psht_mpi.h"
#else
#include "psht.h"
#endif
#include "c_utils.h"

const pshts_cmplx pshts_cmplx_null={0,0};
const pshtd_cmplx pshtd_cmplx_null={0,0};

#define COMPMUL_(a_,b_,c_) \
  { a_.re = b_.re*c_.re - b_.im*c_.im; a_.im = b_.re*c_.im + b_.im*c_.re; }

static void get_chunk_info (int ndata, int *nchunks, int *chunksize)
  {
  static const int chunksize_min=100, nchunks_max=10;
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
  self->norot = (fabs(phi0)<1e-14);
  if (!(self->norot))
    if ((mmax!=self->s_shift-1) || (!FAPPROX(phi0,self->phi0_,1e-12)))
      {
      int m;
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
  return (a->sth<b->sth) ? -1 : (a->sth>b->sth) ? 1 : 0;
  }
static int ringpair_compare (const void *xa, const void *xb)
  {
  const psht_ringpair *a = xa, *b=xb;
  if (a->r1.nph==b->r1.nph)
    return (a->r1.phi0<b->r1.phi0) ? -1 : (a->r1.phi0>b->r1.phi0) ? 1 : 0;
  return (a->r1.nph<b->r1.nph) ? -1 : 1;
  }

static void assert_jobspace (int njobs_now)
  {
  UTIL_ASSERT (njobs_now<psht_maxjobs,
    "Too many jobs added to a libpsht joblist. Exiting ...");
  }

void psht_make_general_alm_info (int lmax, int nm, int stride, const int *mval,
  const ptrdiff_t *mstart, psht_alm_info **alm_info)
  {
  int mi;
  psht_alm_info *info = RALLOC(psht_alm_info,1);
  info->lmax = lmax;
  info->nm = nm;
  info->mval = RALLOC(int,nm);
  info->mvstart = RALLOC(ptrdiff_t,nm);
  info->stride = stride;
  for (mi=0; mi<nm; ++mi)
    {
    info->mval[mi] = mval[mi];
    info->mvstart[mi] = mstart[mi];
    }
  *alm_info = info;
  }

void psht_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, psht_alm_info **alm_info)
  {
  int *mval=RALLOC(int,mmax+1);
  int i;
  for (i=0; i<=mmax; ++i)
    mval[i]=i;
  psht_make_general_alm_info (lmax, mmax+1, stride, mval, mstart, alm_info);
  DEALLOC(mval);
  }

ptrdiff_t psht_alm_index (const psht_alm_info *self, int l, int mi)
  { return self->mvstart[mi]+self->stride*l; }

void psht_destroy_alm_info (psht_alm_info *info)
  {
  DEALLOC (info->mval);
  DEALLOC (info->mvstart);
  DEALLOC (info);
  }

void psht_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, psht_geom_info **geom_info)
  {
  psht_geom_info *info = RALLOC(psht_geom_info,1);
  psht_ringinfo *infos = RALLOC(psht_ringinfo,nrings);

  int pos=0, m;
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
    info->pair[info->npairs].r1=infos[pos];
    if ((pos<nrings-1) && FAPPROX(infos[pos].cth,-infos[pos+1].cth,1e-12))
      {
      info->pair[info->npairs].r2=infos[pos+1];
      ++pos;
      }
    else
      info->pair[info->npairs].r2.nph=-1;
    ++pos;
    ++info->npairs;
    }
  DEALLOC(infos);

  qsort(info->pair,info->npairs,sizeof(psht_ringpair),ringpair_compare);
  }

void psht_destroy_geom_info (psht_geom_info *geom_info)
  {
  DEALLOC (geom_info->pair);
  DEALLOC (geom_info);
  }

static int psht_get_mmax (int *mval, int nm)
  {
  int i;
  int *mcheck=RALLOC(int,nm);
  SET_ARRAY(mcheck,0,nm,0);
  for (i=0; i<nm; ++i)
    {
    int m_cur=mval[i];
    UTIL_ASSERT((m_cur>=0) && (m_cur<nm), "m out of range");
    UTIL_ASSERT(mcheck[m_cur]==0, "duplicate m value");
    mcheck[m_cur]=1;
    }
  DEALLOC(mcheck);
  return nm-1;
  }

#ifdef USE_MPI

typedef struct
  {
  int ntasks;     /* number of tasks */
  int mytask;     /* own task number */
  MPI_Comm comm;  /* communicator to use */

  int *nm;        /* number of m values on every task */
  int *ofs_m;     /* accumulated nm */
  int nmtotal;    /* total number of m values (must be mmax+1) */
  int *mval;      /* array containing all m values of task 0, task 1 etc. */
  int mmax;

  int *npair;     /* number of ring pairs on every task */
  int *ofs_pair;  /* accumulated npair */
  int npairtotal; /* total number of ring pairs */

  double *theta;  /* theta of first ring of every pair on task 0, task 1 etc. */
  int *ispair;    /* is this really a pair? */

  int *almcount, *almdisp, *mapcount, *mapdisp; /* for all2all communication */
  } psht_mpi_info;

static void psht_make_mpi_info (MPI_Comm comm, const psht_alm_info *ainfo,
  const psht_geom_info *ginfo, psht_mpi_info *minfo)
  {
  int i;
  double *theta_tmp;
  int *ispair_tmp;

  minfo->comm = comm;
  MPI_Comm_size (comm, &minfo->ntasks);
  MPI_Comm_rank (comm, &minfo->mytask);

  minfo->nm=RALLOC(int,minfo->ntasks);
  MPI_Allgather ((int *)(&ainfo->nm), 1, MPI_INT, minfo->nm, 1, MPI_INT, comm);
  minfo->ofs_m=RALLOC(int,minfo->ntasks+1);
  minfo->ofs_m[0]=0;
  for (i=1; i<=minfo->ntasks; ++i)
    minfo->ofs_m[i] = minfo->ofs_m[i-1]+minfo->nm[i-1];
  minfo->nmtotal=minfo->ofs_m[minfo->ntasks];
  minfo->mval=RALLOC(int,minfo->nmtotal);
  MPI_Allgatherv(ainfo->mval, ainfo->nm, MPI_INT, minfo->mval, minfo->nm,
    minfo->ofs_m, MPI_INT, comm);

  minfo->mmax=psht_get_mmax(minfo->mval,minfo->nmtotal);

  minfo->npair=RALLOC(int,minfo->ntasks);
  MPI_Allgather ((int *)(&ginfo->npairs), 1, MPI_INT, minfo->npair, 1,
    MPI_INT, comm);
  minfo->ofs_pair=RALLOC(int,minfo->ntasks+1);
  minfo->ofs_pair[0]=0;
  for (i=1; i<=minfo->ntasks; ++i)
    minfo->ofs_pair[i] = minfo->ofs_pair[i-1]+minfo->npair[i-1];
  minfo->npairtotal=minfo->ofs_pair[minfo->ntasks];

  theta_tmp=RALLOC(double,ginfo->npairs);
  ispair_tmp=RALLOC(int,ginfo->npairs);
  for (i=0; i<ginfo->npairs; ++i)
    {
    theta_tmp[i]=ginfo->pair[i].r1.theta;
    ispair_tmp[i]=ginfo->pair[i].r2.nph>0;
    }
  minfo->theta=RALLOC(double,minfo->npairtotal);
  minfo->ispair=RALLOC(int,minfo->npairtotal);
  MPI_Allgatherv(theta_tmp, ginfo->npairs, MPI_DOUBLE, minfo->theta,
    minfo->npair, minfo->ofs_pair, MPI_DOUBLE, comm);
  MPI_Allgatherv(ispair_tmp, ginfo->npairs, MPI_INT, minfo->ispair,
    minfo->npair, minfo->ofs_pair, MPI_INT, comm);
  DEALLOC(theta_tmp);
  DEALLOC(ispair_tmp);

  minfo->almcount=RALLOC(int,minfo->ntasks);
  minfo->almdisp=RALLOC(int,minfo->ntasks+1);
  minfo->mapcount=RALLOC(int,minfo->ntasks);
  minfo->mapdisp=RALLOC(int,minfo->ntasks+1);
  minfo->almdisp[0]=minfo->mapdisp[0]=0;
  for (i=0; i<minfo->ntasks; ++i)
    {
    minfo->almcount[i] = 2*minfo->nm[minfo->mytask]*minfo->npair[i];
    minfo->almdisp[i+1] = minfo->almdisp[i]+minfo->almcount[i];
    minfo->mapcount[i] = 2*minfo->nm[i]*minfo->npair[minfo->mytask];
    minfo->mapdisp[i+1] = minfo->mapdisp[i]+minfo->mapcount[i];
    }
  }

static void psht_destroy_mpi_info (psht_mpi_info *minfo)
  {
  DEALLOC(minfo->nm);
  DEALLOC(minfo->ofs_m);
  DEALLOC(minfo->mval);
  DEALLOC(minfo->npair);
  DEALLOC(minfo->ofs_pair);
  DEALLOC(minfo->theta);
  DEALLOC(minfo->ispair);
  DEALLOC(minfo->almcount);
  DEALLOC(minfo->almdisp);
  DEALLOC(minfo->mapcount);
  DEALLOC(minfo->mapdisp);
  }

static void psht_communicate_alm2map (const psht_mpi_info *minfo,
  pshtd_cmplx **ph)
  {
  int task, th, mi;
  pshtd_cmplx *phas_tmp = RALLOC(pshtd_cmplx,minfo->mapdisp[minfo->ntasks]/2);

  MPI_Alltoallv (*ph,minfo->almcount,minfo->almdisp,MPI_DOUBLE,phas_tmp,
    minfo->mapcount,minfo->mapdisp,MPI_DOUBLE,minfo->comm);

  DEALLOC(*ph);
  ALLOC(*ph,pshtd_cmplx,minfo->npair[minfo->mytask]*minfo->nmtotal);

  for (task=0; task<minfo->ntasks; ++task)
    for (th=0; th<minfo->npair[minfo->mytask]; ++th)
      for (mi=0; mi<minfo->nm[task]; ++mi)
        {
        int m = minfo->mval[mi+minfo->ofs_m[task]];
        (*ph)[th*(minfo->mmax+1) + m] =
          phas_tmp[minfo->mapdisp[task]/2+mi+th*minfo->nm[task]];
        }
  DEALLOC(phas_tmp);
  }

static void psht_communicate_map2alm (const psht_mpi_info *minfo,
  pshtd_cmplx **ph)
  {
  int task, th, mi;
  pshtd_cmplx *phas_tmp = RALLOC(pshtd_cmplx,minfo->mapdisp[minfo->ntasks]/2);

  for (task=0; task<minfo->ntasks; ++task)
    for (th=0; th<minfo->npair[minfo->mytask]; ++th)
      for (mi=0; mi<minfo->nm[task]; ++mi)
        {
        int m = minfo->mval[mi+minfo->ofs_m[task]];
        phas_tmp[minfo->mapdisp[task]/2+mi+th*minfo->nm[task]] =
          (*ph)[th*(minfo->mmax+1) + m];
        }

  DEALLOC(*ph);
  ALLOC(*ph,pshtd_cmplx,minfo->nm[minfo->mytask]*minfo->npairtotal);

  MPI_Alltoallv (phas_tmp,minfo->mapcount,minfo->mapdisp,MPI_DOUBLE,
    *ph,minfo->almcount,minfo->almdisp,MPI_DOUBLE,minfo->comm);

  DEALLOC(phas_tmp);
  }

#endif

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

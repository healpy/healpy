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

/*! \file psht_inc.c
 *  Type-dependent included code for the spherical transform library
 *
 *  Copyright (C) 2006-2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#define COMPMUL_(a_,b_,c_) \
  { a_.re = b_.re*c_.re - b_.im*c_.im; a_.im = b_.re*c_.im + b_.im*c_.re; }

static void X(ringhelper_phase2ring) (ringhelper *self,
  const psht_ringinfo *info, FLT *data, int mmax, const pshtd_cmplx *phase)
  {
  int m;
  int nph = info->nph;
  FLT *ring = data + info->ofs;
  int stride = info->stride;

  ringhelper_update (self, nph, mmax, info->phi0);
  self->work[0]=phase[0];
  SET_ARRAY(self->work,1,nph,pshtd_cmplx_null);

  if (self->norot)
    for (m=1; m<=mmax; ++m)
      {
      int idx1 = m%nph;
      int idx2 = nph-1-((m-1)%nph);
      self->work[idx1].re += phase[m].re; self->work[idx1].im += phase[m].im;
      self->work[idx2].re += phase[m].re; self->work[idx2].im -= phase[m].im;
      }
  else
    for (m=1; m<=mmax; ++m)
      {
      int idx1 = m%nph;
      int idx2 = nph-1-((m-1)%nph);
      pshtd_cmplx tmp;
      COMPMUL_(tmp,phase[m],self->shiftarr[m]);
      self->work[idx1].re += tmp.re; self->work[idx1].im += tmp.im;
      self->work[idx2].re += tmp.re; self->work[idx2].im -= tmp.im;
      }
  real_plan_backward_c (self->plan, &self->work[0].re);
  for (m=0; m<nph; ++m) ring[m*stride] += (FLT)self->work[m].re;
  }

static void X(ringhelper_ring2phase) (ringhelper *self,
  const psht_ringinfo *info, const FLT *data, int mmax, pshtd_cmplx *phase)
  {
  int m;
  int nph = info->nph;
  int maxidx = IMIN(nph-1,mmax);
/* Enable this for traditional Healpix compatibility */
#if 1
  maxidx = mmax;
#endif

  ringhelper_update (self, nph, mmax, -info->phi0);
  for (m=0; m<nph; ++m)
    {
    self->work[m].re = data[info->ofs+m*info->stride]*info->weight;
    self->work[m].im = 0;
    }

  real_plan_forward_c (self->plan, &self->work[0].re);

  if (self->norot)
    for (m=0; m<=maxidx; ++m)
      phase[m] = self->work[m%nph];
  else
    for (m=0; m<=maxidx; ++m)
      COMPMUL_(phase[m],self->work[m%nph],self->shiftarr[m]);

  SET_ARRAY(phase,maxidx+1,mmax+1,pshtd_cmplx_null);
  }

static void X(ringhelper_pair2phase) (ringhelper *self, int mmax,
  const psht_ringpair *pair, const FLT *data,
  pshtd_cmplx *phase1, pshtd_cmplx *phase2)
  {
  if (pair->r1.nph>0)
    X(ringhelper_ring2phase) (self, &(pair->r1), data, mmax, phase1);
  if (pair->r2.nph>0)
    X(ringhelper_ring2phase) (self, &(pair->r2), data, mmax, phase2);
  }

static void X(ringhelper_phase2pair) (ringhelper *self, int mmax,
  const pshtd_cmplx *phase1, const pshtd_cmplx *phase2,
  const psht_ringpair *pair, FLT *data)
  {
  if (pair->r1.nph>0)
    X(ringhelper_phase2ring) (self, &(pair->r1), data, mmax, phase1);
  if (pair->r2.nph>0)
    X(ringhelper_phase2ring) (self, &(pair->r2), data, mmax, phase2);
  }


static void X(fill_map) (const psht_geom_info *ginfo, FLT *map, double value)
  {
  int i,j;
  for (j=0;j<ginfo->npairs;++j)
    {
    for (i=0;i<ginfo->pair[j].r1.nph;++i)
      map[ginfo->pair[j].r1.ofs+i*ginfo->pair[j].r1.stride]=(FLT)value;
    for (i=0;i<ginfo->pair[j].r2.nph;++i)
      map[ginfo->pair[j].r2.ofs+i*ginfo->pair[j].r2.stride]=(FLT)value;
    }
  }

static void X(fill_alm) (const psht_alm_info *ainfo, X(cmplx) *alm,
  X(cmplx) value)
  {
  int l,m;
  for (m=0;m<=ainfo->mmax;++m)
    for (l=m;l<=ainfo->lmax;++l)
      alm[ainfo->mstart[m]+l*ainfo->stride] = value;
  }

static void X(init_output) (X(joblist) *jobs, const psht_geom_info *ginfo,
  const psht_alm_info *alm)
  {
  int ijob,i;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    if (!curjob->add_output)
      switch (curjob->type)
        {
        case MAP2ALM:
          for (i=0; i<curjob->nalm; ++i)
            X(fill_alm) (alm,curjob->alm[i],X(cmplx_null));
          break;
        case ALM2MAP:
        case ALM2MAP_DERIV1:
          for (i=0; i<curjob->nmaps; ++i)
            X(fill_map) (ginfo,curjob->map[i],0);
          break;
        default:
          break;
        }
    }
  }

static void X(alloc_phase) (X(joblist) *jobs, int mmax, int chunksize)
  {
  int ijob,i;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    for (i=0; i<curjob->nmaps; ++i)
      {
      curjob->phas1[i]=RALLOC(pshtd_cmplx,(mmax+1)*chunksize);
      curjob->phas2[i]=RALLOC(pshtd_cmplx,(mmax+1)*chunksize);
      }
    }
  }

static void X(dealloc_phase) (X(joblist) *jobs)
  {
  int ijob,i;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    for (i=0; i<curjob->nmaps; ++i)
      { DEALLOC(curjob->phas1[i]); DEALLOC(curjob->phas2[i]); }
    }
  }

static void X(map2phase) (X(joblist) *jobs, const psht_geom_info *ginfo,
  int mmax, int llim, int ulim)
  {
#pragma omp parallel
{
  ringhelper helper;
  int ith;
  ringhelper_init(&helper);
#pragma omp for schedule(dynamic,1)
  for (ith=llim; ith<ulim; ++ith)
    {
    int ijob,i;
    int dim2 = (ith-llim)*(mmax+1);
    for (ijob=0; ijob<jobs->njobs; ++ijob)
      {
      X(job) *curjob = &jobs->job[ijob];
      switch (curjob->type)
        {
        case MAP2ALM:
          for (i=0; i<curjob->nmaps; ++i)
            X(ringhelper_pair2phase)(&helper,mmax,&ginfo->pair[ith],
              curjob->map[i], &curjob->phas1[i][dim2], &curjob->phas2[i][dim2]);
          break;
        default:
          break;
        }
      }
    }
    ringhelper_destroy(&helper);
} /* end of parallel region */
  }

static void X(alloc_almtmp) (X(joblist) *jobs, int lmax)
  {
  int ijob,i;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    for (i=0; i<curjob->nalm; ++i)
#ifdef PLANCK_HAVE_SSE2
      if (curjob->spin==0)
        curjob->alm_tmp.v[i]=RALLOC(v2df,lmax+1);
      else
        curjob->alm_tmp.v2[i]=RALLOC(v2df2,lmax+1);
#else
      curjob->alm_tmp.c[i]=RALLOC(pshtd_cmplx,lmax+1);
#endif
    }
  }

static void X(dealloc_almtmp) (X(joblist) *jobs)
  {
  int ijob,i;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    for (i=0; i<curjob->nalm; ++i)
#ifdef PLANCK_HAVE_SSE2
      if (curjob->spin==0)
        DEALLOC(curjob->alm_tmp.v[i]);
      else
        DEALLOC(curjob->alm_tmp.v2[i]);
#else
      DEALLOC(curjob->alm_tmp.c[i]);
#endif
    }
  }

static void X(alm2almtmp) (X(joblist) *jobs, int lmax, int m,
  const psht_alm_info *alm)
  {
  int ijob,i,l;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    switch (curjob->type)
      {
      case ALM2MAP:
        {
        const double *norm_l = curjob->norm_l;
        for (i=0; i<curjob->nalm; ++i)
          {
          X(cmplx) *curalm = curjob->alm[i];
          for (l=m; l<=lmax; ++l)
            {
            ptrdiff_t aidx = psht_alm_index(alm,l,m);
#ifdef PLANCK_HAVE_SSE2
            if (curjob->spin==0)
              curjob->alm_tmp.v[i][l] = build_v2df(curalm[aidx].re*norm_l[l],
                                                   curalm[aidx].im*norm_l[l]);
            else
              {
              curjob->alm_tmp.v2[i][l].a = _mm_set1_pd(curalm[aidx].re*norm_l[l]);
              curjob->alm_tmp.v2[i][l].b = _mm_set1_pd(curalm[aidx].im*norm_l[l]);
              }
#else
            curjob->alm_tmp.c[i][l].re = curalm[aidx].re*norm_l[l];
            curjob->alm_tmp.c[i][l].im = curalm[aidx].im*norm_l[l];
#endif
            }
          }
        break;
        }
      case ALM2MAP_DERIV1:
        {
        const double *norm_l = curjob->norm_l;
        for (i=0; i<curjob->nalm; ++i)
          {
          X(cmplx) *curalm = curjob->alm[i];
          for (l=m; l<=lmax; ++l)
            {
            ptrdiff_t aidx = psht_alm_index(alm,l,m);
            double fct = fabs(norm_l[l])*sqrt(l*(l+1.));
#ifdef PLANCK_HAVE_SSE2
            curjob->alm_tmp.v2[i][l].a = _mm_set1_pd(-curalm[aidx].re*fct);
            curjob->alm_tmp.v2[i][l].b = _mm_set1_pd(-curalm[aidx].im*fct);
#else
            curjob->alm_tmp.c[i][l].re = -curalm[aidx].re*fct;
            curjob->alm_tmp.c[i][l].im = -curalm[aidx].im*fct;
#endif
            }
          }
        break;
        }
      case MAP2ALM:
        for (i=0; i<curjob->nalm; ++i)
          {
#ifdef PLANCK_HAVE_SSE2
          if (curjob->spin==0)
            SET_ARRAY(curjob->alm_tmp.v[i],m,lmax+1,_mm_setzero_pd());
          else
            SET_ARRAY(curjob->alm_tmp.v2[i],m,lmax+1,zero_v2df2());
#else
          SET_ARRAY(curjob->alm_tmp.c[i],m,lmax+1,pshtd_cmplx_null);
#endif
          }
        break;
      default:
        break;
      }
    }
  }


#ifdef PLANCK_HAVE_SSE2

#define ALM2MAP_MACRO(px) \
  { \
  v2df ar=_mm_shuffle_pd(almtmp[l],almtmp[l],_MM_SHUFFLE2(0,0)), \
       ai=_mm_shuffle_pd(almtmp[l],almtmp[l],_MM_SHUFFLE2(1,1)); \
  px.a=_mm_add_pd(px.a,_mm_mul_pd(ar,Ylm[l])); \
  px.b=_mm_add_pd(px.b,_mm_mul_pd(ai,Ylm[l])); \
  ++l; \
  }

#define ALM2MAP_SPIN_MACRO(Qx,Qy,Ux,Uy) \
  { \
  const v2df lw = lwx[l].a, lx = lwx[l].b; \
  Qx.a=_mm_add_pd(Qx.a,_mm_mul_pd(almtmpG[l].a,lw)); \
  Uy.b=_mm_sub_pd(Uy.b,_mm_mul_pd(almtmpG[l].a,lx)); \
  Qx.b=_mm_add_pd(Qx.b,_mm_mul_pd(almtmpG[l].b,lw)); \
  Uy.a=_mm_add_pd(Uy.a,_mm_mul_pd(almtmpG[l].b,lx)); \
  Ux.a=_mm_add_pd(Ux.a,_mm_mul_pd(almtmpC[l].a,lw)); \
  Qy.b=_mm_add_pd(Qy.b,_mm_mul_pd(almtmpC[l].a,lx)); \
  Ux.b=_mm_add_pd(Ux.b,_mm_mul_pd(almtmpC[l].b,lw)); \
  Qy.a=_mm_sub_pd(Qy.a,_mm_mul_pd(almtmpC[l].b,lx)); \
  ++l; \
  }

#define ALM2MAP_DERIV1_MACRO(Qx,Uy) \
  { \
  const v2df lw = lwx[l].a, lx = lwx[l].b; \
  Qx.a = _mm_add_pd(Qx.a,_mm_mul_pd(almtmp[l].a,lw)); \
  Uy.b = _mm_sub_pd(Uy.b,_mm_mul_pd(almtmp[l].a,lx)); \
  Qx.b = _mm_add_pd(Qx.b,_mm_mul_pd(almtmp[l].b,lw)); \
  Uy.a = _mm_add_pd(Uy.a,_mm_mul_pd(almtmp[l].b,lx)); \
  ++l; \
  }

#ifdef PLANCK_HAVE_SSE3
#define MAP2ALM_MACRO(px) \
  { \
  const v2df t1_ = _mm_mul_pd(px.a,Ylm[l]), t2_ = _mm_mul_pd(px.b,Ylm[l]); \
  const v2df t5_=_mm_hadd_pd(t1_,t2_);\
  almtmp[l] = _mm_add_pd(almtmp[l],t5_); \
  ++l; \
  }
#else
#define MAP2ALM_MACRO(px) \
  { \
  const v2df t1_ = _mm_mul_pd(px.a,Ylm[l]), t2_ = _mm_mul_pd(px.b,Ylm[l]); \
  const v2df t3_ = _mm_shuffle_pd(t1_,t2_,_MM_SHUFFLE2(0,0)), \
             t4_ = _mm_shuffle_pd(t1_,t2_,_MM_SHUFFLE2(1,1)); \
  const v2df t5_=_mm_add_pd(t3_,t4_);\
  almtmp[l] = _mm_add_pd(almtmp[l],t5_); \
  ++l; \
  }
#endif

#define MAP2ALM_SPIN_MACRO(Qx,Qy,Ux,Uy) \
  { \
  const v2df lw = lwx[l].a, lx = lwx[l].b; \
  almtmpG[l].a=_mm_add_pd(almtmpG[l].a, \
    _mm_sub_pd(_mm_mul_pd(Qx.a,lw),_mm_mul_pd(Uy.b,lx))); \
  almtmpG[l].b=_mm_add_pd(almtmpG[l].b, \
    _mm_add_pd(_mm_mul_pd(Qx.b,lw),_mm_mul_pd(Uy.a,lx))); \
  almtmpC[l].a=_mm_add_pd(almtmpC[l].a, \
    _mm_add_pd(_mm_mul_pd(Ux.a,lw),_mm_mul_pd(Qy.b,lx))); \
  almtmpC[l].b=_mm_add_pd(almtmpC[l].b, \
    _mm_sub_pd(_mm_mul_pd(Ux.b,lw),_mm_mul_pd(Qy.a,lx))); \
  ++l; \
  }

#define EXTRACT(pa,pb,pha,phb,phc,phd) \
  { \
  read_v2df (_mm_add_pd(pa.a,pb.a), &pha->re, &phc->re); \
  read_v2df (_mm_sub_pd(pa.a,pb.a), &phb->re, &phd->re); \
  read_v2df (_mm_add_pd(pa.b,pb.b), &pha->im, &phc->im); \
  read_v2df (_mm_sub_pd(pa.b,pb.b), &phb->im, &phd->im); \
  }

#define COMPOSE(pa,pb,pha,phb,phc,phd) \
  { \
  pa.a = build_v2df(pha.re+phb.re,phc.re+phd.re); \
  pa.b = build_v2df(pha.im+phb.im,phc.im+phd.im); \
  pb.a = build_v2df(pha.re-phb.re,phc.re-phd.re); \
  pb.b = build_v2df(pha.im-phb.im,phc.im-phd.im); \
  }

static void X(inner_loop) (X(joblist) *jobs, const psht_geom_info *ginfo,
  int lmax, int mmax, int llim, int ulim, Ylmgen_C *generator, int m)
  {
  const v2df2 v2df2_zero = zero_v2df2();
  int ith,ijob;
  for (ith=0; ith<ulim-llim; ith+=2)
    {
    pshtd_cmplx dum;
    int dual = (ith+1)<(ulim-llim),
        rpair1 = ginfo->pair[ith+llim].r2.nph>0,
        rpair2 = dual && (ginfo->pair[ith+1+llim].r2.nph>0);
    int phas_idx1 =  ith   *(mmax+1)+m,
        phas_idx2 = (ith+1)*(mmax+1)+m;
    Ylmgen_prepare_sse2(generator,ith,dual?(ith+1):ith,m);

    for (ijob=0; ijob<jobs->njobs; ++ijob)
      {
      X(job) *curjob = &jobs->job[ijob];
      switch (curjob->type)
        {
        case ALM2MAP:
          {
          switch (curjob->spin)
            {
            case 0:
              {
              pshtd_cmplx *ph1 =          &curjob->phas1[0][phas_idx1],
                          *ph2 = rpair1 ? &curjob->phas2[0][phas_idx1] : &dum,
                          *ph3 = dual   ? &curjob->phas1[0][phas_idx2] : &dum,
                          *ph4 = rpair2 ? &curjob->phas2[0][phas_idx2] : &dum;
              *ph1 = *ph2 = *ph3 = *ph4 = pshtd_cmplx_null;
              Ylmgen_recalc_Ylm_sse2 (generator);
              if (generator->firstl[0]<=lmax)
                {
                int l = generator->firstl[0];
                v2df2 p1=v2df2_zero, p2=v2df2_zero;
                const v2df *Ylm = generator->ylm_sse2;
                const v2df *almtmp = curjob->alm_tmp.v[0];

                if ((l-m)&1)
                  ALM2MAP_MACRO(p2)
                for (;l<lmax;)
                  {
                  ALM2MAP_MACRO(p1)
                  ALM2MAP_MACRO(p2)
                  }
                if (l==lmax)
                  ALM2MAP_MACRO(p1)

                EXTRACT (p1,p2,ph1,ph2,ph3,ph4)
                }
              break;
              }

            default:
              {
              pshtd_cmplx *ph1Q =          &curjob->phas1[0][phas_idx1],
                          *ph1U =          &curjob->phas1[1][phas_idx1],
                          *ph2Q = rpair1 ? &curjob->phas2[0][phas_idx1] : &dum,
                          *ph2U = rpair1 ? &curjob->phas2[1][phas_idx1] : &dum,
                          *ph3Q = dual   ? &curjob->phas1[0][phas_idx2] : &dum,
                          *ph3U = dual   ? &curjob->phas1[1][phas_idx2] : &dum,
                          *ph4Q = rpair2 ? &curjob->phas2[0][phas_idx2] : &dum,
                          *ph4U = rpair2 ? &curjob->phas2[1][phas_idx2] : &dum;
              *ph1Q = *ph2Q = *ph1U = *ph2U = *ph3Q = *ph4Q = *ph3U = *ph4U
                = pshtd_cmplx_null;
              Ylmgen_recalc_lambda_wx_sse2 (generator,curjob->spin);
              if (generator->firstl[curjob->spin]<=lmax)
                {
                int l = generator->firstl[curjob->spin];
                v2df2 p1Q=v2df2_zero, p2Q=v2df2_zero,
                      p1U=v2df2_zero, p2U=v2df2_zero;
                const v2df2 *lwx = generator->lambda_wx_sse2[curjob->spin];
                const v2df2 *almtmpG = curjob->alm_tmp.v2[0],
                            *almtmpC = curjob->alm_tmp.v2[1];

                if ((l-m+curjob->spin)&1)
                  ALM2MAP_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                for (;l<lmax;)
                  {
                  ALM2MAP_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                  ALM2MAP_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                  }
                if (l==lmax)
                  ALM2MAP_SPIN_MACRO(p1Q,p2Q,p1U,p2U)

                EXTRACT (p1Q,p2Q,ph1Q,ph2Q,ph3Q,ph4Q)
                EXTRACT (p1U,p2U,ph1U,ph2U,ph3U,ph4U)
                }
              break;
              }
            }
          break;
          }
        case ALM2MAP_DERIV1:
          {
          pshtd_cmplx *ph1Q =          &curjob->phas1[0][phas_idx1],
                      *ph1U =          &curjob->phas1[1][phas_idx1],
                      *ph2Q = rpair1 ? &curjob->phas2[0][phas_idx1] : &dum,
                      *ph2U = rpair1 ? &curjob->phas2[1][phas_idx1] : &dum,
                      *ph3Q = dual   ? &curjob->phas1[0][phas_idx2] : &dum,
                      *ph3U = dual   ? &curjob->phas1[1][phas_idx2] : &dum,
                      *ph4Q = rpair2 ? &curjob->phas2[0][phas_idx2] : &dum,
                      *ph4U = rpair2 ? &curjob->phas2[1][phas_idx2] : &dum;
          *ph1Q = *ph2Q = *ph1U = *ph2U = *ph3Q = *ph4Q = *ph3U = *ph4U
            = pshtd_cmplx_null;
          Ylmgen_recalc_lambda_wx_sse2 (generator,1);
          if (generator->firstl[1]<=lmax)
            {
            int l = generator->firstl[1];
            v2df2 p1Q=v2df2_zero, p2Q=v2df2_zero,
                  p1U=v2df2_zero, p2U=v2df2_zero;
            const v2df2 *lwx = generator->lambda_wx_sse2[1];
            const v2df2 *almtmp = curjob->alm_tmp.v2[0];

            if ((l-m+1)&1)
              ALM2MAP_DERIV1_MACRO(p2Q,p1U)
            for (;l<lmax;)
              {
              ALM2MAP_DERIV1_MACRO(p1Q,p2U)
              ALM2MAP_DERIV1_MACRO(p2Q,p1U)
              }
            if (l==lmax)
              ALM2MAP_DERIV1_MACRO(p1Q,p2U)

            EXTRACT (p1Q,p2Q,ph1Q,ph2Q,ph3Q,ph4Q)
            EXTRACT (p1U,p2U,ph1U,ph2U,ph3U,ph4U)
            }
          break;
          }
        case MAP2ALM:
          {
          switch (curjob->spin)
            {
            case 0:
              {
              Ylmgen_recalc_Ylm_sse2 (generator);
              if (generator->firstl[0]<=lmax)
                {
                int l = generator->firstl[0];
                const v2df *Ylm = generator->ylm_sse2;
                v2df *almtmp = curjob->alm_tmp.v[0];
                pshtd_cmplx
                  ph1 =          curjob->phas1[0][phas_idx1],
                  ph2 = rpair1 ? curjob->phas2[0][phas_idx1] : pshtd_cmplx_null,
                  ph3 = dual   ? curjob->phas1[0][phas_idx2] : pshtd_cmplx_null,
                  ph4 = rpair2 ? curjob->phas2[0][phas_idx2] : pshtd_cmplx_null;
                v2df2 p1,p2;
                COMPOSE (p1,p2,ph1,ph2,ph3,ph4)

                if ((l-m)&1)
                  MAP2ALM_MACRO(p2)
                for (;l<lmax;)
                  {
                  MAP2ALM_MACRO(p1)
                  MAP2ALM_MACRO(p2)
                  }
                if (l==lmax)
                  MAP2ALM_MACRO(p1)
                }
              break;
              }

            default:
              {
              Ylmgen_recalc_lambda_wx_sse2 (generator,curjob->spin);
              if (generator->firstl[curjob->spin]<=lmax)
                {
                int l = generator->firstl[curjob->spin];
                const v2df2 *lwx = generator->lambda_wx_sse2[curjob->spin];
                v2df2 *almtmpG = curjob->alm_tmp.v2[0],
                      *almtmpC = curjob->alm_tmp.v2[1];
                pshtd_cmplx
                  ph1Q =          curjob->phas1[0][phas_idx1],
                  ph1U =          curjob->phas1[1][phas_idx1],
                  ph2Q = rpair1 ? curjob->phas2[0][phas_idx1] :pshtd_cmplx_null,
                  ph2U = rpair1 ? curjob->phas2[1][phas_idx1] :pshtd_cmplx_null,
                  ph3Q = dual   ? curjob->phas1[0][phas_idx2] :pshtd_cmplx_null,
                  ph3U = dual   ? curjob->phas1[1][phas_idx2] :pshtd_cmplx_null,
                  ph4Q = rpair2 ? curjob->phas2[0][phas_idx2] :pshtd_cmplx_null,
                  ph4U = rpair2 ? curjob->phas2[1][phas_idx2] :pshtd_cmplx_null;
                v2df2 p1Q, p2Q, p1U, p2U;
                COMPOSE (p1Q,p2Q,ph1Q,ph2Q,ph3Q,ph4Q)
                COMPOSE (p1U,p2U,ph1U,ph2U,ph3U,ph4U)

                if ((l-m+curjob->spin)&1)
                  MAP2ALM_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                for (;l<lmax;)
                  {
                  MAP2ALM_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                  MAP2ALM_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                  }
                if (l==lmax)
                  MAP2ALM_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                }
              break;
              }
            }
          break;
          }
        }
      }
    }
  }

#else

#define ALM2MAP_MACRO(px) \
  { \
  px.re += almtmp[l].re*Ylm[l]; px.im += almtmp[l].im*Ylm[l]; \
  ++l; \
  }

#define MAP2ALM_MACRO(px) \
  { \
  almtmp[l].re += px.re*Ylm[l]; almtmp[l].im += px.im*Ylm[l]; \
  ++l; \
  }

#define ALM2MAP_SPIN_MACRO(Qx,Qy,Ux,Uy) \
  { \
  const double lw = lwx[l][0], lx = lwx[l][1]; \
  Qx.re+=almtmpG[l].re*lw; Qx.im+=almtmpG[l].im*lw; \
  Ux.re+=almtmpC[l].re*lw; Ux.im+=almtmpC[l].im*lw; \
  Qy.re-=almtmpC[l].im*lx; Qy.im+=almtmpC[l].re*lx; \
  Uy.re+=almtmpG[l].im*lx; Uy.im-=almtmpG[l].re*lx; \
  ++l; \
  }

#define ALM2MAP_DERIV1_MACRO(Qx,Uy) \
  { \
  const double lw = lwx[l][0], lx = lwx[l][1]; \
  Qx.re+=almtmp[l].re*lw; Qx.im+=almtmp[l].im*lw; \
  Uy.re+=almtmp[l].im*lx; Uy.im-=almtmp[l].re*lx; \
  ++l; \
  }

#define MAP2ALM_SPIN_MACRO(Qx,Qy,Ux,Uy) \
  { \
  const double lw = lwx[l][0], lx = lwx[l][1]; \
  almtmpG[l].re += Qx.re*lw - Uy.im*lx; almtmpG[l].im += Qx.im*lw + Uy.re*lx; \
  almtmpC[l].re += Ux.re*lw + Qy.im*lx; almtmpC[l].im += Ux.im*lw - Qy.re*lx; \
  ++l; \
  }

static void X(inner_loop) (X(joblist) *jobs, const psht_geom_info *ginfo,
  int lmax, int mmax, int llim, int ulim, Ylmgen_C *generator, int m)
  {
  int ith,ijob;
  for (ith=0; ith<ulim-llim; ++ith)
    {
    pshtd_cmplx dum;
    int rpair = ginfo->pair[ith+llim].r2.nph>0;
    int phas_idx = ith*(mmax+1)+m;
    Ylmgen_prepare(generator,ith,m);

    for (ijob=0; ijob<jobs->njobs; ++ijob)
      {
      X(job) *curjob = &jobs->job[ijob];
      switch (curjob->type)
        {
        case ALM2MAP:
          {
          switch (curjob->spin)
            {
            case 0:
              {
              pshtd_cmplx *ph1 =         &curjob->phas1[0][phas_idx],
                          *ph2 = rpair ? &curjob->phas2[0][phas_idx] : &dum;
              *ph1 = *ph2 = pshtd_cmplx_null;
              Ylmgen_recalc_Ylm (generator);
              if (generator->firstl[0]<=lmax)
                {
                int l = generator->firstl[0];
                pshtd_cmplx p1=pshtd_cmplx_null,p2=pshtd_cmplx_null;
                const double *Ylm = generator->ylm;
                const pshtd_cmplx *almtmp = curjob->alm_tmp.c[0];

                if ((l-m)&1)
                  ALM2MAP_MACRO(p2)
                for (;l<lmax;)
                  {
                  ALM2MAP_MACRO(p1)
                  ALM2MAP_MACRO(p2)
                  }
                if (l==lmax)
                  ALM2MAP_MACRO(p1)

                ph1->re = p1.re+p2.re; ph1->im = p1.im+p2.im;
                ph2->re = p1.re-p2.re; ph2->im = p1.im-p2.im;
                }
              break;
              }

            default:
              {
              pshtd_cmplx *ph1Q =         &curjob->phas1[0][phas_idx],
                          *ph1U =         &curjob->phas1[1][phas_idx],
                          *ph2Q = rpair ? &curjob->phas2[0][phas_idx] : &dum,
                          *ph2U = rpair ? &curjob->phas2[1][phas_idx] : &dum;
              *ph1Q = *ph2Q = *ph1U = *ph2U = pshtd_cmplx_null;
              Ylmgen_recalc_lambda_wx (generator,curjob->spin);
              if (generator->firstl[curjob->spin]<=lmax)
                {
                int l = generator->firstl[curjob->spin];
                pshtd_cmplx p1Q=pshtd_cmplx_null,p2Q=pshtd_cmplx_null,
                            p1U=pshtd_cmplx_null,p2U=pshtd_cmplx_null;
                ylmgen_dbl2 *lwx = generator->lambda_wx[curjob->spin];
                const pshtd_cmplx *almtmpG = curjob->alm_tmp.c[0],
                                  *almtmpC = curjob->alm_tmp.c[1];

                if ((l-m+curjob->spin)&1)
                  ALM2MAP_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                for (;l<lmax;)
                  {
                  ALM2MAP_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                  ALM2MAP_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                  }
                if (l==lmax)
                  ALM2MAP_SPIN_MACRO(p1Q,p2Q,p1U,p2U)

                ph1Q->re = p1Q.re+p2Q.re; ph1Q->im = p1Q.im+p2Q.im;
                ph1U->re = p1U.re+p2U.re; ph1U->im = p1U.im+p2U.im;
                ph2Q->re = p1Q.re-p2Q.re; ph2Q->im = p1Q.im-p2Q.im;
                ph2U->re = p1U.re-p2U.re; ph2U->im = p1U.im-p2U.im;
                }
              break;
              }
            }
          break;
          }
        case ALM2MAP_DERIV1:
          {
          pshtd_cmplx *ph1Q =         &curjob->phas1[0][phas_idx],
                      *ph1U =         &curjob->phas1[1][phas_idx],
                      *ph2Q = rpair ? &curjob->phas2[0][phas_idx] : &dum,
                      *ph2U = rpair ? &curjob->phas2[1][phas_idx] : &dum;
          *ph1Q = *ph2Q = *ph1U = *ph2U = pshtd_cmplx_null;
          Ylmgen_recalc_lambda_wx (generator,1);
          if (generator->firstl[1]<=lmax)
            {
            int l = generator->firstl[1];
            pshtd_cmplx p1Q=pshtd_cmplx_null,p2Q=pshtd_cmplx_null,
                        p1U=pshtd_cmplx_null,p2U=pshtd_cmplx_null;
            ylmgen_dbl2 *lwx = generator->lambda_wx[1];
            const pshtd_cmplx *almtmp = curjob->alm_tmp.c[0];

            if ((l-m+1)&1)
              ALM2MAP_DERIV1_MACRO(p2Q,p1U)
            for (;l<lmax;)
              {
              ALM2MAP_DERIV1_MACRO(p1Q,p2U)
              ALM2MAP_DERIV1_MACRO(p2Q,p1U)
              }
            if (l==lmax)
              ALM2MAP_DERIV1_MACRO(p1Q,p2U)

            ph1Q->re = p1Q.re+p2Q.re; ph1Q->im = p1Q.im+p2Q.im;
            ph1U->re = p1U.re+p2U.re; ph1U->im = p1U.im+p2U.im;
            ph2Q->re = p1Q.re-p2Q.re; ph2Q->im = p1Q.im-p2Q.im;
            ph2U->re = p1U.re-p2U.re; ph2U->im = p1U.im-p2U.im;
            }
          break;
          }
        case MAP2ALM:
          {
          switch (curjob->spin)
            {
            case 0:
              {
              Ylmgen_recalc_Ylm (generator);
              if (generator->firstl[0]<=lmax)
                {
                int l = generator->firstl[0];
                const double *Ylm = generator->ylm;
                pshtd_cmplx *almtmp = curjob->alm_tmp.c[0];
                const pshtd_cmplx
                  ph1 =         curjob->phas1[0][phas_idx],
                  ph2 = rpair ? curjob->phas2[0][phas_idx] : pshtd_cmplx_null;
                pshtd_cmplx p1,p2;
                p1.re = ph1.re+ph2.re; p1.im = ph1.im+ph2.im;
                p2.re = ph1.re-ph2.re; p2.im = ph1.im-ph2.im;

                if ((l-m)&1)
                  MAP2ALM_MACRO(p2)
                for (;l<lmax;)
                  {
                  MAP2ALM_MACRO(p1)
                  MAP2ALM_MACRO(p2)
                  }
                if (l==lmax)
                  MAP2ALM_MACRO(p1)
                }
              break;
              }

            default:
              {
              Ylmgen_recalc_lambda_wx (generator,curjob->spin);
              if (generator->firstl[curjob->spin]<=lmax)
                {
                int l = generator->firstl[curjob->spin];
                ylmgen_dbl2 *lwx = generator->lambda_wx[curjob->spin];
                pshtd_cmplx *almtmpG = curjob->alm_tmp.c[0],
                            *almtmpC = curjob->alm_tmp.c[1];
                const pshtd_cmplx
                  ph1Q =         curjob->phas1[0][phas_idx],
                  ph1U =         curjob->phas1[1][phas_idx],
                  ph2Q = rpair ? curjob->phas2[0][phas_idx] : pshtd_cmplx_null,
                  ph2U = rpair ? curjob->phas2[1][phas_idx] : pshtd_cmplx_null;
                pshtd_cmplx p1Q,p2Q,p1U,p2U;
                p1Q.re = ph1Q.re+ph2Q.re; p1Q.im = ph1Q.im+ph2Q.im;
                p2Q.re = ph1Q.re-ph2Q.re; p2Q.im = ph1Q.im-ph2Q.im;
                p1U.re = ph1U.re+ph2U.re; p1U.im = ph1U.im+ph2U.im;
                p2U.re = ph1U.re-ph2U.re; p2U.im = ph1U.im-ph2U.im;

                if ((l-m+curjob->spin)&1)
                  MAP2ALM_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                for (;l<lmax;)
                  {
                  MAP2ALM_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                  MAP2ALM_SPIN_MACRO(p2Q,p1Q,p2U,p1U)
                  }
                if (l==lmax)
                  MAP2ALM_SPIN_MACRO(p1Q,p2Q,p1U,p2U)
                }
              break;
              }
            }
          break;
          }
        }
      }
    }
  }

#endif

static void X(almtmp2alm) (X(joblist) *jobs, int lmax, int m,
  const psht_alm_info *alm)
  {
  int ijob,i,l;
  for (ijob=0; ijob<jobs->njobs; ++ijob)
    {
    X(job) *curjob = &jobs->job[ijob];
    switch (curjob->type)
      {
      case MAP2ALM:
        {
        const double *norm_l = curjob->norm_l;
        for (i=0;i<curjob->nalm;++i)
          {
          X(cmplx) *curalm = curjob->alm[i];
          for (l=m; l<=lmax; ++l)
            {
            ptrdiff_t aidx = psht_alm_index(alm,l,m);
#ifdef PLANCK_HAVE_SSE2
            if (curjob->spin==0)
              {
              V2DF t; t.v = curjob->alm_tmp.v[i][l];
              curalm[aidx].re += (FLT)(t.d[0]*norm_l[l]);
              curalm[aidx].im += (FLT)(t.d[1]*norm_l[l]);
              }
            else
              {
              V2DF2 t = to_V2DF2(curjob->alm_tmp.v2[i][l]);
              curalm[aidx].re += (FLT)((t.a.d[0]+t.a.d[1])*norm_l[l]);
              curalm[aidx].im += (FLT)((t.b.d[0]+t.b.d[1])*norm_l[l]);
              }
#else
            curalm[aidx].re += (FLT)(curjob->alm_tmp.c[i][l].re*norm_l[l]);
            curalm[aidx].im += (FLT)(curjob->alm_tmp.c[i][l].im*norm_l[l]);
#endif
            }
          }
        break;
        }
      default:
        break;
      }
    }
  }

static void X(phase2map) (X(joblist) *jobs, const psht_geom_info *ginfo,
  int mmax, int llim, int ulim)
  {
#pragma omp parallel
{
  ringhelper helper;
  int ith;
  ringhelper_init(&helper);
#pragma omp for schedule(dynamic,1)
  for (ith=llim; ith<ulim; ++ith)
    {
    int ijob,i;
    int dim2 = (ith-llim)*(mmax+1);
    for (ijob=0; ijob<jobs->njobs; ++ijob)
      {
      X(job) *curjob = &jobs->job[ijob];
      switch (curjob->type)
        {
        case ALM2MAP:
        case ALM2MAP_DERIV1:
          for (i=0; i<curjob->nmaps; ++i)
            X(ringhelper_phase2pair)(&helper,mmax,&curjob->phas1[i][dim2],
              &curjob->phas2[i][dim2],&ginfo->pair[ith],curjob->map[i]);
          break;
        default:
          break;
        }
      }
    }
  ringhelper_destroy(&helper);
} /* end of parallel region */
  }

void X(execute_jobs) (X(joblist) *joblist, const psht_geom_info *geom_info,
  const psht_alm_info *alm_info)
  {
  int lmax = alm_info->lmax, mmax = alm_info->mmax;
  int nchunks, chunksize, chunk, spinrec=0, ijob;

  for (ijob=0; ijob<joblist->njobs; ++ijob)
    if (joblist->job[ijob].spin<=1) { spinrec=1; break; }

  for (ijob=0; ijob<joblist->njobs; ++ijob)
    joblist->job[ijob].norm_l =
      Ylmgen_get_norm (lmax, joblist->job[ijob].spin, spinrec);

/* clear output arrays if requested */
  X(init_output) (joblist, geom_info, alm_info);

  get_chunk_info(geom_info->npairs,&nchunks,&chunksize);
  X(alloc_phase) (joblist,mmax,chunksize);

/* chunk loop */
  for (chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=IMIN(llim+chunksize,geom_info->npairs);

/* map->phase where necessary */
    X(map2phase) (joblist, geom_info, mmax, llim, ulim);

#pragma omp parallel
{
    int m;
    X(joblist) ljobs = *joblist;
    Ylmgen_C generator;
    double *theta = RALLOC(double,ulim-llim);
    for (m=0; m<ulim-llim; ++m)
      theta[m] = geom_info->pair[m+llim].r1.theta;
    Ylmgen_init (&generator,lmax,mmax,spinrec,1e-30);
    Ylmgen_set_theta (&generator,theta,ulim-llim);
    DEALLOC(theta);
    X(alloc_almtmp)(&ljobs,lmax);

#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
/* alm->alm_tmp where necessary */
      X(alm2almtmp) (&ljobs, lmax, m, alm_info);

/* inner conversion loop */
      X(inner_loop) (&ljobs, geom_info, lmax, mmax, llim, ulim, &generator, m);

/* alm_tmp->alm where necessary */
      X(almtmp2alm) (&ljobs, lmax, m, alm_info);
      }

    Ylmgen_destroy(&generator);
    X(dealloc_almtmp)(&ljobs);
} /* end of parallel region */

/* phase->map where necessary */
    X(phase2map) (joblist, geom_info, mmax, llim, ulim);
    } /* end of chunk loop */

  for (ijob=0; ijob<joblist->njobs; ++ijob)
    DEALLOC(joblist->job[ijob].norm_l);
  X(dealloc_phase) (joblist);
  }

void X(make_joblist) (X(joblist) **joblist)
  {
  *joblist = RALLOC(X(joblist),1);
  (*joblist)->njobs=0;
  }

void X(clear_joblist) (X(joblist) *joblist)
  { joblist->njobs=0; }

void X(destroy_joblist) (X(joblist) *joblist)
  { DEALLOC(joblist); }

static void X(addjob) (X(joblist) *joblist, psht_jobtype type, int spin,
  int add_output, int nalm, int nmaps, X(cmplx) *alm0, X(cmplx) *alm1,
  X(cmplx) *alm2, FLT *map0, FLT *map1, FLT *map2)
  {
  assert_jobspace(joblist->njobs);
  joblist->job[joblist->njobs].type = type;
  joblist->job[joblist->njobs].spin = spin;
  joblist->job[joblist->njobs].norm_l = NULL;
  joblist->job[joblist->njobs].alm[0] = alm0;
  joblist->job[joblist->njobs].alm[1] = alm1;
  joblist->job[joblist->njobs].alm[2] = alm2;
  joblist->job[joblist->njobs].map[0] = map0;
  joblist->job[joblist->njobs].map[1] = map1;
  joblist->job[joblist->njobs].map[2] = map2;
  joblist->job[joblist->njobs].add_output = add_output;
  joblist->job[joblist->njobs].nmaps = nmaps;
  joblist->job[joblist->njobs].nalm = nalm;
  ++joblist->njobs;
  }

void X(add_job_alm2map) (X(joblist) *joblist, const X(cmplx) *alm, FLT *map,
  int add_output)
  {
  X(addjob) (joblist, ALM2MAP, 0, add_output, 1, 1,
    (X(cmplx) *)alm, NULL, NULL, map, NULL, NULL);
  }
void X(add_job_map2alm) (X(joblist) *joblist, const FLT *map, X(cmplx) *alm,
  int add_output)
  {
  X(addjob) (joblist, MAP2ALM, 0, add_output, 1, 1,
    alm, NULL, NULL, (FLT *)map, NULL, NULL);
  }
void X(add_job_alm2map_spin) (X(joblist) *joblist, const X(cmplx) *alm1,
  const X(cmplx) *alm2, FLT *map1, FLT *map2, int spin, int add_output)
  {
  UTIL_ASSERT((spin>0)&&(spin<=Ylmgen_maxspin()),
    "bad spin in add_job_alm2map_spin()");
  X(addjob) (joblist, ALM2MAP, spin, add_output, 2, 2,
    (X(cmplx) *)alm1, (X(cmplx) *)alm2, NULL, map1, map2, NULL);
  }
void X(add_job_alm2map_pol) (X(joblist) *joblist, const X(cmplx) *almT,
  const X(cmplx) *almG, const X(cmplx) *almC, FLT *mapT, FLT *mapQ, FLT *mapU,
  int add_output)
  {
  X(add_job_alm2map) (joblist, almT, mapT, add_output);
  X(add_job_alm2map_spin) (joblist, almG, almC, mapQ, mapU, 2, add_output);
  }
void X(add_job_map2alm_spin) (X(joblist) *joblist, const FLT *map1,
  const FLT *map2, X(cmplx) *alm1, X(cmplx) *alm2, int spin, int add_output)
  {
  UTIL_ASSERT((spin>0)&&(spin<=Ylmgen_maxspin()),
    "bad spin in add_job_map2alm_spin()");
  X(addjob) (joblist, MAP2ALM, spin, add_output, 2, 2,
    alm1, alm2, NULL, (FLT *)map1, (FLT *)map2, NULL);
  }
void X(add_job_map2alm_pol) (X(joblist) *joblist, const FLT *mapT,
  const FLT *mapQ, const FLT *mapU, X(cmplx) *almT, X(cmplx) *almG,
  X(cmplx) *almC, int add_output)
  {
  X(add_job_map2alm) (joblist, mapT, almT, add_output);
  X(add_job_map2alm_spin) (joblist, mapQ, mapU, almG, almC, 2, add_output);
  }
void X(add_job_alm2map_deriv1) (X(joblist) *joblist, const X(cmplx) *alm,
  FLT *mapdtheta, FLT *mapdphi, int add_output)
  {
  X(addjob) (joblist, ALM2MAP_DERIV1, 1, add_output, 1, 2,
    (X(cmplx) *)alm, NULL, NULL, mapdtheta, mapdphi, NULL);
  }

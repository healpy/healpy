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

/*! \file sharp_core_inc2.c
 *  Type-dependent code for the computational core
 *
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

static void Z(alm2map_kernel) (const Tb cth, Y(Tbri) * restrict p1,
  Y(Tbri) * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, const dcmplx * restrict alm,
  int l, int lmax NJ1)
  {
if (njobs>1)
  {
  while (l<lmax-2)
    {
    Tb lam_3, lam_4;
    Tv r0=vload(rf[l].f[0]),r1=vload(rf[l].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_3.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    r0=vload(rf[l+1].f[0]);r1=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_4.v[i] = vsub(vmul(vmul(cth.v[i],lam_3.v[i]),r0),vmul(lam_2.v[i],r1));
    r0=vload(rf[l+2].f[0]);r1=vload(rf[l+2].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_4.v[i]),r0),vmul(lam_3.v[i],r1));
    for (int j=0; j<njobs; ++j)
      {
      Tv ar2=vload(creal(alm[njobs*l+j])),
         ai2=vload(cimag(alm[njobs*l+j])),
         ar4=vload(creal(alm[njobs*(l+2)+j])),
         ai4=vload(cimag(alm[njobs*(l+2)+j]));
      for (int i=0; i<nvec; ++i)
        {
        vfmaaeq(p1[j].r.v[i],lam_2.v[i],ar2,lam_4.v[i],ar4);
        vfmaaeq(p1[j].i.v[i],lam_2.v[i],ai2,lam_4.v[i],ai4);
        }
      Tv ar3=vload(creal(alm[njobs*(l+1)+j])),
         ai3=vload(cimag(alm[njobs*(l+1)+j])),
         ar1=vload(creal(alm[njobs*(l+3)+j])),
         ai1=vload(cimag(alm[njobs*(l+3)+j]));
      for (int i=0; i<nvec; ++i)
        {
        vfmaaeq(p2[j].r.v[i],lam_3.v[i],ar3,lam_1.v[i],ar1);
        vfmaaeq(p2[j].i.v[i],lam_3.v[i],ai3,lam_1.v[i],ai1);
        }
      }
    r0=vload(rf[l+3].f[0]);r1=vload(rf[l+3].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_4.v[i],r1));
    l+=4;
    }
  }
  while (l<lmax)
    {
    Tv r0=vload(rf[l].f[0]),r1=vload(rf[l].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    for (int j=0; j<njobs; ++j)
      {
      Tv ar=vload(creal(alm[njobs*l+j])),
         ai=vload(cimag(alm[njobs*l+j]));
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(p1[j].r.v[i],lam_2.v[i],ar);
        vfmaeq(p1[j].i.v[i],lam_2.v[i],ai);
        }
      ar=vload(creal(alm[njobs*(l+1)+j]));
      ai=vload(cimag(alm[njobs*(l+1)+j]));
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(p2[j].r.v[i],lam_1.v[i],ar);
        vfmaeq(p2[j].i.v[i],lam_1.v[i],ai);
        }
      }
    r0=vload(rf[l+1].f[0]);r1=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_2.v[i],r1));
    l+=2;
    }
  if (l==lmax)
    {
    for (int j=0; j<njobs; ++j)
      {
      Tv ar=vload(creal(alm[njobs*l+j])),ai=vload(cimag(alm[njobs*l+j]));
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(p1[j].r.v[i],lam_2.v[i],ar);
        vfmaeq(p1[j].i.v[i],lam_2.v[i],ai);
        }
      }
    }
  }

static void Z(map2alm_kernel) (const Tb cth, const Y(Tbri) * restrict p1,
  const Y(Tbri) * restrict p2, Tb lam_1, Tb lam_2,
  const sharp_ylmgen_dbl2 * restrict rf, dcmplx * restrict alm, int l, int lmax
  NJ1)
  {
  while (l<lmax)
    {
    Tv r0=vload(rf[l].f[0]),r1=vload(rf[l].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    for (int j=0; j<njobs; ++j)
      {
      Tv tr1=vzero, ti1=vzero, tr2=vzero, ti2=vzero;
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(tr1,lam_2.v[i],p1[j].r.v[i]);
        vfmaeq(ti1,lam_2.v[i],p1[j].i.v[i]);
        }
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(tr2,lam_1.v[i],p2[j].r.v[i]);
        vfmaeq(ti2,lam_1.v[i],p2[j].i.v[i]);
        }
      vhsum_cmplx2(tr1,ti1,tr2,ti2,&alm[l*njobs+j],&alm[(l+1)*njobs+j]);
      }
    r0=vload(rf[l+1].f[0]);r1=vload(rf[l+1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_2.v[i],r1));
    l+=2;
    }
  if (l==lmax)
    {
    for (int j=0; j<njobs; ++j)
      {
      Tv tre=vzero, tim=vzero;
      for (int i=0; i<nvec; ++i)
        {
        vfmaeq(tre,lam_2.v[i],p1[j].r.v[i]);
        vfmaeq(tim,lam_2.v[i],p1[j].i.v[i]);
        }
      alm[l*njobs+j]+=vhsum_cmplx(tre,tim);
      }
    }
  }

static void Z(calc_alm2map) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbri) * restrict p1,
  Y(Tbri) * restrict p2 NJ1)
  {
  int l,lmax=gen->lmax;
  Tb lam_1,lam_2,scale;
  Y(iter_to_ieee) (sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * (4+4*njobs)*VLEN*nvec;

  Tb corfac;
  Y(getCorfac)(scale,&corfac,gen->cf);
  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scale,sharp_minscale);
  while (!full_ieee)
    {
    for (int j=0; j<njobs; ++j)
      {
      Tv ar=vload(creal(alm[njobs*l+j])),ai=vload(cimag(alm[njobs*l+j]));
      for (int i=0; i<nvec; ++i)
        {
        Tv tmp=vmul(lam_2.v[i],corfac.v[i]);
        vfmaeq(p1[j].r.v[i],tmp,ar);
        vfmaeq(p1[j].i.v[i],tmp,ai);
        }
      }
    if (++l>lmax) break;
    Tv r0=vload(rf[l-1].f[0]),r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    for (int j=0; j<njobs; ++j)
      {
      Tv ar=vload(creal(alm[njobs*l+j])),ai=vload(cimag(alm[njobs*l+j]));
      for (int i=0; i<nvec; ++i)
        {
        Tv tmp=vmul(lam_1.v[i],corfac.v[i]);
        vfmaeq(p2[j].r.v[i],tmp,ar);
        vfmaeq(p2[j].i.v[i],tmp,ai);
        }
      }
    if (++l>lmax) break;
    r0=vload(rf[l-1].f[0]); r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_2.v[i],r1));
    if (Y(rescale)(&lam_1,&lam_2,&scale))
      {
      Y(getCorfac)(scale,&corfac,gen->cf);
      full_ieee = Y(TballGe)(scale,sharp_minscale);
      }
    }
  if (l>lmax) return;

  Y(Tbmuleq)(&lam_1,corfac); Y(Tbmuleq)(&lam_2,corfac);
  Z(alm2map_kernel) (cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax NJ2);
  }

static void Z(calc_map2alm) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, const Y(Tbri) * restrict p1,
  const Y(Tbri) * restrict p2 NJ1)
  {
  int lmax=gen->lmax;
  Tb lam_1,lam_2,scale;
  int l=gen->m;
  Y(iter_to_ieee) (sth,cth,&l,&lam_1,&lam_2,&scale,gen);
  job->opcnt += (l-gen->m) * 4*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * (4+4*njobs)*VLEN*nvec;

  const sharp_ylmgen_dbl2 * restrict rf = gen->rf;
  Tb corfac;
  Y(getCorfac)(scale,&corfac,gen->cf);
  dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scale,sharp_minscale);
  while (!full_ieee)
    {
    for (int j=0; j<njobs; ++j)
      {
      Tv tre=vzero, tim=vzero;
      for (int i=0; i<nvec; ++i)
        {
        Tv tmp=vmul(lam_2.v[i],corfac.v[i]);
        vfmaeq(tre,tmp,p1[j].r.v[i]);
        vfmaeq(tim,tmp,p1[j].i.v[i]);
        }
      alm[l*njobs+j]+=vhsum_cmplx(tre,tim);
      }
    if (++l>lmax) return;
    Tv r0=vload(rf[l-1].f[0]),r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_1.v[i] = vsub(vmul(vmul(cth.v[i],lam_2.v[i]),r0),vmul(lam_1.v[i],r1));
    for (int j=0; j<njobs; ++j)
      {
      Tv tre=vzero, tim=vzero;
      for (int i=0; i<nvec; ++i)
        {
        Tv tmp=vmul(lam_1.v[i],corfac.v[i]);
        vfmaeq(tre,tmp,p2[j].r.v[i]);
        vfmaeq(tim,tmp,p2[j].i.v[i]);
        }
      alm[l*njobs+j]+=vhsum_cmplx(tre,tim);
      }
    if (++l>lmax) return;
    r0=vload(rf[l-1].f[0]); r1=vload(rf[l-1].f[1]);
    for (int i=0; i<nvec; ++i)
      lam_2.v[i] = vsub(vmul(vmul(cth.v[i],lam_1.v[i]),r0),vmul(lam_2.v[i],r1));
    if (Y(rescale)(&lam_1,&lam_2,&scale))
      {
      Y(getCorfac)(scale,&corfac,gen->cf);
      full_ieee = Y(TballGe)(scale,sharp_minscale);
      }
    }

  Y(Tbmuleq)(&lam_1,corfac); Y(Tbmuleq)(&lam_2,corfac);
  Z(map2alm_kernel) (cth, p1, p2, lam_1, lam_2, rf, alm, l, lmax NJ2);
  }

static inline void Z(saddstep) (Y(Tbqu) * restrict px, Y(Tbqu) * restrict py,
  const Tb rxp, const Tb rxm, const dcmplx * restrict alm NJ1)
  {
  for (int j=0; j<njobs; ++j)
    {
    Tv agr=vload(creal(alm[2*j])), agi=vload(cimag(alm[2*j])),
       acr=vload(creal(alm[2*j+1])), aci=vload(cimag(alm[2*j+1]));
    for (int i=0; i<nvec; ++i)
      {
      Tv lw=vadd(rxp.v[i],rxm.v[i]);
      vfmaeq(px[j].qr.v[i],agr,lw);
      vfmaeq(px[j].qi.v[i],agi,lw);
      vfmaeq(px[j].ur.v[i],acr,lw);
      vfmaeq(px[j].ui.v[i],aci,lw);
      }
    for (int i=0; i<nvec; ++i)
      {
      Tv lx=vsub(rxm.v[i],rxp.v[i]);
      vfmseq(py[j].qr.v[i],aci,lx);
      vfmaeq(py[j].qi.v[i],acr,lx);
      vfmaeq(py[j].ur.v[i],agi,lx);
      vfmseq(py[j].ui.v[i],agr,lx);
      }
    }
  }

static inline void Z(saddstepb) (Y(Tbqu) * restrict p1, Y(Tbqu) * restrict p2,
  const Tb r1p, const Tb r1m, const Tb r2p, const Tb r2m,
  const dcmplx * restrict alm1, const dcmplx * restrict alm2 NJ1)
  {
  for (int j=0; j<njobs; ++j)
    {
    Tv agr1=vload(creal(alm1[2*j])), agi1=vload(cimag(alm1[2*j])),
       acr1=vload(creal(alm1[2*j+1])), aci1=vload(cimag(alm1[2*j+1]));
    Tv agr2=vload(creal(alm2[2*j])), agi2=vload(cimag(alm2[2*j])),
       acr2=vload(creal(alm2[2*j+1])), aci2=vload(cimag(alm2[2*j+1]));
    for (int i=0; i<nvec; ++i)
      {
      Tv lw1=vadd(r2p.v[i],r2m.v[i]);
      Tv lx2=vsub(r1m.v[i],r1p.v[i]);
      vfmaseq(p1[j].qr.v[i],agr1,lw1,aci2,lx2);
      vfmaaeq(p1[j].qi.v[i],agi1,lw1,acr2,lx2);
      vfmaaeq(p1[j].ur.v[i],acr1,lw1,agi2,lx2);
      vfmaseq(p1[j].ui.v[i],aci1,lw1,agr2,lx2);
      }
    for (int i=0; i<nvec; ++i)
      {
      Tv lx1=vsub(r2m.v[i],r2p.v[i]);
      Tv lw2=vadd(r1p.v[i],r1m.v[i]);
      vfmaseq(p2[j].qr.v[i],agr2,lw2,aci1,lx1);
      vfmaaeq(p2[j].qi.v[i],agi2,lw2,acr1,lx1);
      vfmaaeq(p2[j].ur.v[i],acr2,lw2,agi1,lx1);
      vfmaseq(p2[j].ui.v[i],aci2,lw2,agr1,lx1);
      }
    }
  }

static inline void Z(saddstep2) (const Y(Tbqu) * restrict px,
  const Y(Tbqu) * restrict py, const Tb * restrict rxp,
  const Tb * restrict rxm, dcmplx * restrict alm NJ1)
  {
  for (int j=0; j<njobs; ++j)
    {
    Tv agr=vzero, agi=vzero, acr=vzero, aci=vzero;
    for (int i=0; i<nvec; ++i)
      {
      Tv lw=vadd(rxp->v[i],rxm->v[i]);
      vfmaeq(agr,px[j].qr.v[i],lw);
      vfmaeq(agi,px[j].qi.v[i],lw);
      vfmaeq(acr,px[j].ur.v[i],lw);
      vfmaeq(aci,px[j].ui.v[i],lw);
      }
    for (int i=0; i<nvec; ++i)
      {
      Tv lx=vsub(rxm->v[i],rxp->v[i]);
      vfmseq(agr,py[j].ui.v[i],lx);
      vfmaeq(agi,py[j].ur.v[i],lx);
      vfmaeq(acr,py[j].qi.v[i],lx);
      vfmseq(aci,py[j].qr.v[i],lx);
      }
    vhsum_cmplx2(agr,agi,acr,aci,&alm[2*j],&alm[2*j+1]);
    }
  }

static void Z(alm2map_spin_kernel) (Tb cth, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm, int l,
  int lmax NJ1)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec2p.v[i])),
                        vmul(fx2,rec1p.v[i]));
      rec1m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec2m.v[i])),
                        vmul(fx2,rec1m.v[i]));
      }
    Z(saddstepb)(p1,p2,rec1p,rec1m,rec2p,rec2m,&alm[2*njobs*l],
      &alm[2*njobs*(l+1)] NJ2);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec1p.v[i])),
                        vmul(fx2,rec2p.v[i]));
      rec2m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec1m.v[i])),
                        vmul(fx2,rec2m.v[i]));
      }
    l+=2;
    }
  if (l==lmax)
    Z(saddstep)(p1, p2, rec2p, rec2m, &alm[2*njobs*l] NJ2);
  }

static void Z(map2alm_spin_kernel) (Tb cth, const Y(Tbqu) * restrict p1,
  const Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, dcmplx * restrict alm, int l, int lmax
  NJ1)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec2p.v[i])),
                        vmul(fx2,rec1p.v[i]));
      rec1m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec2m.v[i])),
                        vmul(fx2,rec1m.v[i]));
      }
    Z(saddstep2)(p1, p2, &rec2p, &rec2m, &alm[2*njobs*l] NJ2);
    Z(saddstep2)(p2, p1, &rec1p, &rec1m, &alm[2*njobs*(l+1)] NJ2);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec1p.v[i])),
                        vmul(fx2,rec2p.v[i]));
      rec2m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec1m.v[i])),
                        vmul(fx2,rec2m.v[i]));
      }
    l+=2;
    }
  if (l==lmax)
    Z(saddstep2)(p1, p2, &rec2p, &rec2m, &alm[2*njobs*l] NJ2);
  }

static void Z(calc_alm2map_spin) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2 NJ1)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * (12+16*njobs)*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Z(saddstep)(p1, p2, Y(Tbprod)(rec2p,corfacp), Y(Tbprod)(rec2m,corfacm),
      &alm[2*njobs*l] NJ2);
    if (++l>lmax) break;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    Z(saddstep)(p2, p1, Y(Tbprod)(rec1p,corfacp), Y(Tbprod)(rec1m,corfacm),
      &alm[2*njobs*l] NJ2);
    if (++l>lmax) break;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  if (l>lmax) return;

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Z(alm2map_spin_kernel) (cth, p1, p2, rec1p, rec1m, rec2p, rec2m, fx, alm, l,
    lmax NJ2);
  }

static void Z(calc_map2alm_spin) (Tb cth, Tb sth,
  const sharp_Ylmgen_C * restrict gen, sharp_job *job,
  const Y(Tbqu) * restrict p1, const Y(Tbqu) * restrict p2 NJ1)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * (12+16*njobs)*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Tb t1=Y(Tbprod)(rec2p,corfacp), t2=Y(Tbprod)(rec2m,corfacm);
    Z(saddstep2)(p1, p2, &t1, &t2, &alm[2*njobs*l] NJ2);
    if (++l>lmax) return;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    t1=Y(Tbprod)(rec1p,corfacp); t2=Y(Tbprod)(rec1m,corfacm);
    Z(saddstep2)(p2, p1, &t1, &t2, &alm[2*njobs*l] NJ2);
    if (++l>lmax) return;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Z(map2alm_spin_kernel)(cth,p1,p2,rec1p,rec1m,rec2p,rec2m,fx,alm,l,lmax NJ2);
  }

static inline void Z(saddstep_d) (Y(Tbqu) * restrict px, Y(Tbqu) * restrict py,
  const Tb rxp, const Tb rxm, const dcmplx * restrict alm NJ1)
  {
  for (int j=0; j<njobs; ++j)
    {
    Tv ar=vload(creal(alm[j])), ai=vload(cimag(alm[j]));
    for (int i=0; i<nvec; ++i)
      {
      Tv lw=vadd(rxp.v[i],rxm.v[i]);
      vfmaeq(px[j].qr.v[i],ar,lw);
      vfmaeq(px[j].qi.v[i],ai,lw);
      }
    for (int i=0; i<nvec; ++i)
      {
      Tv lx=vsub(rxm.v[i],rxp.v[i]);
      vfmaeq(py[j].ur.v[i],ai,lx);
      vfmseq(py[j].ui.v[i],ar,lx);
      }
    }
  }

static void Z(alm2map_deriv1_kernel) (Tb cth, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2, Tb rec1p, Tb rec1m, Tb rec2p, Tb rec2m,
  const sharp_ylmgen_dbl3 * restrict fx, const dcmplx * restrict alm, int l,
  int lmax NJ1)
  {
  while (l<lmax)
    {
    Tv fx0=vload(fx[l+1].f[0]),fx1=vload(fx[l+1].f[1]),
       fx2=vload(fx[l+1].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec1p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec2p.v[i])),
                        vmul(fx2,rec1p.v[i]));
      rec1m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec2m.v[i])),
                        vmul(fx2,rec1m.v[i]));
      }
    Z(saddstep_d)(p1,p2,rec2p,rec2m,&alm[njobs*l] NJ2);
    Z(saddstep_d)(p2,p1,rec1p,rec1m,&alm[njobs*(l+1)] NJ2);
    fx0=vload(fx[l+2].f[0]);fx1=vload(fx[l+2].f[1]);
    fx2=vload(fx[l+2].f[2]);
    for (int i=0; i<nvec; ++i)
      {
      rec2p.v[i] = vsub(vmul(vsub(cth.v[i],fx1),vmul(fx0,rec1p.v[i])),
                        vmul(fx2,rec2p.v[i]));
      rec2m.v[i] = vsub(vmul(vadd(cth.v[i],fx1),vmul(fx0,rec1m.v[i])),
                        vmul(fx2,rec2m.v[i]));
      }
    l+=2;
    }
  if (l==lmax)
    Z(saddstep_d)(p1, p2, rec2p, rec2m, &alm[njobs*l] NJ2);
  }

static void Z(calc_alm2map_deriv1) (const Tb cth, const Tb sth,
  const sharp_Ylmgen_C *gen, sharp_job *job, Y(Tbqu) * restrict p1,
  Y(Tbqu) * restrict p2 NJ1)
  {
  int l, lmax=gen->lmax;
  Tb rec1p, rec1m, rec2p, rec2m, scalem, scalep;
  Y(iter_to_ieee_spin)
    (cth,sth,&l,&rec1p,&rec1m,&rec2p,&rec2m,&scalep,&scalem,gen);
  job->opcnt += (l-gen->m) * 10*VLEN*nvec;
  if (l>lmax) return;
  job->opcnt += (lmax+1-l) * (12+8*njobs)*VLEN*nvec;

  const sharp_ylmgen_dbl3 * restrict fx = gen->fx;
  Tb corfacp,corfacm;
  Y(getCorfac)(scalep,&corfacp,gen->cf);
  Y(getCorfac)(scalem,&corfacm,gen->cf);
  const dcmplx * restrict alm=job->almtmp;
  int full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
  while (!full_ieee)
    {
    Z(saddstep_d)(p1, p2, Y(Tbprod)(rec2p,corfacp), Y(Tbprod)(rec2m,corfacm),
      &alm[njobs*l] NJ2);
    if (++l>lmax) break;
    Y(rec_step)(&rec1p,&rec1m,&rec2p,&rec2m,cth,fx[l]);
    Z(saddstep_d)(p2, p1, Y(Tbprod)(rec1p,corfacp), Y(Tbprod)(rec1m,corfacm),
      &alm[njobs*l] NJ2);
    if (++l>lmax) break;
    Y(rec_step)(&rec2p,&rec2m,&rec1p,&rec1m,cth,fx[l]);
    if (Y(rescale)(&rec1p,&rec2p,&scalep) | Y(rescale)(&rec1m,&rec2m,&scalem))
      {
      Y(getCorfac)(scalep,&corfacp,gen->cf);
      Y(getCorfac)(scalem,&corfacm,gen->cf);
      full_ieee = Y(TballGe)(scalep,sharp_minscale)
               && Y(TballGe)(scalem,sharp_minscale);
      }
    }

  if (l>lmax) return;

  Y(Tbmuleq)(&rec1p,corfacp); Y(Tbmuleq)(&rec2p,corfacp);
  Y(Tbmuleq)(&rec1m,corfacm); Y(Tbmuleq)(&rec2m,corfacm);
  Z(alm2map_deriv1_kernel) (cth, p1, p2, rec1p, rec1m, rec2p, rec2m, fx, alm, l,
    lmax NJ2);
  }


#define VZERO(var) do { memset(&(var),0,sizeof(var)); } while(0)

static void Z(inner_loop) (sharp_job *job, const int *ispair,
  const double *cth_, const double *sth_, int llim, int ulim,
  sharp_Ylmgen_C *gen, int mi, const int *mlim NJ1)
  {
  const int nval=nvec*VLEN;
  const int m = job->ainfo->mval[mi];
  sharp_Ylmgen_prepare (gen, m);

  switch (job->type)
    {
    case SHARP_ALM2MAP:
    case SHARP_ALM2MAP_DERIV1:
      {
      if (job->spin==0)
        {
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tburi) p1[njobs],p2[njobs]; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;

          int skip=1;
          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            }
          if (!skip)
            Z(calc_alm2map) (cth.b,sth.b,gen,job,&p1[0].b,&p2[0].b NJ2);

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot<ulim-llim)
              {
              for (int j=0; j<njobs; ++j)
                {
                int phas_idx = itot*job->s_th + mi*job->s_m + 2*j;
                complex double r1 = p1[j].s.r[i] + p1[j].s.i[i]*_Complex_I,
                               r2 = p2[j].s.r[i] + p2[j].s.i[i]*_Complex_I;
                job->phase[phas_idx] = r1+r2;
                if (ispair[itot])
                  job->phase[phas_idx+1] = r1-r2;
                }
              }
            }
          }
        }
      else
        {
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tbuqu) p1[njobs],p2[njobs]; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            }
          if (!skip)
            (job->type==SHARP_ALM2MAP) ?
              Z(calc_alm2map_spin  )
                (cth.b,sth.b,gen,job,&p1[0].b,&p2[0].b NJ2) :
              Z(calc_alm2map_deriv1)
                (cth.b,sth.b,gen,job,&p1[0].b,&p2[0].b NJ2);

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot<ulim-llim)
              {
              for (int j=0; j<njobs; ++j)
                {
                int phas_idx = itot*job->s_th + mi*job->s_m + 4*j;
                complex double q1 = p1[j].s.qr[i] + p1[j].s.qi[i]*_Complex_I,
                               q2 = p2[j].s.qr[i] + p2[j].s.qi[i]*_Complex_I,
                               u1 = p1[j].s.ur[i] + p1[j].s.ui[i]*_Complex_I,
                               u2 = p2[j].s.ur[i] + p2[j].s.ui[i]*_Complex_I;
                job->phase[phas_idx] = q1+q2;
                job->phase[phas_idx+2] = u1+u2;
                if (ispair[itot])
                  {
                  dcmplx *phQ = &(job->phase[phas_idx+1]),
                         *phU = &(job->phase[phas_idx+3]);
                  *phQ = q1-q2;
                  *phU = u1-u2;
                  if ((gen->mhi-gen->m+gen->s)&1)
                    { *phQ=-(*phQ); *phU=-(*phU); }
                  }
                }
              }
            }
          }
        }
      break;
      }
    case SHARP_MAP2ALM:
      {
      if (job->spin==0)
        {
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tburi) p1[njobs], p2[njobs]; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            if ((i+ith<ulim-llim)&&(mlim[itot]>=m))
              {
              for (int j=0; j<njobs; ++j)
                {
                int phas_idx = itot*job->s_th + mi*job->s_m + 2*j;
                dcmplx ph1=job->phase[phas_idx];
                dcmplx ph2=ispair[itot] ? job->phase[phas_idx+1] : 0.;
                p1[j].s.r[i]=creal(ph1+ph2); p1[j].s.i[i]=cimag(ph1+ph2);
                p2[j].s.r[i]=creal(ph1-ph2); p2[j].s.i[i]=cimag(ph1-ph2);
                }
              }
            }
          if (!skip)
            Z(calc_map2alm)(cth.b,sth.b,gen,job,&p1[0].b,&p2[0].b NJ2);
          }
        }
      else
        {
        for (int ith=0; ith<ulim-llim; ith+=nval)
          {
          Y(Tbuqu) p1[njobs], p2[njobs]; VZERO(p1); VZERO(p2);
          Y(Tbu) cth, sth;
          int skip=1;

          for (int i=0; i<nval; ++i)
            {
            int itot=i+ith;
            if (itot>=ulim-llim) itot=ulim-llim-1;
            if (mlim[itot]>=m) skip=0;
            cth.s[i]=cth_[itot]; sth.s[i]=sth_[itot];
            if (i+ith<ulim-llim)
              {
              for (int j=0; j<njobs; ++j)
                {
                int phas_idx = itot*job->s_th + mi*job->s_m + 4*j;
                dcmplx p1Q=job->phase[phas_idx],
                       p1U=job->phase[phas_idx+2],
                       p2Q=ispair[itot] ? job->phase[phas_idx+1]:0.,
                       p2U=ispair[itot] ? job->phase[phas_idx+3]:0.;
                if ((gen->mhi-gen->m+gen->s)&1)
                  { p2Q=-p2Q; p2U=-p2U; }
                p1[j].s.qr[i]=creal(p1Q+p2Q); p1[j].s.qi[i]=cimag(p1Q+p2Q);
                p1[j].s.ur[i]=creal(p1U+p2U); p1[j].s.ui[i]=cimag(p1U+p2U);
                p2[j].s.qr[i]=creal(p1Q-p2Q); p2[j].s.qi[i]=cimag(p1Q-p2Q);
                p2[j].s.ur[i]=creal(p1U-p2U); p2[j].s.ui[i]=cimag(p1U-p2U);
                }
              }
            }
          if (!skip)
            Z(calc_map2alm_spin) (cth.b,sth.b,gen,job,&p1[0].b,&p2[0].b NJ2);
          }
        }
      break;
      }
    default:
      {
      UTIL_FAIL("must not happen");
      break;
      }
    }
  }

#undef VZERO

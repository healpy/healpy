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
  fftpack.c : A set of FFT routines in C.
  Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber
  (Version 4, 1985).

  C port by Martin Reinecke (2010)
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fftpack.h"

#define WA(x,i) wa[(i)+(x)*ido]
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define PM(a,b,c,d) { a=c+d; b=c-d; }
#define PMC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
#define ADDC(a,b,c) { a.r=b.r+c.r; a.i=b.i+c.i; }
#define SCALEC(a,b) { a.r*=b; a.i*=b; }
#define CONJFLIPC(a) { double tmp_=a.r; a.r=-a.i; a.i=tmp_; }
/* (a+ib) = conj(c+id) * (e+if) */
#define MULPM(a,b,c,d,e,f) { a=c*e+d*f; b=c*f-d*e; }

typedef struct {
  double r,i;
} cmplx;

#define CONCAT(a,b) a ## b

#define X(arg) CONCAT(passb,arg)
#define BACKWARD
#include "fftpack_inc.c"
#undef BACKWARD
#undef X

#define X(arg) CONCAT(passf,arg)
#include "fftpack_inc.c"
#undef X

#undef CC
#undef CH
#define CC(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+cdim*(c))]

static void radf2 (size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=2;
  size_t i, k, ic;
  double ti2, tr2;

  for (k=0; k<l1; k++)
    PM (CH(0,0,k),CH(ido-1,1,k),CC(0,k,0),CC(0,k,1))
  if ((ido&1)==0)
    for (k=0; k<l1; k++)
      {
      CH(    0,1,k) = -CC(ido-1,k,1);
      CH(ido-1,0,k) =  CC(ido-1,k,0);
      }
  if (ido<=2) return;
  for (k=0; k<l1; k++)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      MULPM (tr2,ti2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      PM (CH(i-1,0,k),CH(ic-1,1,k),CC(i-1,k,0),tr2)
      PM (CH(i  ,0,k),CH(ic  ,1,k),ti2,CC(i  ,k,0))
      }
  }

static void radf3(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=3;
  static const double taur=-0.5, taui=0.86602540378443864676;
  size_t i, k, ic;
  double ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;

  for (k=0; k<l1; k++)
    {
    cr2=CC(0,k,1)+CC(0,k,2);
    CH(0,0,k) = CC(0,k,0)+cr2;
    CH(0,2,k) = taui*(CC(0,k,2)-CC(0,k,1));
    CH(ido-1,1,k) = CC(0,k,0)+taur*cr2;
    }
  if (ido==1) return;
  for (k=0; k<l1; k++)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2))
      cr2=dr2+dr3;
      ci2=di2+di3;
      CH(i-1,0,k) = CC(i-1,k,0)+cr2;
      CH(i  ,0,k) = CC(i  ,k,0)+ci2;
      tr2 = CC(i-1,k,0)+taur*cr2;
      ti2 = CC(i  ,k,0)+taur*ci2;
      tr3 = taui*(di2-di3);
      ti3 = taui*(dr3-dr2);
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr3)
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti3,ti2)
      }
  }

static void radf4(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=4;
  static const double hsqt2=0.70710678118654752440;
  size_t i, k, ic;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

  for (k=0; k<l1; k++)
    {
    PM (tr1,CH(0,2,k),CC(0,k,3),CC(0,k,1))
    PM (tr2,CH(ido-1,1,k),CC(0,k,0),CC(0,k,2))
    PM (CH(0,0,k),CH(ido-1,3,k),tr2,tr1)
    }
  if ((ido&1)==0)
    for (k=0; k<l1; k++)
      {
      ti1=-hsqt2*(CC(ido-1,k,1)+CC(ido-1,k,3));
      tr1= hsqt2*(CC(ido-1,k,1)-CC(ido-1,k,3));
      PM (CH(ido-1,0,k),CH(ido-1,2,k),CC(ido-1,k,0),tr1)
      PM (CH(    0,3,k),CH(    0,1,k),ti1,CC(ido-1,k,2))
      }
  if (ido<=2) return;
  for (k=0; k<l1; k++)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      MULPM(cr2,ci2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      MULPM(cr3,ci3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2))
      MULPM(cr4,ci4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3))
      PM(tr1,tr4,cr4,cr2)
      PM(ti1,ti4,ci2,ci4)
      PM(tr2,tr3,CC(i-1,k,0),cr3)
      PM(ti2,ti3,CC(i  ,k,0),ci3)
      PM(CH(i-1,0,k),CH(ic-1,3,k),tr2,tr1)
      PM(CH(i  ,0,k),CH(ic  ,3,k),ti1,ti2)
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr3,ti4)
      PM(CH(i  ,2,k),CH(ic  ,1,k),tr4,ti3)
      }
  }

static void radf5(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=5;
  static const double tr11= 0.3090169943749474241, ti11=0.95105651629515357212,
                      tr12=-0.8090169943749474241, ti12=0.58778525229247312917;
  size_t i, k, ic;
  double ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3,
         dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;

  for (k=0; k<l1; k++)
    {
    PM (cr2,ci5,CC(0,k,4),CC(0,k,1))
    PM (cr3,ci4,CC(0,k,3),CC(0,k,2))
    CH(0,0,k)=CC(0,k,0)+cr2+cr3;
    CH(ido-1,1,k)=CC(0,k,0)+tr11*cr2+tr12*cr3;
    CH(0,2,k)=ti11*ci5+ti12*ci4;
    CH(ido-1,3,k)=CC(0,k,0)+tr12*cr2+tr11*cr3;
    CH(0,4,k)=ti12*ci5-ti11*ci4;
    }
  if (ido==1) return;
  for (k=0; k<l1;++k)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1))
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2))
      MULPM (dr4,di4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3))
      MULPM (dr5,di5,WA(3,i-2),WA(3,i-1),CC(i-1,k,4),CC(i,k,4))
      PM(cr2,ci5,dr5,dr2)
      PM(ci2,cr5,di2,di5)
      PM(cr3,ci4,dr4,dr3)
      PM(ci3,cr4,di3,di4)
      CH(i-1,0,k)=CC(i-1,k,0)+cr2+cr3;
      CH(i  ,0,k)=CC(i  ,k,0)+ci2+ci3;
      tr2=CC(i-1,k,0)+tr11*cr2+tr12*cr3;
      ti2=CC(i  ,k,0)+tr11*ci2+tr12*ci3;
      tr3=CC(i-1,k,0)+tr12*cr2+tr11*cr3;
      ti3=CC(i  ,k,0)+tr12*ci2+tr11*ci3;
      MULPM(tr5,tr4,cr5,cr4,ti11,ti12)
      MULPM(ti5,ti4,ci5,ci4,ti11,ti12)
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr5)
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti5,ti2)
      PM(CH(i-1,4,k),CH(ic-1,3,k),tr3,tr4)
      PM(CH(i  ,4,k),CH(ic  ,3,k),ti4,ti3)
      }
  }

#undef CH
#undef CC
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define C1(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define C2(a,b) cc[(a)+idl1*(b)]
#define CH2(a,b) ch[(a)+idl1*(b)]
static void radfg(size_t ido, size_t ip, size_t l1, size_t idl1,
  double *cc, double *ch, const double *wa)
  {
  const size_t cdim=ip;
  static const double twopi=6.28318530717958647692;
  size_t idij, ipph, i, j, k, l, j2, ic, jc, lc, ik;
  double ai1, ai2, ar1, ar2, arg;
  double *csarr;
  size_t aidx;

  ipph=(ip+1)/ 2;
  if(ido!=1)
    {
    memcpy(ch,cc,idl1*sizeof(double));

    for(j=1; j<ip; j++)
      for(k=0; k<l1; k++)
        {
        CH(0,k,j)=C1(0,k,j);
        idij=(j-1)*ido+1;
        for(i=2; i<ido; i+=2,idij+=2)
          MULPM(CH(i-1,k,j),CH(i,k,j),wa[idij-1],wa[idij],C1(i-1,k,j),C1(i,k,j))
        }

    for(j=1,jc=ip-1; j<ipph; j++,jc--)
      for(k=0; k<l1; k++)
        for(i=2; i<ido; i+=2)
          {
          PM(C1(i-1,k,j),C1(i  ,k,jc),CH(i-1,k,jc),CH(i-1,k,j ))
          PM(C1(i  ,k,j),C1(i-1,k,jc),CH(i  ,k,j ),CH(i  ,k,jc))
          }
    }
  else
    memcpy(cc,ch,idl1*sizeof(double));

  for(j=1,jc=ip-1; j<ipph; j++,jc--)
    for(k=0; k<l1; k++)
      PM(C1(0,k,j),C1(0,k,jc),CH(0,k,jc),CH(0,k,j))

  csarr=RALLOC(double,2*ip);
  arg=twopi / ip;
  csarr[0]=1.;
  csarr[1]=0.;
  csarr[2]=csarr[2*ip-2]=cos(arg);
  csarr[3]=sin(arg); csarr[2*ip-1]=-csarr[3];
  for (i=2; i<=ip/2; ++i)
    {
    csarr[2*i]=csarr[2*ip-2*i]=cos(i*arg);
    csarr[2*i+1]=sin(i*arg);
    csarr[2*ip-2*i+1]=-csarr[2*i+1];
    }
  for(l=1,lc=ip-1; l<ipph; l++,lc--)
    {
    ar1=csarr[2*l];
    ai1=csarr[2*l+1];
    for(ik=0; ik<idl1; ik++)
      {
      CH2(ik,l)=C2(ik,0)+ar1*C2(ik,1);
      CH2(ik,lc)=ai1*C2(ik,ip-1);
      }
    aidx=2*l;
    for(j=2,jc=ip-2; j<ipph; j++,jc--)
      {
      aidx+=2*l;
      if (aidx>=2*ip) aidx-=2*ip;
      ar2=csarr[aidx];
      ai2=csarr[aidx+1];
      for(ik=0; ik<idl1; ik++)
        {
        CH2(ik,l )+=ar2*C2(ik,j );
        CH2(ik,lc)+=ai2*C2(ik,jc);
        }
      }
    }
  DEALLOC(csarr);

  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      CH2(ik,0)+=C2(ik,j);

  for(k=0; k<l1; k++)
    memcpy(&CC(0,0,k),&CH(0,k,0),ido*sizeof(double));
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      {
      CC(ido-1,j2-1,k) = CH(0,k,j );
      CC(0    ,j2  ,k) = CH(0,k,jc);
      }
    }
  if(ido==1) return;

  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        PM (CC(i-1,j2,k),CC(ic-1,j2-1,k),CH(i-1,k,j ),CH(i-1,k,jc))
        PM (CC(i  ,j2,k),CC(ic  ,j2-1,k),CH(i  ,k,jc),CH(i  ,k,j ))
        }
    }
  }

#undef CC
#undef CH
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]

static void radb2(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=2;
  size_t i, k, ic;
  double ti2, tr2;

  for (k=0; k<l1; k++)
    PM (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(ido-1,1,k))
  if ((ido&1)==0)
    for (k=0; k<l1; k++)
      {
      CH(ido-1,k,0) =  2*CC(ido-1,0,k);
      CH(ido-1,k,1) = -2*CC(0    ,1,k);
      }
  if (ido<=2) return;
  for (k=0; k<l1;++k)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      PM (CH(i-1,k,0),tr2,CC(i-1,0,k),CC(ic-1,1,k))
      PM (ti2,CH(i  ,k,0),CC(i  ,0,k),CC(ic  ,1,k))
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ti2,tr2)
      }
  }

static void radb3(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=3;
  static const double taur=-0.5, taui=0.86602540378443864676;
  size_t i, k, ic;
  double ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

  for (k=0; k<l1; k++)
    {
    tr2=2*CC(ido-1,1,k);
    cr2=CC(0,0,k)+taur*tr2;
    CH(0,k,0)=CC(0,0,k)+tr2;
    ci3=2*taui*CC(0,2,k);
    PM (CH(0,k,2),CH(0,k,1),cr2,ci3);
    }
  if (ido==1) return;
  for (k=0; k<l1; k++)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      tr2=CC(i-1,2,k)+CC(ic-1,1,k);
      ti2=CC(i  ,2,k)-CC(ic  ,1,k);
      cr2=CC(i-1,0,k)+taur*tr2;
      ci2=CC(i  ,0,k)+taur*ti2;
      CH(i-1,k,0)=CC(i-1,0,k)+tr2;
      CH(i  ,k,0)=CC(i  ,0,k)+ti2;
      cr3=taui*(CC(i-1,2,k)-CC(ic-1,1,k));
      ci3=taui*(CC(i  ,2,k)+CC(ic  ,1,k));
      PM(dr3,dr2,cr2,ci3)
      PM(di2,di3,ci2,cr3)
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2)
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3)
      }
  }

static void radb4(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=4;
  static const double sqrt2=1.41421356237309504880;
  size_t i, k, ic;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

  for (k=0; k<l1; k++)
    {
    PM (tr2,tr1,CC(0,0,k),CC(ido-1,3,k))
    tr3=2*CC(ido-1,1,k);
    tr4=2*CC(0,2,k);
    PM (CH(0,k,0),CH(0,k,2),tr2,tr3)
    PM (CH(0,k,3),CH(0,k,1),tr1,tr4)
    }
  if ((ido&1)==0)
    for (k=0; k<l1; k++)
      {
      PM (ti1,ti2,CC(0    ,3,k),CC(0    ,1,k))
      PM (tr2,tr1,CC(ido-1,0,k),CC(ido-1,2,k))
      CH(ido-1,k,0)=tr2+tr2;
      CH(ido-1,k,1)=sqrt2*(tr1-ti1);
      CH(ido-1,k,2)=ti2+ti2;
      CH(ido-1,k,3)=-sqrt2*(tr1+ti1);
      }
  if (ido<=2) return;
  for (k=0; k<l1;++k)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      PM (tr2,tr1,CC(i-1,0,k),CC(ic-1,3,k))
      PM (ti1,ti2,CC(i  ,0,k),CC(ic  ,3,k))
      PM (tr4,ti3,CC(i  ,2,k),CC(ic  ,1,k))
      PM (tr3,ti4,CC(i-1,2,k),CC(ic-1,1,k))
      PM (CH(i-1,k,0),cr3,tr2,tr3)
      PM (CH(i  ,k,0),ci3,ti2,ti3)
      PM (cr4,cr2,tr1,tr4)
      PM (ci2,ci4,ti1,ti4)
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ci2,cr2)
      MULPM (CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),ci3,cr3)
      MULPM (CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),ci4,cr4)
      }
  }

static void radb5(size_t ido, size_t l1, const double *cc, double *ch,
  const double *wa)
  {
  const size_t cdim=5;
  static const double tr11= 0.3090169943749474241, ti11=0.95105651629515357212,
                      tr12=-0.8090169943749474241, ti12=0.58778525229247312917;
  size_t i, k, ic;
  double ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4,
         ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  for (k=0; k<l1; k++)
    {
    ti5=2*CC(0,2,k);
    ti4=2*CC(0,4,k);
    tr2=2*CC(ido-1,1,k);
    tr3=2*CC(ido-1,3,k);
    CH(0,k,0)=CC(0,0,k)+tr2+tr3;
    cr2=CC(0,0,k)+tr11*tr2+tr12*tr3;
    cr3=CC(0,0,k)+tr12*tr2+tr11*tr3;
    MULPM(ci5,ci4,ti5,ti4,ti11,ti12)
    PM(CH(0,k,4),CH(0,k,1),cr2,ci5)
    PM(CH(0,k,3),CH(0,k,2),cr3,ci4)
    }
  if (ido==1) return;
  for (k=0; k<l1;++k)
    for (i=2; i<ido; i+=2)
      {
      ic=ido-i;
      PM(tr2,tr5,CC(i-1,2,k),CC(ic-1,1,k))
      PM(ti5,ti2,CC(i  ,2,k),CC(ic  ,1,k))
      PM(tr3,tr4,CC(i-1,4,k),CC(ic-1,3,k))
      PM(ti4,ti3,CC(i  ,4,k),CC(ic  ,3,k))
      CH(i-1,k,0)=CC(i-1,0,k)+tr2+tr3;
      CH(i  ,k,0)=CC(i  ,0,k)+ti2+ti3;
      cr2=CC(i-1,0,k)+tr11*tr2+tr12*tr3;
      ci2=CC(i  ,0,k)+tr11*ti2+tr12*ti3;
      cr3=CC(i-1,0,k)+tr12*tr2+tr11*tr3;
      ci3=CC(i  ,0,k)+tr12*ti2+tr11*ti3;
      MULPM(cr5,cr4,tr5,tr4,ti11,ti12)
      MULPM(ci5,ci4,ti5,ti4,ti11,ti12)
      PM(dr4,dr3,cr3,ci4)
      PM(di3,di4,ci3,cr4)
      PM(dr5,dr2,cr2,ci5)
      PM(di2,di5,ci2,cr5)
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2)
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3)
      MULPM(CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),di4,dr4)
      MULPM(CH(i,k,4),CH(i-1,k,4),WA(3,i-2),WA(3,i-1),di5,dr5)
      }
  }

static void radbg(size_t ido, size_t ip, size_t l1, size_t idl1,
  double *cc, double *ch, const double *wa)
  {
  const size_t cdim=ip;
  static const double twopi=6.28318530717958647692;
  size_t idij, ipph, i, j, k, l, j2, ic, jc, lc, ik;
  double ai1, ai2, ar1, ar2, arg;
  double *csarr;
  size_t aidx;

  ipph=(ip+1)/ 2;
  for(k=0; k<l1; k++)
    memcpy(&CH(0,k,0),&CC(0,0,k),ido*sizeof(double));
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      {
      CH(0,k,j )=2*CC(ido-1,j2-1,k);
      CH(0,k,jc)=2*CC(0    ,j2  ,k);
      }
    }

  if(ido!=1)
    for(j=1,jc=ip-1; j<ipph; j++,jc--)
      for(k=0; k<l1; k++)
        for(i=2; i<ido; i+=2)
          {
          ic=ido-i;
          PM (CH(i-1,k,j ),CH(i-1,k,jc),CC(i-1,2*j,k),CC(ic-1,2*j-1,k))
          PM (CH(i  ,k,jc),CH(i  ,k,j ),CC(i  ,2*j,k),CC(ic  ,2*j-1,k))
          }

  csarr=RALLOC(double,2*ip);
  arg=twopi/ip;
  csarr[0]=1.;
  csarr[1]=0.;
  csarr[2]=csarr[2*ip-2]=cos(arg);
  csarr[3]=sin(arg); csarr[2*ip-1]=-csarr[3];
  for (i=2; i<=ip/2; ++i)
    {
    csarr[2*i]=csarr[2*ip-2*i]=cos(i*arg);
    csarr[2*i+1]=sin(i*arg);
    csarr[2*ip-2*i+1]=-csarr[2*i+1];
    }
  for(l=1; l<ipph; l++)
    {
    lc=ip-l;
    ar1=csarr[2*l];
    ai1=csarr[2*l+1];
    for(ik=0; ik<idl1; ik++)
      {
      C2(ik,l)=CH2(ik,0)+ar1*CH2(ik,1);
      C2(ik,lc)=ai1*CH2(ik,ip-1);
      }
    aidx=2*l;
    for(j=2; j<ipph; j++)
      {
      jc=ip-j;
      aidx+=2*l;
      if (aidx>=2*ip) aidx-=2*ip;
      ar2=csarr[aidx];
      ai2=csarr[aidx+1];
      for(ik=0; ik<idl1; ik++)
        {
        C2(ik,l )+=ar2*CH2(ik,j );
        C2(ik,lc)+=ai2*CH2(ik,jc);
        }
      }
    }
  DEALLOC(csarr);

  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      CH2(ik,0)+=CH2(ik,j);

  for(j=1,jc=ip-1; j<ipph; j++,jc--)
    for(k=0; k<l1; k++)
      PM (CH(0,k,jc),CH(0,k,j),C1(0,k,j),C1(0,k,jc))

  if(ido==1)
    return;
  for(j=1,jc=ip-1; j<ipph; j++,jc--)
    for(k=0; k<l1; k++)
      for(i=2; i<ido; i+=2)
        {
        PM (CH(i-1,k,jc),CH(i-1,k,j ),C1(i-1,k,j),C1(i  ,k,jc))
        PM (CH(i  ,k,j ),CH(i  ,k,jc),C1(i  ,k,j),C1(i-1,k,jc))
        }
  memcpy(cc,ch,idl1*sizeof(double));

  for(j=1; j<ip; j++)
    for(k=0; k<l1; k++)
      {
      C1(0,k,j)=CH(0,k,j);
      idij=(j-1)*ido+1;
      for(i=2; i<ido; i+=2,idij+=2)
        MULPM (C1(i,k,j),C1(i-1,k,j),wa[idij-1],wa[idij],CH(i,k,j),CH(i-1,k,j))
      }
  }

#undef CC
#undef CH
#undef PM
#undef MULPM


/*----------------------------------------------------------------------
   cfftf1, cfftb1, cfftf, cfftb, cffti1, cffti. Complex FFTs.
  ----------------------------------------------------------------------*/

static void cfft1(size_t n, cmplx c[], cmplx ch[], const cmplx wa[],
  const size_t ifac[], int isign)
  {
  size_t k1, l1=1, nf=ifac[1], iw=0;
  cmplx *p1=c, *p2=ch;

  for(k1=0; k1<nf; k1++)
    {
    size_t ip=ifac[k1+2];
    size_t l2=ip*l1;
    size_t ido = n/l2;
    if(ip==4)
      (isign>0) ? passb4(ido, l1, p1, p2, wa+iw)
                : passf4(ido, l1, p1, p2, wa+iw);
    else if(ip==2)
      (isign>0) ? passb2(ido, l1, p1, p2, wa+iw)
                : passf2(ido, l1, p1, p2, wa+iw);
    else if(ip==3)
      (isign>0) ? passb3(ido, l1, p1, p2, wa+iw)
                : passf3(ido, l1, p1, p2, wa+iw);
    else if(ip==5)
      (isign>0) ? passb5(ido, l1, p1, p2, wa+iw)
                : passf5(ido, l1, p1, p2, wa+iw);
    else if(ip==6)
      (isign>0) ? passb6(ido, l1, p1, p2, wa+iw)
                : passf6(ido, l1, p1, p2, wa+iw);
    else
      (isign>0) ? passbg(ido, ip, l1, p1, p2, wa+iw)
                : passfg(ido, ip, l1, p1, p2, wa+iw);
    SWAP(p1,p2,cmplx *);
    l1=l2;
    iw+=(ip-1)*ido;
    }
  if (p1!=c)
    memcpy (c,p1,n*sizeof(cmplx));
  }

void cfftf(size_t n, double c[], double wsave[])
  {
  if (n!=1)
    cfft1(n, (cmplx*)c, (cmplx*)wsave, (cmplx*)(wsave+2*n),
          (size_t*)(wsave+4*n),-1);
  }

void cfftb(size_t n, double c[], double wsave[])
  {
  if (n!=1)
    cfft1(n, (cmplx*)c, (cmplx*)wsave, (cmplx*)(wsave+2*n),
          (size_t*)(wsave+4*n),+1);
  }

static void factorize (size_t n, const size_t *pf, size_t npf, size_t *ifac)
  {
  size_t nl=n, nf=0, ntry=0, j=0, i;

startloop:
  j++;
  ntry = (j<=npf) ? pf[j-1] : ntry+2;
  do
    {
    size_t nq=nl / ntry;
    size_t nr=nl-ntry*nq;
    if (nr!=0)
      goto startloop;
    nf++;
    ifac[nf+1]=ntry;
    nl=nq;
    if ((ntry==2) && (nf!=1))
      {
      for (i=nf+1; i>2; --i)
        ifac[i]=ifac[i-1];
      ifac[2]=2;
      }
    }
  while(nl!=1);
  ifac[0]=n;
  ifac[1]=nf;
  }

static void cffti1(size_t n, double wa[], size_t ifac[])
  {
  static const size_t ntryh[5]={4,6,3,2,5};
  static const double twopi=6.28318530717958647692;
  size_t j, k, fi;

  double argh=twopi/n;
  size_t i=0, l1=1;
  factorize (n,ntryh,5,ifac);
  for(k=1; k<=ifac[1]; k++)
    {
    size_t ip=ifac[k+1];
    size_t ido=n/(l1*ip);
    for(j=1; j<ip; j++)
      {
      size_t is = i;
      double argld=j*l1*argh;
      wa[i  ]=1;
      wa[i+1]=0;
      for(fi=1; fi<=ido; fi++)
        {
        double arg=fi*argld;
        i+=2;
        wa[i  ]=cos(arg);
        wa[i+1]=sin(arg);
        }
      if(ip>6)
        {
        wa[is  ]=wa[i  ];
        wa[is+1]=wa[i+1];
        }
      }
    l1*=ip;
    }
  }

void cffti(size_t n, double wsave[])
  { if (n!=1) cffti1(n, wsave+2*n,(size_t*)(wsave+4*n)); }


/*----------------------------------------------------------------------
   rfftf1, rfftb1, rfftf, rfftb, rffti1, rffti. Real FFTs.
  ----------------------------------------------------------------------*/

static void rfftf1(size_t n, double c[], double ch[], const double wa[],
  const size_t ifac[])
  {
  size_t k1, l1=n, nf=ifac[1], iw=n-1;
  double *p1=ch, *p2=c;

  for(k1=1; k1<=nf;++k1)
    {
    size_t ip=ifac[nf-k1+2];
    size_t ido=n / l1;
    l1 /= ip;
    iw-=(ip-1)*ido;
    SWAP (p1,p2,double *);
    if(ip==4)
      radf4(ido, l1, p1, p2, wa+iw);
    else if(ip==2)
      radf2(ido, l1, p1, p2, wa+iw);
    else if(ip==3)
      radf3(ido, l1, p1, p2, wa+iw);
    else if(ip==5)
      radf5(ido, l1, p1, p2, wa+iw);
    else
      {
      if (ido==1)
        SWAP (p1,p2,double *);
      radfg(ido, ip, l1, ido*l1, p1, p2, wa+iw);
      SWAP (p1,p2,double *);
      }
    }
  if (p1==c)
    memcpy (c,ch,n*sizeof(double));
  }

static void rfftb1(size_t n, double c[], double ch[], const double wa[],
  const size_t ifac[])
  {
  size_t k1, l1=1, nf=ifac[1], iw=0;
  double *p1=c, *p2=ch;

  for(k1=1; k1<=nf; k1++)
    {
    size_t ip = ifac[k1+1],
           ido= n/(ip*l1);
    if(ip==4)
      radb4(ido, l1, p1, p2, wa+iw);
    else if(ip==2)
      radb2(ido, l1, p1, p2, wa+iw);
    else if(ip==3)
      radb3(ido, l1, p1, p2, wa+iw);
    else if(ip==5)
      radb5(ido, l1, p1, p2, wa+iw);
    else
      {
      radbg(ido, ip, l1, ido*l1, p1, p2, wa+iw);
      if (ido!=1)
        SWAP (p1,p2,double *);
      }
    SWAP (p1,p2,double *);
    l1*=ip;
    iw+=(ip-1)*ido;
    }
  if (p1!=c)
    memcpy (c,ch,n*sizeof(double));
  }

void rfftf(size_t n, double r[], double wsave[])
  { if(n!=1) rfftf1(n, r, wsave, wsave+n,(size_t*)(wsave+2*n)); }

void rfftb(size_t n, double r[], double wsave[])
  { if(n!=1) rfftb1(n, r, wsave, wsave+n,(size_t*)(wsave+2*n)); }

static void rffti1(size_t n, double wa[], size_t ifac[])
  {
  static const size_t ntryh[4]={4,2,3,5};
  static const double twopi=6.28318530717958647692;
  size_t i, j, k, fi;

  double argh=twopi/n;
  size_t is=0, l1=1;
  factorize (n,ntryh,4,ifac);
  for (k=1; k<ifac[1]; k++)
    {
    size_t ip=ifac[k+1],
           ido=n/(l1*ip);
    for (j=1; j<ip; ++j)
      {
      double argld=j*l1*argh;
      for(i=is,fi=1; i<=ido+is-3; i+=2,++fi)
        {
        double arg=fi*argld;
        wa[i  ]=cos(arg);
        wa[i+1]=sin(arg);
        }
      is+=ido;
      }
    l1*=ip;
    }
  }

void rffti(size_t n, double wsave[])
  { if (n!=1) rffti1(n, wsave+n,(size_t*)(wsave+2*n)); }

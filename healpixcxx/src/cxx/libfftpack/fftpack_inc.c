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

#ifdef BACKWARD
#define PSIGN +
#define PMSIGNC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
/* a = b*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
#else
#define PSIGN -
#define PMSIGNC(a,b,c,d) { a.r=c.r-d.r; a.i=c.i-d.i; b.r=c.r+d.r; b.i=c.i+d.i; }
/* a = conj(b)*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r+b.i*c.i; a.i=b.r*c.i-b.i*c.r; }
#endif

static void X(2) (size_t ido, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=2;
  size_t k,i;
  cmplx t;
  if (ido==1)
    for (k=0;k<l1;++k)
      PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
  else
    for (k=0;k<l1;++k)
      for (i=0;i<ido;++i)
        {
        PMC (CH(i,k,0),t,CC(i,0,k),CC(i,1,k))
        MULPMSIGNC (CH(i,k,1),WA(0,i),t)
        }
  }

static void X(3)(size_t ido, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=3;
  static const double taur=-0.5, taui= PSIGN 0.86602540378443864676;
  size_t i, k;
  cmplx c2, c3, d2, d3, t2;

  if (ido==1)
    for (k=0; k<l1; ++k)
      {
      PMC (t2,c3,CC(0,1,k),CC(0,2,k))
      ADDC (CH(0,k,0),t2,CC(0,0,k))
      SCALEC(t2,taur)
      ADDC(c2,CC(0,0,k),t2)
      SCALEC(c3,taui)
      CONJFLIPC(c3)
      PMC(CH(0,k,1),CH(0,k,2),c2,c3)
      }
  else
    for (k=0; k<l1; ++k)
      for (i=0; i<ido; ++i)
        {
        PMC (t2,c3,CC(i,1,k),CC(i,2,k))
        ADDC (CH(i,k,0),t2,CC(i,0,k))
        SCALEC(t2,taur)
        ADDC(c2,CC(i,0,k),t2)
        SCALEC(c3,taui)
        CONJFLIPC(c3)
        PMC(d2,d3,c2,c3)
        MULPMSIGNC(CH(i,k,1),WA(0,i),d2)
        MULPMSIGNC(CH(i,k,2),WA(1,i),d3)
        }
  }

static void X(4)(size_t ido, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=4;
  size_t i, k;
  cmplx c2, c3, c4, t1, t2, t3, t4;

  if (ido==1)
    for (k=0; k<l1; ++k)
      {
      PMC(t2,t1,CC(0,0,k),CC(0,2,k))
      PMC(t3,t4,CC(0,1,k),CC(0,3,k))
      CONJFLIPC(t4)
      PMC(CH(0,k,0),CH(0,k,2),t2,t3)
      PMSIGNC (CH(0,k,1),CH(0,k,3),t1,t4)
      }
  else
    for (k=0; k<l1; ++k)
      for (i=0; i<ido; ++i)
        {
        PMC(t2,t1,CC(i,0,k),CC(i,2,k))
        PMC(t3,t4,CC(i,1,k),CC(i,3,k))
        CONJFLIPC(t4)
        PMC(CH(i,k,0),c3,t2,t3)
        PMSIGNC (c2,c4,t1,t4)
        MULPMSIGNC (CH(i,k,1),WA(0,i),c2)
        MULPMSIGNC (CH(i,k,2),WA(1,i),c3)
        MULPMSIGNC (CH(i,k,3),WA(2,i),c4)
        }
  }

static void X(5)(size_t ido, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=5;
  static const double tr11= 0.3090169943749474241,
                      ti11= PSIGN 0.95105651629515357212,
                      tr12=-0.8090169943749474241,
                      ti12= PSIGN 0.58778525229247312917;
  size_t i, k;
  cmplx c2, c3, c4, c5, d2, d3, d4, d5, t2, t3, t4, t5;

  if (ido==1)
    for (k=0; k<l1; ++k)
      {
      PMC (t2,t5,CC(0,1,k),CC(0,4,k))
      PMC (t3,t4,CC(0,2,k),CC(0,3,k))
      CH(0,k,0).r=CC(0,0,k).r+t2.r+t3.r;
      CH(0,k,0).i=CC(0,0,k).i+t2.i+t3.i;
      c2.r=CC(0,0,k).r+tr11*t2.r+tr12*t3.r;
      c2.i=CC(0,0,k).i+tr11*t2.i+tr12*t3.i;
      c3.r=CC(0,0,k).r+tr12*t2.r+tr11*t3.r;
      c3.i=CC(0,0,k).i+tr12*t2.i+tr11*t3.i;
      c5.r=ti11*t5.r+ti12*t4.r;
      c5.i=ti11*t5.i+ti12*t4.i;
      c4.r=ti12*t5.r-ti11*t4.r;
      c4.i=ti12*t5.i-ti11*t4.i;
      CONJFLIPC(c5)
      PMC(CH(0,k,1),CH(0,k,4),c2,c5)
      CONJFLIPC(c4)
      PMC(CH(0,k,2),CH(0,k,3),c3,c4)
      }
  else
    for (k=0; k<l1; ++k)
      for (i=0; i<ido; ++i)
        {
        PMC (t2,t5,CC(i,1,k),CC(i,4,k))
        PMC (t3,t4,CC(i,2,k),CC(i,3,k))
        CH(i,k,0).r=CC(i,0,k).r+t2.r+t3.r;
        CH(i,k,0).i=CC(i,0,k).i+t2.i+t3.i;
        c2.r=CC(i,0,k).r+tr11*t2.r+tr12*t3.r;
        c2.i=CC(i,0,k).i+tr11*t2.i+tr12*t3.i;
        c3.r=CC(i,0,k).r+tr12*t2.r+tr11*t3.r;
        c3.i=CC(i,0,k).i+tr12*t2.i+tr11*t3.i;
        c5.r=ti11*t5.r+ti12*t4.r;
        c5.i=ti11*t5.i+ti12*t4.i;
        c4.r=ti12*t5.r-ti11*t4.r;
        c4.i=ti12*t5.i-ti11*t4.i;
        CONJFLIPC(c5)
        PMC(d2,d5,c2,c5)
        CONJFLIPC(c4)
        PMC(d3,d4,c3,c4)
        MULPMSIGNC (CH(i,k,1),WA(0,i),d2)
        MULPMSIGNC (CH(i,k,2),WA(1,i),d3)
        MULPMSIGNC (CH(i,k,3),WA(2,i),d4)
        MULPMSIGNC (CH(i,k,4),WA(3,i),d5)
        }
  }

static void X(6)(size_t ido, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=6;
  static const double taui= PSIGN 0.86602540378443864676;
  cmplx ta1,ta2,ta3,a0,a1,a2,tb1,tb2,tb3,b0,b1,b2,d1,d2,d3,d4,d5;
  size_t i, k;

  if (ido==1)
    for (k=0; k<l1; ++k)
      {
      PMC(ta1,ta3,CC(0,2,k),CC(0,4,k))
      ta2.r = CC(0,0,k).r - .5*ta1.r;
      ta2.i = CC(0,0,k).i - .5*ta1.i;
      SCALEC(ta3,taui)
      ADDC(a0,CC(0,0,k),ta1)
      CONJFLIPC(ta3)
      PMC(a1,a2,ta2,ta3)
      PMC(tb1,tb3,CC(0,5,k),CC(0,1,k))
      tb2.r = CC(0,3,k).r - .5*tb1.r;
      tb2.i = CC(0,3,k).i - .5*tb1.i;
      SCALEC(tb3,taui)
      ADDC(b0,CC(0,3,k),tb1)
      CONJFLIPC(tb3)
      PMC(b1,b2,tb2,tb3)
      PMC(CH(0,k,0),CH(0,k,3),a0,b0)
      PMC(CH(0,k,4),CH(0,k,1),a1,b1)
      PMC(CH(0,k,2),CH(0,k,5),a2,b2)
      }
  else
    for (k=0; k<l1; ++k)
      for (i=0; i<ido; ++i)
        {
        PMC(ta1,ta3,CC(i,2,k),CC(i,4,k))
        ta2.r = CC(i,0,k).r - .5*ta1.r;
        ta2.i = CC(i,0,k).i - .5*ta1.i;
        SCALEC(ta3,taui)
        ADDC(a0,CC(i,0,k),ta1)
        CONJFLIPC(ta3)
        PMC(a1,a2,ta2,ta3)
        PMC(tb1,tb3,CC(i,5,k),CC(i,1,k))
        tb2.r = CC(i,3,k).r - .5*tb1.r;
        tb2.i = CC(i,3,k).i - .5*tb1.i;
        SCALEC(tb3,taui)
        ADDC(b0,CC(i,3,k),tb1)
        CONJFLIPC(tb3)
        PMC(b1,b2,tb2,tb3)
        PMC(CH(i,k,0),d3,a0,b0)
        PMC(d4,d1,a1,b1)
        PMC(d2,d5,a2,b2)
        MULPMSIGNC (CH(i,k,1),WA(0,i),d1)
        MULPMSIGNC (CH(i,k,2),WA(1,i),d2)
        MULPMSIGNC (CH(i,k,3),WA(2,i),d3)
        MULPMSIGNC (CH(i,k,4),WA(3,i),d4)
        MULPMSIGNC (CH(i,k,5),WA(4,i),d5)
        }
  }

static void X(g)(size_t ido, size_t ip, size_t l1, const cmplx *cc, cmplx *ch,
  const cmplx *wa)
  {
  const size_t cdim=ip;
  cmplx *tarr=RALLOC(cmplx,2*ip);
  cmplx *ccl=tarr, *wal=tarr+ip;
  size_t i,j,k,l,jc,lc;
  size_t ipph = (ip+1)/2;

  for (i=1; i<ip; ++i)
    wal[i]=wa[ido*(i-1)];
  for (k=0; k<l1; ++k)
    for (i=0; i<ido; ++i)
      {
      cmplx s=CC(i,0,k);
      ccl[0] = CC(i,0,k);
      for(j=1,jc=ip-1; j<ipph; ++j,--jc)
        {
        PMC (ccl[j],ccl[jc],CC(i,j,k),CC(i,jc,k))
        ADDC (s,s,ccl[j])
        }
      CH(i,k,0) = s;
      for (j=1, jc=ip-1; j<=ipph; ++j,--jc)
        {
        cmplx abr=ccl[0], abi={0.,0.};
        size_t iang=0;
        for (l=1,lc=ip-1; l<ipph; ++l,--lc)
          {
          iang+=j;
          if (iang>ip) iang-=ip;
          abr.r += ccl[l ].r*wal[iang].r;
          abr.i += ccl[l ].i*wal[iang].r;
          abi.r += ccl[lc].r*wal[iang].i;
          abi.i += ccl[lc].i*wal[iang].i;
          }
#ifndef BACKWARD
          { abi.i=-abi.i; abi.r=-abi.r; }
#endif
        CONJFLIPC(abi)
        PMC(CH(i,k,j),CH(i,k,jc),abr,abi)
        }
      }

  DEALLOC(tarr);

  if (ido==1) return;

  for (j=1; j<ip; ++j)
    for (k=0; k<l1; ++k)
      {
      size_t idij=(j-1)*ido+1;
      for(i=1; i<ido; ++i, ++idij)
        {
        cmplx t=CH(i,k,j);
        MULPMSIGNC (CH(i,k,j),wa[idij],t)
        }
      }
  }

#undef PSIGN
#undef PMSIGNC
#undef MULPMSIGNC

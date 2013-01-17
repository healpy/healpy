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

/*! \file sharp_core.c
 *  Computational core
 *
 *  Copyright (C) 2012-2013 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <complex.h>
#include <math.h>
#include <string.h>
#include "sharp_vecsupport.h"
#include "sharp_complex_hacks.h"
#include "sharp_ylmgen_c.h"
#include "sharp.h"
#include "sharp_core.h"
#include "c_utils.h"

typedef complex double dcmplx;

// must be in the range [0;6]
#define MAXJOB_SPECIAL 2

#define XCONCAT2(a,b) a##_##b
#define CONCAT2(a,b) XCONCAT2(a,b)
#define XCONCAT3(a,b,c) a##_##b##_##c
#define CONCAT3(a,b,c) XCONCAT3(a,b,c)

#define nvec 1
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 2
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 3
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 4
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 5
#include "sharp_core_inchelper.c"
#undef nvec

#define nvec 6
#include "sharp_core_inchelper.c"
#undef nvec

void inner_loop (sharp_job *job, const int *ispair,const double *cth,
  const double *sth, int llim, int ulim, sharp_Ylmgen_C *gen, int mi,
  const int *mlim)
  {
  int njobs=job->ntrans, nv=job->flags&SHARP_NVMAX;
  if (njobs<=MAXJOB_SPECIAL)
    {
    switch (njobs*16+nv)
      {
#if ((MAXJOB_SPECIAL>=1)&&(SHARP_MAXTRANS>=1))
      case 0x11:
        CONCAT3(inner_loop,1,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x12:
        CONCAT3(inner_loop,2,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x13:
        CONCAT3(inner_loop,3,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x14:
        CONCAT3(inner_loop,4,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x15:
        CONCAT3(inner_loop,5,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x16:
        CONCAT3(inner_loop,6,1) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
#if ((MAXJOB_SPECIAL>=2)&&(SHARP_MAXTRANS>=2))
      case 0x21:
        CONCAT3(inner_loop,1,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x22:
        CONCAT3(inner_loop,2,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x23:
        CONCAT3(inner_loop,3,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x24:
        CONCAT3(inner_loop,4,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x25:
        CONCAT3(inner_loop,5,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x26:
        CONCAT3(inner_loop,6,2) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
#if ((MAXJOB_SPECIAL>=3)&&(SHARP_MAXTRANS>=3))
      case 0x31:
        CONCAT3(inner_loop,1,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x32:
        CONCAT3(inner_loop,2,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x33:
        CONCAT3(inner_loop,3,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x34:
        CONCAT3(inner_loop,4,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x35:
        CONCAT3(inner_loop,5,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x36:
        CONCAT3(inner_loop,6,3) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
#if ((MAXJOB_SPECIAL>=4)&&(SHARP_MAXTRANS>=4))
      case 0x41:
        CONCAT3(inner_loop,1,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x42:
        CONCAT3(inner_loop,2,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x43:
        CONCAT3(inner_loop,3,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x44:
        CONCAT3(inner_loop,4,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x45:
        CONCAT3(inner_loop,5,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x46:
        CONCAT3(inner_loop,6,4) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
#if ((MAXJOB_SPECIAL>=5)&&(SHARP_MAXTRANS>=5))
      case 0x51:
        CONCAT3(inner_loop,1,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x52:
        CONCAT3(inner_loop,2,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x53:
        CONCAT3(inner_loop,3,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x54:
        CONCAT3(inner_loop,4,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x55:
        CONCAT3(inner_loop,5,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x56:
        CONCAT3(inner_loop,6,5) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
#if ((MAXJOB_SPECIAL>=6)&&(SHARP_MAXTRANS>=6))
      case 0x61:
        CONCAT3(inner_loop,1,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x62:
        CONCAT3(inner_loop,2,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x63:
        CONCAT3(inner_loop,3,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x64:
        CONCAT3(inner_loop,4,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x65:
        CONCAT3(inner_loop,5,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
      case 0x66:
        CONCAT3(inner_loop,6,6) (job, ispair,cth,sth,llim,ulim,gen,mi,mlim);
        return;
#endif
      }
    }
#if (SHARP_MAXTRANS>MAXJOB_SPECIAL)
  else
    {
    switch (nv)
      {
      case 1:
        CONCAT2(inner_loop,1)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      case 2:
        CONCAT2(inner_loop,2)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      case 3:
        CONCAT2(inner_loop,3)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      case 4:
        CONCAT2(inner_loop,4)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      case 5:
        CONCAT2(inner_loop,5)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      case 6:
        CONCAT2(inner_loop,6)
          (job, ispair,cth,sth,llim,ulim,gen,mi,mlim,job->ntrans);
        return;
      }
    }
#endif
  UTIL_FAIL("Incorrect vector parameters");
  }

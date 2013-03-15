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

/*! \file sharp_ylmgen_c.h
 *  Code for efficient calculation of Y_lm(phi=0,theta)
 *
 *  Copyright (C) 2005-2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef SHARP_YLMGEN_C_H
#define SHARP_YLMGEN_C_H

#ifdef __cplusplus
extern "C" {
#endif

enum { sharp_minscale=0, sharp_limscale=1, sharp_maxscale=1 };
static const double sharp_fbig=0x1p+800,sharp_fsmall=0x1p-800;
static const double sharp_ftol=0x1p-60;
static const double sharp_fbighalf=0x1p+400;

typedef struct { double f[2]; } sharp_ylmgen_dbl2;
typedef struct { double f[3]; } sharp_ylmgen_dbl3;

typedef struct
  {
/* for public use; immutable during lifetime */
  int lmax, mmax, s;
  double *cf;

/* for public use; will typically change after call to Ylmgen_prepare() */
  int m;

/* used if s==0 */
  double *mfac;
  sharp_ylmgen_dbl2 *rf;

/* used if s!=0 */
  int sinPow, cosPow, preMinus_p, preMinus_m;
  double *prefac;
  int *fscale;
  sharp_ylmgen_dbl3 *fx;

/* internal usage only */
/* used if s==0 */
  double *root, *iroot;

/* used if s!=0 */
  double *flm1, *flm2, *inv;
  int mlo, mhi;
  } sharp_Ylmgen_C;

/*! Creates a generator which will calculate helper data for Y_lm calculation
    up to \a l=l_max and \a m=m_max. */
void sharp_Ylmgen_init (sharp_Ylmgen_C *gen, int l_max, int m_max, int spin);

/*! Deallocates a generator previously initialised by Ylmgen_init(). */
void sharp_Ylmgen_destroy (sharp_Ylmgen_C *gen);

/*! Prepares the object for the calculation at \a m. */
void sharp_Ylmgen_prepare (sharp_Ylmgen_C *gen, int m);

/*! Returns a pointer to an array with \a lmax+1 entries containing
    normalisation factors that must be applied to Y_lm values computed for
    \a spin. The array must be deallocated (using free()) by the user. */
double *sharp_Ylmgen_get_norm (int lmax, int spin);

/*! Returns a pointer to an array with \a lmax+1 entries containing
    normalisation factors that must be applied to Y_lm values computed for
    first derivatives. The array must be deallocated (using free()) by the
    user. */
double *sharp_Ylmgen_get_d1norm (int lmax);

#ifdef __cplusplus
}
#endif

#endif

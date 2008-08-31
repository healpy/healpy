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

/*! \file ls_fft.h
 *  Interface for the LevelS FFT package.
 *
 *  Copyright (C) 2004 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_LS_FFT_H
#define PLANCK_LS_FFT_H

#ifdef __cplusplus
extern "C" {
#endif

/*!\defgroup fftgroup FFT interface
This package is intended to calculate one-dimensional real or complex FFTs
with high accuracy and good efficiency even for lengths containing large
prime factors.
The code is written in C, but a Fortran wrapper exists as well.

Before any FFT is executed, a plan must be generated for it. Plan creation
is designed to be fast, so that there is no significant overhead if the
plan is only used once or a few times.

The main component of the code is a C port of Swarztrauber's FFTPACK
(http://www.netlib.org/fftpack/), which was originally done by Pekka Janhunen
and reformatted by Joerg Arndt.

I added a few digits to the floating-point constants to achieve higher
precision, split the complex transforms into separate routines for forward
and backward transforms (to increase performance a bit), and replaced
the iterative twiddling factor calculation in radfg() and radbg() by an exact
calculation of the factors.

Since FFTPACK becomes quite slow for FFT lengths with large prime factors
(in the worst case of prime lengths it reaches O(n*n) complexity), I
implemented Bluestein's algorithm, which computes a FFT of length n by
several FFTs of length n2>=2*n and a convolution. Since n2 can be chosen
to be highly composite, this algorithm is more efficient if n has large
prime factors. The longer FFTs themselves are then computed using the FFTPACK
routines.
Bluestein's algorithm was implemented according to the description at
<a href="http://en.wikipedia.org/wiki/Bluestein%27s_FFT_algorithm">Wikipedia</a>.

\b Thread-safety:
All routines can be called concurrently; all information needed by ls_fft
is stored in the plan variable. However, using the same plan variable on
multiple threads simultaneously is not supported and will lead to data
corruption.
*/
/*! \{ */

typedef struct
  {
  double *work;
  int length;
  int bluestein;
  } complex_plan_i;

/*! The opaque handle type for complex-FFT plans. */
typedef complex_plan_i * complex_plan;

/*! Returns a plan for a complex FFT with \a length elements. */
complex_plan make_complex_plan (int length);
/*! Destroys a plan for a complex FFT. */
void kill_complex_plan (complex_plan plan);
/*! Computes a complex forward FFT on \a data, using \a plan.
    \a Data has the form r0, i0, r1, i1, ..., r[length-1], i[length-1]. */
void complex_plan_forward (complex_plan plan, double *data);
/*! Computes a complex backward FFT on \a data, using \a plan.
    \a Data has the form r0, i0, r1, i1, ..., r[length-1], i[length-1]. */
void complex_plan_backward (complex_plan plan, double *data);

typedef struct
  {
  double *work;
  int length;
  int bluestein;
  } real_plan_i;

/*! The opaque handle type for real-FFT plans. */
typedef real_plan_i * real_plan;

/*! Returns a plan for a real FFT with \a length elements. */
real_plan make_real_plan (int length);
/*! Destroys a plan for a real FFT. */
void kill_real_plan (real_plan plan);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a data has the form r0, r1, ..., r[length-1];
    - on exit, it has the form r0, r1, i1, r2, i2, ...
      (a total of \a length values). */
void real_plan_forward_fftpack (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTPACK storage scheme:
    - on entry, \a data has the form r0, r1, i1, r2, i2, ...
    (a total of \a length values);
    - on exit, it has the form r0, r1, ..., r[length-1]. */
void real_plan_backward_fftpack (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming the FFTW halfcomplex storage scheme:
    - on entry, \a data has the form r0, r1, ..., r[length-1];
    - on exit, it has the form r0, r1, r2, ..., i2, i1. */
void real_plan_forward_fftw (real_plan plan, double *data);
/*! Computes a real backward FFT on \a data, using \a plan
    and assuming the FFTW halfcomplex storage scheme:
    - on entry, \a data has the form r0, r1, r2, ..., i2, i1.
    - on exit, it has the form r0, r1, ..., r[length-1]. */
void real_plan_backward_fftw (real_plan plan, double *data);
/*! Computes a real forward FFT on \a data, using \a plan
    and assuming a full-complex storage scheme:
    - on entry, \a data has the form r0, [ignored], r1, [ignored], ...,
      r[length-1], [ignored];
    - on exit, it has the form r0, i0, r1, i1, ..., r[length-1], i[length-1].
    */
void real_plan_forward_c (real_plan plan, double *data);
/*! Computes a real backward FFT on \a data, using \a plan
    and assuming a full-complex storage scheme:
    - on entry, \a data has the form r0, i0, r1, i1, ...,
      r[length-1], i[length-1];
    - on exit, it has the form r0, 0, r1, 0, ..., r[length-1], 0. */
void real_plan_backward_c (real_plan plan, double *data);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif

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
  fftpack.h : function declarations for fftpack.c
  Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber
  (Version 4, 1985).

  Pekka Janhunen 23.2.1995

  (reformatted by joerg arndt)

  reformatted and slightly enhanced by Martin Reinecke (2004)
 */

#ifndef PLANCK_FFTPACK_H
#define PLANCK_FFTPACK_H

#include "c_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! forward complex transform */
void cfftf(size_t N, double complex_data[], double wrk[]);
/*! backward complex transform */
void cfftb(size_t N, double complex_data[], double wrk[]);
/*! initializer for complex transforms */
void cffti(size_t N, double wrk[]);

/*! forward real transform */
void rfftf(size_t N, double data[], double wrk[]);
/*! backward real transform */
void rfftb(size_t N, double data[], double wrk[]);
/*! initializer for real transforms */
void rffti(size_t N, double wrk[]);

#ifdef __cplusplus
}
#endif

#endif

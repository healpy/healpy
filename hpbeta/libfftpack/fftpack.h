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

#ifndef PLANCK_FFTPACK_H
#define PLANCK_FFTPACK_H

#ifdef __cplusplus
extern "C" {
#endif

/*! forward complex transform */
void cfftf(int N, double complex_data[], double wrk[]);
/*! backward complex transform */
void cfftb(int N, double complex_data[], double wrk[]);
/*! initializer for complex transforms */
void cffti(int N, double wrk[]);

/*! forward real transform */
void rfftf(int N, double data[], double wrk[]);             
/*! backward real transform */
void rfftb(int N, double data[], double wrk[]);
/*! initializer for real transforms */
void rffti(int N, double wrk[]);

#ifdef __cplusplus
}
#endif

#endif

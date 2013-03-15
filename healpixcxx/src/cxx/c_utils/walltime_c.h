/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file walltime_c.h
 *  Functionality for reading wall clock time
 *
 *  Copyright (C) 2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_WALLTIME_C_H
#define PLANCK_WALLTIME_C_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Returns an approximation of the current wall time (in seconds).
    The first available of the following timers will be used:
    <ul>
    <li> \a omp_get_wtime(), if OpenMP is available
    <li> \a MPI_Wtime(), if MPI is available
    <li> \a gettimeofday() otherwise
    </ul>
    \note Only useful for measuring time differences.
    \note This function has an execution time between 10 and 100 nanoseconds. */
double wallTime(void);

#ifdef __cplusplus
}
#endif

#endif

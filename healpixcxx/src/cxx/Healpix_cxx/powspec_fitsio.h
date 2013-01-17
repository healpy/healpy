/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file powspec_fitsio.h
 *  Copyright (C) 2003-2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef POWSPEC_FITSIO_H
#define POWSPEC_FITSIO_H

#include <string>
class fitshandle;

class PowSpec;

/*! \defgroup powspec_fitsio_group FITS-based I/O of power spectra */
/*! \{ */

/*! Opens the FITS file \a filename, jumps to HDU 2, and reads \a nspecs
    columns into \a powspec. \a nspecs must be 1, 4, or 6. */
void read_powspec_from_fits (const std::string &infile,
  PowSpec &powspec, int nspecs, int lmax);
/*! Inserts a new binary table into \a out, which contains \a nspecs columns
    of FITS type TDOUBLE, and writes the components of \a powspec into it.
    \a nspecs must be 1, 4, or 6. */
void write_powspec_to_fits (fitshandle &out,
  const PowSpec &powspec, int nspecs);

/*! Creates a new FITS file called \a outfile, inserts a binary table,
    which contains \a nspecs columns of FITS type TDOUBLE, and writes the
    components of \a powspec into it. \a nspecs must be 1, 4, or 6. */
void write_powspec_to_fits (const std::string &outfile,
  const PowSpec &powspec, int nspecs);

/*! \} */

#endif

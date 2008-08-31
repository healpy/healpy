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

/*
 *  Copyright (C) 2003, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef HEALPIX_DATA_IO_H
#define HEALPIX_DATA_IO_H

#include <string>

class paramfile;
class simparams;
template<typename T> class arr;

void read_weight_ring (const std::string &dir, int nside,
  arr<double> &weight);
void get_ring_weights (paramfile &params, simparams &par, int nside,
  arr<double> &weight);

void read_pixwin (const std::string &dir, int nside, arr<double> &temp);
void read_pixwin (const std::string &dir, int nside, arr<double> &temp,
  arr<double> &pol);

void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin);
void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin, arr<double> &pixwin_pol);

#endif

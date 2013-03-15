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
 *  Copyright (C) 2005-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "alm_fitsio.h"
#include "fitshandle.h"
#include "alm_powspec_tools.h"
#include "trafos.h"
#include "announce.h"
#include "string_utils.h"
#include "lsconstants.h"

using namespace std;

namespace {

Trafo maketrafo (int num)
  {
  switch (num)
    {
    case  1: return Trafo(2000,2000,Equatorial,Galactic);
    case  2: return Trafo(2000,2000,Galactic,Equatorial);
    case  3: return Trafo(2000,2000,Equatorial,Ecliptic);
    case  4: return Trafo(2000,2000,Ecliptic,Equatorial);
    case  5: return Trafo(2000,2000,Ecliptic,Galactic);
    case  6: return Trafo(2000,2000,Galactic,Ecliptic);
    case  7: return Trafo(1950,1950,Equatorial,Galactic);
    case  8: return Trafo(1950,1950,Galactic,Equatorial);
    case  9: return Trafo(1950,1950,Equatorial,Ecliptic);
    case 10: return Trafo(1950,1950,Ecliptic,Equatorial);
    case 11: return Trafo(1950,1950,Ecliptic,Galactic);
    case 12: return Trafo(1950,1950,Galactic,Ecliptic);
    default: planck_fail("Unsupported transformation "+dataToString(num));
    }
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup("rotalm_cxx", (argc==5)||(argc==7),
    "Usage: rotalm_cxx <infile> <outfile> <itransform> <pol>\n"
    "or   : rotalm_cxx <infile> <outfile> <psi> <theta> <phi> <pol>\n\n"
    "itransform: 1: Equatorial (2000) -> Galactic   (2000)\n"
    "            2: Galactic   (2000) -> Equatorial (2000)\n"
    "            3: Equatorial (2000) -> Ecliptic   (2000)\n"
    "            4: Ecliptic   (2000) -> Equatorial (2000)\n"
    "            5: Ecliptic   (2000) -> Galactic   (2000)\n"
    "            6: Galactic   (2000) -> Ecliptic   (2000)\n"
    "            7: Equatorial (1950) -> Galactic   (1950)\n"
    "            8: Galactic   (1950) -> Equatorial (1950)\n"
    "            9: Equatorial (1950) -> Ecliptic   (1950)\n"
    "           10: Ecliptic   (1950) -> Equatorial (1950)\n"
    "           11: Ecliptic   (1950) -> Galactic   (1950)\n"
    "           12: Galactic   (1950) -> Ecliptic   (1950)\n\n"
    "psi, theta, phi: Euler angles (in degrees)\n\n"
    "pol: T or F\n");

  string infile  = argv[1];
  string outfile = argv[2];
  bool polarisation = stringToData<bool>(argv[argc-1]);
  rotmatrix rm;
  if (argc==5)
    {
    int trafo = stringToData<int>(argv[3]);
    Trafo tr(maketrafo(trafo));
    rm=tr.Matrix();
    }
  else
    {
    rm.Make_CPAC_Euler_Matrix(degr2rad*stringToData<double>(argv[5]),
                              degr2rad*stringToData<double>(argv[4]),
                              degr2rad*stringToData<double>(argv[3]));
    }

  Alm<xcomplex<double> > almT,almG,almC;

  if (!polarisation)
    {
    int lmax, dummy;
    get_almsize (infile,lmax,dummy);
    read_Alm_from_fits (infile, almT, lmax, lmax);
    rotate_alm(almT,rm);
    write_Alm_to_fits (outfile,almT,lmax,lmax,PLANCK_FLOAT32);
    }
  else
    {
    int lmax, dummy;
    get_almsize_pol (infile,lmax,dummy);
    read_Alm_from_fits (infile, almT, almG, almC, lmax, lmax);
    rotate_alm(almT,almG,almC,rm);
    write_Alm_to_fits (outfile,almT,almG,almC,lmax,lmax,PLANCK_FLOAT32);
    }

PLANCK_DIAGNOSIS_END
  }

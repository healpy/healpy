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
 *  Copyright (C) 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm.h"
#include "alm_fitsio.h"
#include "fitshandle.h"
#include "alm_powspec_tools.h"
#include "trafos.h"

using namespace std;

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
    default: throw Message_error ("Unsupported transformation");
    }
  }

int main(int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  announce ("rotalm_cxx");

  if (argc!=5)
    {
    cout << "Usage: rotalm_cxx <infile> <outfile> <itransform> <pol>"
         << endl
         << "Transform 1: Equatorial (2000) -> Galactic   (2000)" << endl
         << "          2: Galactic   (2000) -> Equatorial (2000)" << endl
         << "          3: Equatorial (2000) -> Ecliptic   (2000)" << endl
         << "          4: Ecliptic   (2000) -> Equatorial (2000)" << endl
         << "          5: Ecliptic   (2000) -> Galactic   (2000)" << endl
         << "          6: Galactic   (2000) -> Ecliptic   (2000)" << endl
         << "          7: Equatorial (1950) -> Galactic   (1950)" << endl
         << "          8: Galactic   (1950) -> Equatorial (1950)" << endl
         << "          9: Equatorial (1950) -> Ecliptic   (1950)" << endl
         << "         10: Ecliptic   (1950) -> Equatorial (1950)" << endl
         << "         11: Ecliptic   (1950) -> Galactic   (1950)" << endl
         << "         12: Galactic   (1950) -> Ecliptic   (1950)" << endl
         << endl
         << "pol: T or F" << endl << endl;
    throw Message_error();
    }

  string infile  = argv[1];
  string outfile = argv[2];
  int trafo = stringToData<int>(argv[3]);
  bool polarisation = stringToData<bool>(argv[4]);

  Trafo tr(maketrafo(trafo));

  fitshandle out;
  out.create (outfile);

  Alm<xcomplex<double> > almT,almG,almC;

  if (!polarisation)
    {
    int lmax, dummy;
    get_almsize (infile,lmax,dummy);
    read_Alm_from_fits (infile, almT, lmax, lmax);
    rotate_alm(almT,tr.Matrix());
    write_Alm_to_fits (out,almT,lmax,lmax,TFLOAT);
    }
  else
    {
    int lmax, dummy;
    get_almsize_pol (infile,lmax,dummy);
    read_Alm_from_fits (infile, almT, lmax, lmax, 2);
    read_Alm_from_fits (infile, almG, lmax, lmax, 3);
    read_Alm_from_fits (infile, almC, lmax, lmax, 4);
    rotate_alm(almT,almG,almC,tr.Matrix());
    write_Alm_to_fits (out,almT,lmax,lmax,TDOUBLE);
    write_Alm_to_fits (out,almG,lmax,lmax,TDOUBLE);
    write_Alm_to_fits (out,almC,lmax,lmax,TDOUBLE);
    }

PLANCK_DIAGNOSIS_END
  }

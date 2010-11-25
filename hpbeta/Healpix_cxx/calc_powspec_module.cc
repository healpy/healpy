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
 *  Copyright (C) 2004-2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "fitshandle.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_powspec_tools.h"
#include "levels_facilities.h"

using namespace std;

int calc_powspec_module (int argc, const char **argv)
  {
  announce ("calc_powspec");
  planck_assert (argc==3||argc==4,
    "usage: calc_powspec <almfile1> [<almfile2>] <powspec_file>");

  if (argc==3)
    {
    int lmax,mmax;
    get_almsize(argv[1],lmax,mmax,2);
    Alm<xcomplex<float> > alm;
    read_Alm_from_fits (argv[1],alm,lmax,mmax,2);
    PowSpec powspec;
    extract_powspec (alm,powspec);
    write_powspec_to_fits (argv[2],powspec,1);
    }
  else
    {
    int lmax,mmax;
    get_almsize(argv[1],lmax,mmax,2);
    Alm<xcomplex<float> > alm1;
    read_Alm_from_fits (argv[1],alm1,lmax,mmax,2);
    get_almsize(argv[2],lmax,mmax,2);
    Alm<xcomplex<float> > alm2;
    read_Alm_from_fits (argv[2],alm2,lmax,mmax,2);
    PowSpec powspec;
    extract_crosspowspec (alm1,alm2,powspec);
    write_powspec_to_fits (argv[3],powspec,1);
    }

  return 0;
  }

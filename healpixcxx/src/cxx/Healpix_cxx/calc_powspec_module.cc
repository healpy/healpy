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
 *  Copyright (C) 2004-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "fitshandle.h"
#include "paramfile.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_powspec_tools.h"
#include "levels_facilities.h"
#include "announce.h"

using namespace std;

int calc_powspec_module (int argc, const char **argv)
  {
  module_startup ("calc_powspec",argc,argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  bool pol=params.find<bool>("pol",false);
  string alm1=params.find<string>("alm1");
  string ps=params.find<string>("ps");

  if (!params.param_present("alm2"))
    {
    int lmax,mmax;
    pol ? get_almsize_pol(alm1,lmax,mmax)
        : get_almsize    (alm1,lmax,mmax);
    Alm<xcomplex<float> > almT, almG, almC;
    pol ? read_Alm_from_fits (alm1,almT,almG,almC,lmax,mmax)
        : read_Alm_from_fits (alm1,almT,lmax,mmax);
    PowSpec powspec;
    pol ? extract_powspec (almT,almG,almC,powspec)
        : extract_powspec (almT,powspec);
    write_powspec_to_fits (ps,powspec,pol ? 6 : 1);
    }
  else
    {
    planck_assert(!pol, "polarisation not supported for cross-powerspectra");
    int lmax,mmax;
    get_almsize(alm1,lmax,mmax);
    Alm<xcomplex<float> > Alm1;
    read_Alm_from_fits (alm1,Alm1,lmax,mmax);
    string alm2=params.find<string>("alm2");
    get_almsize(alm2,lmax,mmax);
    Alm<xcomplex<float> > Alm2;
    read_Alm_from_fits (alm2,Alm2,lmax,mmax);
    PowSpec powspec;
    extract_crosspowspec (Alm1,Alm2,powspec);
    write_powspec_to_fits (ps,powspec,1);
    }

  return 0;
  }

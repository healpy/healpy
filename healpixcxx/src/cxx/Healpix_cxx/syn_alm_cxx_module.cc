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
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "paramfile.h"
#include "planck_rng.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_powspec_tools.h"
#include "fitshandle.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

namespace {

template<typename T> void syn_alm_cxx (paramfile &params)
  {
  int nlmax = params.template find<int>("nlmax");
  int nmmax = params.template find<int>("nmmax",nlmax);
  planck_assert(nmmax<=nlmax,"nmmax must not be larger than nlmax");
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int rand_seed = params.template find<int>("rand_seed");
  double fwhm = arcmin2rad*params.template find<double>("fwhm_arcmin",0.);
  bool polarisation = params.template find<bool>("polarisation");

  PowSpec powspec;
  int nspecs = polarisation ? 4 : 1;
  read_powspec_from_fits (infile, powspec, nspecs, nlmax);
  powspec.smoothWithGauss(fwhm);

  planck_rng rng(rand_seed);

  if (polarisation)
    {
    Alm<xcomplex<T> >
      almT(nlmax,nmmax), almG(nlmax,nmmax), almC(nlmax,nmmax);
    create_alm_pol (powspec, almT, almG, almC, rng);
    write_Alm_to_fits(outfile,almT,almG,almC,nlmax,nmmax,planckType<T>());
    }
  else
    {
    Alm<xcomplex<T> > almT(nlmax,nmmax);
    create_alm (powspec, almT, rng);
    write_Alm_to_fits(outfile,almT,nlmax,nmmax,planckType<T>());
    }
  }

} // unnamed namespace

int syn_alm_cxx_module (int argc, const char **argv)
  {
  module_startup ("syn_alm_cxx", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  bool dp = params.find<bool> ("double_precision",false);
  dp ? syn_alm_cxx<double>(params) : syn_alm_cxx<float>(params);
  return 0;
  }

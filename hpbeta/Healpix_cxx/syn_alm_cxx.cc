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
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "simparams.h"
#include "planck_rng.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "fitshandle.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_powspec_tools.h"

using namespace std;

template<typename T> void syn_alm_cxx (paramfile &params, simparams &par)
  {
  par.add_comment("-----------------------");
  par.add_comment(" ** syn_alm_cxx 1.0 **");
  par.add_comment("-----------------------");

  int nlmax = params.template find<int>("nlmax");
  par.add("nlmax","NLMAX",nlmax,"maximum l of the alms");
  int nmmax = params.template find<int>("nmmax",nlmax);
  par.add("nmmax","NMMAX",nmmax,"maximum m of the alms");
  planck_assert(nmmax<=nlmax,"nmmax must not be larger than nlmax");
  string infile = params.template find<string>("infile");
  par.add_source_file (infile, 2);
  string outfile = params.template find<string>("outfile");
  int rand_seed = params.template find<int>("rand_seed");
  par.add("rand_seed","RANDSEED",rand_seed,"random number seed");
  double fwhm = params.template find<double>("fwhm_arcmin",0.);
  par.add("fwhm","FWHM",fwhm,"[arcmin] Gaussian smoothing FWHM");
  fwhm *= degr2rad/60;
  bool polarisation = params.template find<bool>("polarisation");

  PowSpec powspec;
  int nspecs = polarisation ? 4 : 1;
  read_powspec_from_fits (infile, powspec, nspecs, nlmax);
  powspec.Smooth_with_Gauss(fwhm);

  planck_rng rng(rand_seed);

  if (polarisation)
    {
    Alm<xcomplex<T> >
      almT(nlmax,nmmax), almG(nlmax,nmmax), almC(nlmax,nmmax);
    create_alm_pol (powspec, almT, almG, almC, rng);

    fitshandle out;
    out.create(outfile);
    write_Alm_to_fits(out,almT,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    out.add_key("PDMTYPE",string("MAPALM"),"Planck data model type");
    par.add_keys(out);
    write_Alm_to_fits(out,almG,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    out.add_key("PDMTYPE",string("MAPALM"),"Planck data model type");
    par.add_keys(out);
    write_Alm_to_fits(out,almC,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    out.add_key("PDMTYPE",string("MAPALM"),"Planck data model type");
    par.add_keys(out);
    }
  else
    {
    Alm<xcomplex<T> > almT(nlmax,nmmax);
    create_alm (powspec, almT, rng);

    fitshandle out;
    out.create(outfile);
    write_Alm_to_fits(out,almT,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    out.add_key("PDMTYPE",string("MAPALM"),"Planck data model type");
    par.add_keys(out);
    }
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("syn_alm_cxx", argc, argv, 2, "<parameter file>");
  paramfile params (argv[1]);
  simparams par;

  bool dp = params.find<bool> ("double_precision",false);
  dp ? syn_alm_cxx<double>(params,par) : syn_alm_cxx<float>(params,par);
PLANCK_DIAGNOSIS_END
  }

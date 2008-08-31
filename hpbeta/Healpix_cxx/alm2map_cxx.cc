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
#include "healpix_data_io.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "fitshandle.h"

using namespace std;

template<typename T> void alm2map_cxx (paramfile &params, simparams &par)
  {
  par.add_comment("-----------------------");
  par.add_comment(" ** alm2map_cxx 1.0 **");
  par.add_comment("-----------------------");

  int nlmax = params.template find<int>("nlmax");
  par.add("nlmax","NLMAX",nlmax,"maximum l of the alms");
  int nmmax = params.template find<int>("nmmax",nlmax);
  par.add("nmmax","NMMAX",nmmax,"maximum m of the alms");
  planck_assert(nmmax<=nlmax,"nmmax must not be larger than nlmax");
  string infile = params.template find<string>("infile");
  par.add_source_file (infile, 2);
  string outfile = params.template find<string>("outfile");
  int nside = params.template find<int>("nside");
  double fwhm_arcmin = params.template find<double>("fwhm_arcmin",0);

  arr<double> temp, pol;
  get_pixwin (params,par,nlmax,nside,temp,pol);

  bool deriv = params.template find<bool>("derivatives",false);
  if (deriv)
    {
    Alm<xcomplex<T> > alm;
    read_Alm_from_fits(infile,alm,nlmax,nmmax,2);
    if (fwhm_arcmin>0) smooth_with_Gauss (alm, fwhm_arcmin);
    Healpix_Map<T> map(nside,RING,SET_NSIDE),
                   mapdth(nside,RING,SET_NSIDE),
                   mapdph(nside,RING,SET_NSIDE);
    alm.ScaleL(temp);

    double offset = alm(0,0).real()/sqrt(fourpi);
    alm(0,0) = 0;
    alm2map_der1(alm,map,mapdth,mapdph);
    for (int m=0; m<map.Npix(); ++m) map[m]+=offset;
    fitshandle out;
    out.create (outfile);
    write_Healpix_map_to_fits (out,map,mapdth,mapdph,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    return;
    }

  bool polarisation = params.template find<bool>("polarisation");
  if (!polarisation)
    {
    Alm<xcomplex<T> > alm;
    read_Alm_from_fits(infile,alm,nlmax,nmmax,2);
    if (fwhm_arcmin>0) smooth_with_Gauss (alm, fwhm_arcmin);
    Healpix_Map<T> map(nside,RING,SET_NSIDE);
    alm.ScaleL(temp);

    double offset = alm(0,0).real()/sqrt(fourpi);
    alm(0,0) = 0;
    alm2map(alm,map);
    map.add(offset);
    fitshandle out;
    out.create (outfile);
    write_Healpix_map_to_fits (out,map,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    }
  else
    {
    Alm<xcomplex<T> > almT, almG, almC;
    read_Alm_from_fits(infile,almT,nlmax,nmmax,2);
    read_Alm_from_fits(infile,almG,nlmax,nmmax,3);
    read_Alm_from_fits(infile,almC,nlmax,nmmax,4);
    if (fwhm_arcmin>0) smooth_with_Gauss (almT, almG, almC, fwhm_arcmin);
    Healpix_Map<T> mapT(nside,RING,SET_NSIDE), mapQ(nside,RING,SET_NSIDE),
                   mapU(nside,RING,SET_NSIDE);
    almT.ScaleL(temp);
    almG.ScaleL(pol); almC.ScaleL(pol);

    double offset = almT(0,0).real()/sqrt(fourpi);
    almT(0,0) = 0;
    alm2map_pol(almT,almG,almC,mapT,mapQ,mapU);
    mapT.add(offset);
    fitshandle out;
    out.create (outfile);
    write_Healpix_map_to_fits (out,mapT,mapQ,mapU,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    }
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("alm2map_cxx", argc, argv, 2, "<parameter file>");
  paramfile params (argv[1]);
  simparams par;

  bool dp = params.find<bool> ("double_precision",false);
  dp ? alm2map_cxx<double>(params,par) : alm2map_cxx<float>(params,par);
PLANCK_DIAGNOSIS_END
  }

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
#include "healpix_data_io.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "alm_powspec_tools.h"
#include "fitshandle.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"

using namespace std;

namespace {

template<typename T> void mult_alm (paramfile &params)
  {
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int nside_pixwin_in = params.template find<int>("nside_pixwin_in",0);
  planck_assert (nside_pixwin_in>=0,"nside_pixwin_in must be >= 0");
  int nside_pixwin_out = params.template find<int>("nside_pixwin_out",0);
  planck_assert (nside_pixwin_out>=0,"nside_pixwin_out must be >= 0");
  string cl_in = params.template find<string>("cl_in","");
  string cl_out = params.template find<string>("cl_out","");
  double fwhm_in = arcmin2rad*params.template find<double>("fwhm_arcmin_in",0);
  planck_assert (fwhm_in>=0,"fwhm_arcmin_in must be >= 0");
  double fwhm_out = arcmin2rad*params.template find<double>("fwhm_arcmin_out",0);
  planck_assert (fwhm_out>=0,"fwhm_arcmin_out must be >= 0");

  string datadir;
  if ((nside_pixwin_in>0) || (nside_pixwin_out>0))
    datadir = params.template find<string>("healpix_data");

  bool polarisation = params.template find<bool>("polarisation");
  if (!polarisation)
    {
    int nlmax, nmmax;
    get_almsize(infile, nlmax, nmmax, 2);
    Alm<xcomplex<T> > alm;
    read_Alm_from_fits(infile,alm,nlmax,nmmax,2);
    if (fwhm_in>0) smoothWithGauss (alm, -fwhm_in);
    arr<double> temp(nlmax+1);
    PowSpec tps;
    if (nside_pixwin_in>0)
      {
      read_pixwin(datadir,nside_pixwin_in,temp);
      for (int l=0; l<=nlmax; ++l)
        temp[l] = 1/temp[l];
      alm.ScaleL (temp);
      }
    if (cl_in!="")
      {
      read_powspec_from_fits (cl_in,tps,1,alm.Lmax());
      for (int l=0; l<=nlmax; ++l)
        temp[l] = 1./sqrt(tps.tt(l));
      alm.ScaleL (temp);
      }
    if (nside_pixwin_out>0)
      {
      read_pixwin(datadir,nside_pixwin_out,temp);
      alm.ScaleL (temp);
      }
    if (cl_out!="")
      {
      read_powspec_from_fits (cl_out,tps,1,alm.Lmax());
      for (int l=0; l<=nlmax; ++l)
        temp[l] = sqrt(tps.tt(l));
      alm.ScaleL (temp);
      }
    if (fwhm_out>0) smoothWithGauss (alm, fwhm_out);
    write_Alm_to_fits (outfile,alm,nlmax,nmmax,planckType<T>());
    }
  else
    {
    int nlmax, nmmax;
    get_almsize_pol(infile, nlmax, nmmax);
    Alm<xcomplex<T> > almT, almG, almC;
    read_Alm_from_fits(infile,almT,almG,almC,nlmax,nmmax,2);
    if (fwhm_in>0) smoothWithGauss (almT, almG, almC, -fwhm_in);
    arr<double> temp(nlmax+1), pol(nlmax+1);
    if (nside_pixwin_in>0)
      {
      read_pixwin(datadir,nside_pixwin_in,temp,pol);
      for (int l=0; l<=nlmax; ++l)
        { temp[l] = 1/temp[l]; if (pol[l]!=0) pol[l] = 1/pol[l]; }
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }
    if (cl_in!="")
      planck_fail ("power spectra not (yet) supported with polarisation");
    if (nside_pixwin_out>0)
      {
      read_pixwin(datadir,nside_pixwin_out,temp,pol);
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }
    if (cl_out!="")
      planck_fail ("power spectra not (yet) supported with polarisation");
    if (fwhm_out>0) smoothWithGauss (almT, almG, almC, fwhm_out);
    write_Alm_to_fits (outfile,almT,almG,almC,nlmax,nmmax,planckType<T>());
    }
  }

} // unnamed namespace

int mult_alm_module (int argc, const char **argv)
  {
  module_startup ("mult_alm", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  bool dp = params.find<bool> ("double_precision",false);
  dp ? mult_alm<double>(params) : mult_alm<float>(params);

  return 0;
  }

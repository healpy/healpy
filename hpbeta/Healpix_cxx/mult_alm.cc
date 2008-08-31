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
#include "alm_powspec_tools.h"
#include "fitshandle.h"

using namespace std;

template<typename T> void mult_alm (paramfile &params, simparams &par)
  {
  par.add_comment("----------------");
  par.add_comment(" ** mult_alm **");
  par.add_comment("----------------");

  string infile = params.template find<string>("infile");
  par.add_source_file (infile, 2);
  string outfile = params.template find<string>("outfile");
  int nside_pixwin_in = params.template find<int>("nside_pixwin_in",0);
  planck_assert (nside_pixwin_in>=0,"nside_pixwin_in must be >= 0");
  int nside_pixwin_out = params.template find<int>("nside_pixwin_out",0);
  planck_assert (nside_pixwin_out>=0,"nside_pixwin_out must be >= 0");
  double fwhm_arcmin_in = params.template find<double>("fwhm_arcmin_in",0);
  planck_assert (fwhm_arcmin_in>=0,"fwhm_arcmin_in must be >= 0");
  double fwhm_arcmin_out = params.template find<double>("fwhm_arcmin_out",0);
  planck_assert (fwhm_arcmin_out>=0,"fwhm_arcmin_out must be >= 0");

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
    if (fwhm_arcmin_in>0) smooth_with_Gauss (alm, -fwhm_arcmin_in);
    arr<double> temp(nlmax+1);
    if (nside_pixwin_in>0)
      {
      read_pixwin(datadir,nside_pixwin_in,temp);
      for (int l=0; l<=nlmax; ++l)
        temp[l] = 1/temp[l];
      alm.ScaleL (temp);
      }
    if (nside_pixwin_out>0)
      {
      read_pixwin(datadir,nside_pixwin_out,temp);
      alm.ScaleL (temp);
      }
    if (fwhm_arcmin_out>0) smooth_with_Gauss (alm, fwhm_arcmin_out);
    fitshandle out;
    out.create (outfile);
    write_Alm_to_fits (out,alm,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    }
  else
    {
    int nlmax, nmmax;
    get_almsize_pol(infile, nlmax, nmmax);
    Alm<xcomplex<T> > almT, almG, almC;
    read_Alm_from_fits(infile,almT,nlmax,nmmax,2);
    read_Alm_from_fits(infile,almG,nlmax,nmmax,3);
    read_Alm_from_fits(infile,almC,nlmax,nmmax,4);
    if (fwhm_arcmin_in>0)
      smooth_with_Gauss (almT, almG, almC, -fwhm_arcmin_in);
    arr<double> temp(nlmax+1), pol(nlmax+1);
    if (nside_pixwin_in>0)
      {
      read_pixwin(datadir,nside_pixwin_in,temp,pol);
      for (int l=0; l<=nlmax; ++l)
        { temp[l] = 1/temp[l]; pol[l] = 1/pol[l]; }
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }
    if (nside_pixwin_out>0)
      {
      read_pixwin(datadir,nside_pixwin_out,temp,pol);
      almT.ScaleL(temp); almG.ScaleL(pol); almC.ScaleL(pol);
      }
    if (fwhm_arcmin_out>0)
      smooth_with_Gauss (almT, almG, almC, fwhm_arcmin_out);
    fitshandle out;
    out.create (outfile);
    write_Alm_to_fits (out,almT,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    write_Alm_to_fits (out,almG,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    write_Alm_to_fits (out,almC,nlmax,nmmax,FITSUTIL<T>::DTYPE);
    par.add_keys(out);
    }
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("mult_alm", argc, argv, 2, "<parameter file>");
  paramfile params (argv[1]);
  simparams par;

  bool dp = params.find<bool> ("double_precision",false);
  dp ? mult_alm<double>(params,par) : mult_alm<float>(params,par);
PLANCK_DIAGNOSIS_END
  }

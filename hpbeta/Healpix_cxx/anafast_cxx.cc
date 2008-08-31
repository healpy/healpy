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
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "fitshandle.h"

using namespace std;

template<typename T> void anafast_cxx (paramfile &params, simparams &par)
  {
  int nlmax = params.template find<int>("nlmax");
  int nmmax = params.template find<int>("nmmax",nlmax);
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile","");
  string outfile_alms = params.template find<string>("outfile_alms","");
  planck_assert ((outfile!="") || (outfile_alms!=""),
    "no output specified, nothing done");
  bool polarisation = params.template find<bool>("polarisation");
  int num_iter = params.template find<int>("iter_order",0);

  if (!polarisation)
    {
    Healpix_Map<T> map;
    read_Healpix_map_from_fits(infile,map,1,2);
    arr<double> weight;
    get_ring_weights (params,par,map.Nside(),weight);

    Alm<xcomplex<T> > alm(nlmax,nmmax);
    double avg=map.average();
    map.add(-avg);
    if (map.Scheme()==NEST) map.swap_scheme();
    map2alm_iter(map,alm,num_iter,weight);

    alm(0,0) += avg*sqrt(fourpi);

    if (outfile!="")
      {
      PowSpec powspec;
      extract_powspec (alm,powspec);
      fitshandle out;
      out.create (outfile);
      write_powspec_to_fits (out,powspec,1);
      }
    if (outfile_alms!="")
      {
      fitshandle out;
      out.create (outfile_alms);
      write_Alm_to_fits (out,alm,alm.Lmax(),alm.Mmax(),FITSUTIL<T>::DTYPE);
      }
    }
  else
    {
    Healpix_Map<T> mapT, mapQ, mapU;
    read_Healpix_map_from_fits(infile,mapT,1,2);
    read_Healpix_map_from_fits(infile,mapQ,2,2);
    read_Healpix_map_from_fits(infile,mapU,3,2);
    arr<double> weight;
    get_ring_weights (params,par,mapT.Nside(),weight);

    Alm<xcomplex<T> > almT(nlmax,nmmax), almG(nlmax,nmmax), almC(nlmax,nmmax);
    double avg=mapT.average();
    mapT.add(-avg);
    if (mapT.Scheme()==NEST) mapT.swap_scheme();
    if (mapQ.Scheme()==NEST) mapQ.swap_scheme();
    if (mapU.Scheme()==NEST) mapU.swap_scheme();
    map2alm_pol_iter
      (mapT,mapQ,mapU,almT,almG,almC,num_iter,weight);

    almT(0,0) += avg*sqrt(fourpi);

    if (outfile!="")
      {
      PowSpec powspec;
      extract_powspec (almT,almG,almC,powspec);
      fitshandle out;
      out.create (outfile);
      write_powspec_to_fits (out,powspec,4);
      }
    if (outfile_alms!="")
      {
      fitshandle out;
      out.create (outfile_alms);
      write_Alm_to_fits (out,almT,almT.Lmax(),almT.Mmax(),FITSUTIL<T>::DTYPE);
      write_Alm_to_fits (out,almG,almG.Lmax(),almG.Mmax(),FITSUTIL<T>::DTYPE);
      write_Alm_to_fits (out,almC,almC.Lmax(),almC.Mmax(),FITSUTIL<T>::DTYPE);
      }
    }
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("anafast_cxx", argc, argv, 2, "<parameter file>");
  paramfile params (argv[1]);
  simparams par;

  bool dp = params.find<bool> ("double_precision",false);
  dp ? anafast_cxx<double>(params,par) : anafast_cxx<float>(params,par);
PLANCK_DIAGNOSIS_END
  }

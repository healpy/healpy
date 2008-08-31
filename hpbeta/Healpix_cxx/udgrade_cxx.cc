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
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "cxxutils.h"
#include "paramfile.h"
#include "simparams.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"

using namespace std;

template<typename T> void udgrade_cxx (paramfile &params, simparams &par)
  {
  par.add_comment("-------------------");
  par.add_comment(" ** udgrade_cxx **");
  par.add_comment("-------------------");

  string infile = params.template find<string>("infile");
  par.add_source_file (infile, 2);
  string outfile = params.template find<string>("outfile");
  int order = Healpix_Base::nside2order (params.template find<int>("nside"));
  bool polarisation = params.template find<bool>("polarisation",false);
  bool pessimistic = params.template find<bool>("pessimistic",false);

  if (!polarisation)
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_fits(infile,inmap,1,2);
    Healpix_Map<T> outmap (order, inmap.Scheme());

    outmap.Import(inmap,pessimistic);
    fitshandle out;
    out.create(outfile);
    write_Healpix_map_to_fits (out,outmap,FITSUTIL<T>::DTYPE);
    }
  else
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_fits(infile,inmap,1,2);
    Healpix_Map<T> outmapT (order, inmap.Scheme()),
                   outmapQ (order, inmap.Scheme()),
                   outmapU (order, inmap.Scheme());

    outmapT.Import(inmap,pessimistic);
    read_Healpix_map_from_fits(infile,inmap,2,2);
    outmapQ.Import(inmap,pessimistic);
    read_Healpix_map_from_fits(infile,inmap,3,2);
    outmapU.Import(inmap,pessimistic);
    fitshandle out;
    out.create(outfile);
    write_Healpix_map_to_fits (out,outmapT,outmapQ,outmapU,FITSUTIL<T>::DTYPE);
    }
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  module_startup ("udgrade_cxx", argc, argv, 2, "<parameter file>");
  paramfile params (argv[1]);
  simparams par;

  bool dp = params.find<bool> ("double_precision",false);
  dp ? udgrade_cxx<double>(params,par) : udgrade_cxx<float>(params,par);
PLANCK_DIAGNOSIS_END
  }

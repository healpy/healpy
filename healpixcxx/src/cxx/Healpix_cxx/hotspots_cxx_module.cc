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
 *  Copyright (C) 2003-2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <fstream>
#include "paramfile.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"
#include "levels_facilities.h"
#include "announce.h"

using namespace std;

int hotspots_cxx_module (int argc, const char **argv)
  {
  module_startup ("hotspots_cxx", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  string infile = params.find<string>("infile");
  string mapfile = params.find<string>("outmap","");
  bool have_mapfile = mapfile!="";
  string minfile = params.find<string>("minfile","");
  bool have_minfile = minfile!="";
  string maxfile = params.find<string>("maxfile","");
  bool have_maxfile = maxfile!="";
  planck_assert (have_mapfile || have_minfile || have_maxfile,
    "no output file specified");

  Healpix_Map<float> inmap;
  read_Healpix_map_from_fits(infile,inmap,1,2);
  Healpix_Map<float> outmap;
  if (have_mapfile) outmap.Set(inmap.Order(),inmap.Scheme());

  ofstream minout, maxout;
  if (have_minfile) minout.open(minfile.c_str());
  if (have_maxfile) maxout.open(maxfile.c_str());

  fix_arr<int,8> nb;
// FIXME: This should be parallelized
  for (int m=0; m<inmap.Npix(); ++m)
    {
    float value = inmap[m];
    if (!approx<double>(value, Healpix_undef))
      {
      inmap.neighbors(m,nb);
      bool ismax=true, ismin=true;
      for (tsize n=0; n<nb.size(); ++n)
        {
        if (nb[n] >=0)
          {
          float nbval = inmap[nb[n]];
          if (!approx<double>(nbval, Healpix_undef))
            {
            if (nbval>=value) ismax=false;
            if (nbval<=value) ismin=false;
            }
          }
        }
      if (have_mapfile)
        outmap[m] = float((ismax||ismin) ? value : Healpix_undef);
      if (have_minfile && ismin) minout << m << " " << value << endl;
      if (have_maxfile && ismax) maxout << m << " " << value << endl;
      }
    }
  if (have_mapfile)
    write_Healpix_map_to_fits (mapfile,outmap,PLANCK_FLOAT32);

  return 0;
  }

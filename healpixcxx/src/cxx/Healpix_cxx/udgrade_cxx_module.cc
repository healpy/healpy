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
 *  Copyright (C) 2003-2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "paramfile.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"
#include "levels_facilities.h"
#include "geom_utils.h"
#include "xcomplex.h"
#include "announce.h"

using namespace std;

namespace {

double alpha (const pointing &p1, const pointing &p2)
  {
  vec3 v1(p1), v2(p2);
  vec3 dir(v2-v1);
  return orientation (v2, dir) - orientation (v1,dir);
  }

template<typename T> void degrade_pol (const Healpix_Map<T> &q1,
  const Healpix_Map<T> &u1, Healpix_Map<T> &q2, Healpix_Map<T> &u2,
  bool pessimistic)
  {
  planck_assert(q1.conformable(u1) && q2.conformable(u2), "map mismatch");
  planck_assert(q2.Nside()<q1.Nside(),"this is no degrade");
  int fact = q1.Nside()/q2.Nside();
  planck_assert (q1.Nside()==q2.Nside()*fact,
    "the larger Nside must be a multiple of the smaller one");

  int minhits = pessimistic ? fact*fact : 1;
#pragma omp parallel
{
  int m, npix=q2.Npix();
#pragma omp for schedule (static)
  for (m=0; m<npix; ++m)
    {
    int x,y,f;
    q2.pix2xyf(m,x,y,f);
    int hits = 0;
    xcomplex<double> sum = 0;
    for (int j=fact*y; j<fact*(y+1); ++j)
      for (int i=fact*x; i<fact*(x+1); ++i)
        {
        int opix = q1.xyf2pix(i,j,f);
        if (!(approx<double>(q1[opix],Healpix_undef)
            ||approx<double>(u1[opix],Healpix_undef)))
          {
          ++hits;
          xcomplex<double> val(q1[opix],u1[opix]);
          double ang=alpha(q2.pix2ang(m),q1.pix2ang(opix));
          xcomplex<double> mul(cos(2*ang),sin(2*ang));
          sum += val*mul;
          }
        }
    q2[m] = T((hits<minhits) ? Healpix_undef : sum.re/hits);
    u2[m] = T((hits<minhits) ? Healpix_undef : sum.im/hits);
    }
}
  }

template<typename T> void udgrade_cxx (paramfile &params)
  {
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int nside = params.template find<int>("nside");
  bool polarisation = params.template find<bool>("polarisation",false);
  bool pessimistic = params.template find<bool>("pessimistic",false);

  if (!polarisation)
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_fits(infile,inmap,1,2);
    Healpix_Map<T> outmap (nside,inmap.Scheme(),SET_NSIDE);

    outmap.Import(inmap,pessimistic);
    write_Healpix_map_to_fits (outfile,outmap,planckType<T>());
    }
  else
    {
    Healpix_Map<T> inmap;
    read_Healpix_map_from_fits(infile,inmap,1,2);
    Healpix_Map<T> outmapT (nside,inmap.Scheme(),SET_NSIDE),
                   outmapQ (nside,inmap.Scheme(),SET_NSIDE),
                   outmapU (nside,inmap.Scheme(),SET_NSIDE);
//    planck_assert(inmap.Nside()<=outmapT.Nside(),
//      "degrading not supported for polarised maps");

    outmapT.Import(inmap,pessimistic);
    if ((outmapQ.Nside()<inmap.Nside())
      && params.template find<bool>("parallel_transport",true))
      {
cout << "Experimental: polarised degrade with parallel transport" << endl;
      read_Healpix_map_from_fits(infile,inmap,2,2);
      Healpix_Map<T> inmap2;
      read_Healpix_map_from_fits(infile,inmap2,3,2);
      degrade_pol (inmap, inmap2, outmapQ, outmapU, pessimistic);
      }
    else
      {
cout << "WARNING: polarised degrade without parallel transport" << endl;
      read_Healpix_map_from_fits(infile,inmap,2,2);
      outmapQ.Import(inmap,pessimistic);
      read_Healpix_map_from_fits(infile,inmap,3,2);
      outmapU.Import(inmap,pessimistic);
      }
    write_Healpix_map_to_fits (outfile,outmapT,outmapQ,outmapU,planckType<T>());
    }
  }

} // unnamed namespace

int udgrade_cxx_module (int argc, const char **argv)
  {
  module_startup ("udgrade_cxx", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  bool dp = params.find<bool> ("double_precision",false);
  dp ? udgrade_cxx<double>(params) : udgrade_cxx<float>(params);

  return 0;
  }

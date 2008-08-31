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
 *  Copyright (C) 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "cxxutils.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"

using namespace std;

template<typename Iterator> typename iterator_traits<Iterator>::value_type
  median(Iterator first, Iterator last)
  {
  Iterator mid = first+(last-first-1)/2;
  nth_element(first,mid,last);
  if ((last-first)&1) return *mid;
  return 0.5*((*mid)+(*min_element(mid+1,last)));
  }

int main (int argc, const char **argv)
  {
  module_startup ("median_filter", argc, argv, 4,
                  "<input map> <output map> <radius in arcmin>");

  double radius = stringToData<double>(argv[3])/60*degr2rad;

  Healpix_Map<float> inmap;
  read_Healpix_map_from_fits(argv[1],inmap,1,2);
  Healpix_Map<float> outmap (inmap.Nside(), inmap.Scheme(), SET_NSIDE);

  vector<int> list;
  vector<float> list2;
  for (int m=0; m<inmap.Npix(); ++m)
    {
    inmap.query_disc(inmap.pix2ang(m),radius,list);
    list2.resize(list.size());
    for (unsigned int i=0; i<list.size(); ++i)
      list2[i] = inmap[list[i]];
    outmap[m]=median(list2.begin(),list2.end());
    }

  fitshandle out;
  out.create(argv[2]);
  write_Healpix_map_to_fits (out,outmap,TFLOAT);
  }

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

#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"

using namespace std;

template<typename T> void read_Healpix_map_from_fits
  (fitshandle &inp, Healpix_Map<T> &map, int colnum)
  {
  string ordering;
  inp.get_key ("ORDERING", ordering);
  arr<T> myarr;
  inp.read_entire_column (colnum, myarr);
  map.Set (myarr, ordering=="RING" ? RING : NEST);
  }

template void read_Healpix_map_from_fits (fitshandle &inp,
  Healpix_Map<float> &map, int colnum);
template void read_Healpix_map_from_fits (fitshandle &inp,
  Healpix_Map<double> &map, int colnum);
template void read_Healpix_map_from_fits (fitshandle &inp,
  Healpix_Map<int> &map, int colnum);


template<typename T> void read_Healpix_map_from_fits
  (const string &filename, Healpix_Map<T> &map, int colnum, int hdunum)
  {
  fitshandle inp;
  inp.open (filename);
  inp.goto_hdu (hdunum);
  read_Healpix_map_from_fits (inp,map,colnum);
  }

template void read_Healpix_map_from_fits (const string &filename,
  Healpix_Map<float> &map, int colnum, int hdunum);
template void read_Healpix_map_from_fits (const string &filename,
  Healpix_Map<double> &map, int colnum, int hdunum);
template void read_Healpix_map_from_fits (const string &filename,
  Healpix_Map<int> &map, int colnum, int hdunum);


void prepare_Healpix_fitsmap
  (fitshandle &out, const Healpix_Base &base, int datatype,
  const arr<string> &colname)
  {
  vector<fitscolumn> cols;
  int repcount = healpix_repcount (base.Npix());
  for (int m=0; m<colname.size(); ++m)
    cols.push_back (fitscolumn (colname[m],"unknown",repcount, datatype));
  out.insert_bintab(cols);
  out.add_key ("PIXTYPE",string("HEALPIX"),"HEALPIX pixelisation");
  string ordering = (base.Scheme()==RING) ? "RING" : "NESTED";
  out.add_key ("ORDERING",ordering,
               "Pixel ordering scheme, either RING or NESTED");
  out.add_key ("NSIDE",base.Nside(),"Resolution parameter for HEALPIX");
  out.add_key ("FIRSTPIX",0,"First pixel # (0 based)");
  out.add_key ("LASTPIX",base.Npix()-1,"Last pixel # (0 based)");
  out.add_key ("INDXSCHM",string("IMPLICIT"),"Indexing: IMPLICIT or EXPLICIT");
  }

template<typename T> void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<T> &map, int datatype)
  {
  arr<string> colname(1);
  colname[0] = "signal";
  prepare_Healpix_fitsmap (out, map, datatype, colname);
  out.write_column(1,map.Map());
  }

template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<float> &map, int datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<double> &map, int datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<int> &map, int datatype);


template<typename T> void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, int datatype)
  {
  arr<string> colname(3);
  colname[0] = "signal";
  colname[1] = "Q-pol";
  colname[2] = "U-pol";
  prepare_Healpix_fitsmap (out, mapT, datatype, colname);
  out.write_column(1,mapT.Map());
  out.write_column(2,mapQ.Map());
  out.write_column(3,mapU.Map());
  }

template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<float> &mapT,
   const Healpix_Map<float> &mapQ, const Healpix_Map<float> &mapU,
   int datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ, const Healpix_Map<double> &mapU,
   int datatype);

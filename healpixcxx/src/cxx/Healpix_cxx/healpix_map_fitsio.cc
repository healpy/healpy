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

#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"
#include "share_utils.h"
#include "string_utils.h"

using namespace std;

namespace {

unsigned int healpix_repcount (tsize npix)
  {
  if (npix<1024) return 1;
  if ((npix%1024)==0) return 1024;
  return isqrt (npix/12);
  }

} // unnamed namespace

template<typename T> void read_Healpix_map_from_fits
  (fitshandle &inp, Healpix_Map<T> &map, int colnum)
  {
  arr<T> myarr;
  inp.read_entire_column (colnum, myarr);
  int64 nside = inp.get_key<int>("NSIDE");
  planck_assert (int64(myarr.size())==12*nside*nside,
    string("mismatch between number of map pixels ("
    +dataToString(myarr.size())+") and Nside ("+dataToString(nside)+")"));
  map.Set (myarr, string2HealpixScheme(inp.get_key<string>("ORDERING")));
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

template<typename T> void read_Healpix_map_from_fits
  (fitshandle &inp, Healpix_Map<T> &mapT, Healpix_Map<T> &mapQ,
  Healpix_Map<T> &mapU)
  {
  int64 nside = inp.get_key<int>("NSIDE");
  Healpix_Ordering_Scheme scheme
    = string2HealpixScheme(inp.get_key<string>("ORDERING"));
  mapT.SetNside(nside,scheme);
  mapQ.SetNside(nside,scheme);
  mapU.SetNside(nside,scheme);
  planck_assert (multiequal(int64(mapT.Npix()),inp.nelems(1),inp.nelems(2),
    inp.nelems(3)), "mismatch between number of map pixels and Nside");
  chunkMaker cm(mapT.Npix(),inp.efficientChunkSize(1));
  uint64 offset,ppix;
  while(cm.getNext(offset,ppix))
    {
    inp.read_column_raw(1,&mapT[offset],ppix,offset);
    inp.read_column_raw(2,&mapQ[offset],ppix,offset);
    inp.read_column_raw(3,&mapU[offset],ppix,offset);
    }
  }

template void read_Healpix_map_from_fits (fitshandle &inp,
  Healpix_Map<float> &mapT, Healpix_Map<float> &mapQ,
  Healpix_Map<float> &mapU);
template void read_Healpix_map_from_fits (fitshandle &inp,
  Healpix_Map<double> &mapT, Healpix_Map<double> &mapQ,
  Healpix_Map<double> &mapU);

template<typename T> void read_Healpix_map_from_fits
  (const string &filename, Healpix_Map<T> &mapT, Healpix_Map<T> &mapQ,
  Healpix_Map<T> &mapU, int hdunum)
  {
  fitshandle inp;
  inp.open(filename);
  inp.goto_hdu(hdunum);
  read_Healpix_map_from_fits (inp,mapT,mapQ,mapU);
  }

template void read_Healpix_map_from_fits (const string &filename,
  Healpix_Map<float> &mapT, Healpix_Map<float> &mapQ,
  Healpix_Map<float> &mapU, int hdunum);
template void read_Healpix_map_from_fits (const string &filename,
  Healpix_Map<double> &mapT, Healpix_Map<double> &mapQ,
  Healpix_Map<double> &mapU, int hdunum);

void prepare_Healpix_fitsmap
  (fitshandle &out, const Healpix_Base &base, PDT datatype,
  const arr<string> &colname)
  {
  vector<fitscolumn> cols;
  int repcount = healpix_repcount (base.Npix());
  for (tsize m=0; m<colname.size(); ++m)
    cols.push_back (fitscolumn (colname[m],"unknown",repcount, datatype));
  out.insert_bintab(cols);
  out.set_key ("PIXTYPE",string("HEALPIX"),"HEALPIX pixelisation");
  string ordering = (base.Scheme()==RING) ? "RING" : "NESTED";
  out.set_key ("ORDERING",ordering,
                  "Pixel ordering scheme, either RING or NESTED");
  out.set_key ("NSIDE",base.Nside(),"Resolution parameter for HEALPIX");
  out.set_key ("FIRSTPIX",0,"First pixel # (0 based)");
  out.set_key ("LASTPIX",base.Npix()-1,"Last pixel # (0 based)");
  out.set_key ("INDXSCHM",string("IMPLICIT"),
                  "Indexing: IMPLICIT or EXPLICIT");
  }

template<typename T> void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<T> &map, PDT datatype)
  {
  arr<string> colname(1);
  colname[0] = "signal";
  prepare_Healpix_fitsmap (out, map, datatype, colname);
  out.write_column(1,map.Map());
  }

template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<float> &map, PDT datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<double> &map, PDT datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<int> &map, PDT datatype);


template<typename T> void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ, const Healpix_Map<T> &mapU, PDT datatype)
  {
  arr<string> colname(3);
  colname[0] = "signal";
  colname[1] = "Q-pol";
  colname[2] = "U-pol";
  prepare_Healpix_fitsmap (out, mapT, datatype, colname);
  chunkMaker cm(mapT.Npix(),out.efficientChunkSize(1));
  uint64 offset,ppix;
  while(cm.getNext(offset,ppix))
    {
    out.write_column_raw(1,&mapT[offset],ppix,offset);
    out.write_column_raw(2,&mapQ[offset],ppix,offset);
    out.write_column_raw(3,&mapU[offset],ppix,offset);
    }
  }

template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<float> &mapT,
   const Healpix_Map<float> &mapQ, const Healpix_Map<float> &mapU,
   PDT datatype);
template void write_Healpix_map_to_fits
  (fitshandle &out, const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ, const Healpix_Map<double> &mapU,
   PDT datatype);

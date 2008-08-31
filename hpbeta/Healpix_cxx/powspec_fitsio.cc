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

#include "powspec_fitsio.h"
#include "powspec.h"
#include "fitshandle.h"

using namespace std;

void read_powspec_from_fits (const string &infile, PowSpec &powspec,
  int nspecs, int lmax)
  {
  planck_assert ((nspecs==1)||(nspecs==4), "wrong number of spectra");
  fitshandle inp;
  inp.open(infile);
  inp.goto_hdu(2);
  if (inp.key_present("PDMTYPE"))    // official Planck format
    {
    inp.assert_pdmtype("POWERSPEC");
    arr<double> tmp;
    arr<int> itmp;
    if ((inp.coltype(1)==TINT)||(inp.coltype(1)==TLONG))
      inp.read_entire_column(1,itmp);
    else
      {
      cerr << "Warning: column containing l values is not of integer type!"
           << endl;
      inp.read_entire_column(1,tmp);
      itmp.alloc(tmp.size());
      for (int m=0; m<tmp.size(); ++m) itmp[m] = nearest<int>(tmp[m]);
      }

    if (nspecs==1)
      {
      arr<double> tt(lmax+1);
      tt.fill(0);
      inp.read_entire_column(2,tmp);
      for (int m=0; m<itmp.size(); ++m)
        if (itmp[m]<=lmax) tt[itmp[m]] = tmp[m];
      powspec.Set(tt);
      }
    else
      {
      arr<double> tt(lmax+1), gg(lmax+1), cc(lmax+1), tg(lmax+1);
      tt.fill(0); gg.fill(0); cc.fill(0); tg.fill(0);
      inp.read_entire_column(2,tmp);
      for (int m=0; m<itmp.size(); ++m)
        if (itmp[m]<=lmax) tt[itmp[m]] = tmp[m];
      inp.read_entire_column(3,tmp);
      for (int m=0; m<itmp.size(); ++m)
        if (itmp[m]<=lmax) gg[itmp[m]] = tmp[m];
      inp.read_entire_column(4,tmp);
      for (int m=0; m<itmp.size(); ++m)
        if (itmp[m]<=lmax) cc[itmp[m]] = tmp[m];
      inp.read_entire_column(5,tmp);
      for (int m=0; m<itmp.size(); ++m)
        if (itmp[m]<=lmax) tg[itmp[m]] = tmp[m];
      powspec.Set(tt,gg,cc,tg);
      }
    }
  else
    {
    int lmax_file = inp.nelems(1)-1;
    if (lmax_file<lmax)
      cerr << "warning: lmax in file smaller than expected; padding with 0."
           << endl;
    int lmax_read = min (lmax,lmax_file);
    if (nspecs==1)
      {
      arr<double> tt(lmax+1);
      inp.read_column_raw (1,&tt[0],lmax_read+1);
      for (int l=lmax_read+1; l<=lmax; ++l)
        tt[l]=0;
      powspec.Set(tt);
      }
    else
      {
      arr<double> tt(lmax+1), gg(lmax+1), cc(lmax+1), tg(lmax+1);
      inp.read_column_raw (1,&tt[0],lmax_read+1);
      inp.read_column_raw (2,&gg[0],lmax_read+1);
      inp.read_column_raw (3,&cc[0],lmax_read+1);
      inp.read_column_raw (4,&tg[0],lmax_read+1);
      for (int l=lmax_read+1; l<=lmax; ++l)
        tt[l]=gg[l]=cc[l]=tg[l]=0;
      powspec.Set(tt,gg,cc,tg);
      }
    }
  }

void write_powspec_to_fits (fitshandle &out,
  const PowSpec &powspec, int nspecs)
  {
  planck_assert ((nspecs==1)||(nspecs==4), "wrong number of spectra");
  vector<fitscolumn> cols;
  cols.push_back(fitscolumn("l","multipole order",1,TINT32BIT));
  cols.push_back(fitscolumn("Temperature C_l","Kelvin-squared",1,TDOUBLE));
  if (nspecs>1)
    {
    cols.push_back(fitscolumn("E-mode C_l","Kelvin-squared",1,TDOUBLE));
    cols.push_back(fitscolumn("B-mode C_l","Kelvin-squared",1,TDOUBLE));
    cols.push_back(fitscolumn("T-E cross-corr.","Kelvin-squared",1,TDOUBLE));
    }
  out.insert_bintab(cols);
  out.add_key("PDMTYPE",string("POWERSPEC"),"Planck data model type");
  arr<int> tmparr(powspec.Lmax()+1);
  for (int l=0; l<=powspec.Lmax(); ++l) tmparr[l]=l;
  out.write_column(1,tmparr);
  out.write_column(2,powspec.tt());
  if (nspecs>1)
    {
    out.write_column(3,powspec.gg());
    out.write_column(4,powspec.cc());
    out.write_column(5,powspec.tg());
    }
  }

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

#include "powspec_fitsio.h"
#include "powspec.h"
#include "fitshandle.h"

using namespace std;

void read_powspec_from_fits (const string &infile, PowSpec &powspec,
  int nspecs, int lmax)
  {
  planck_assert ((nspecs==1)||(nspecs==4)||(nspecs==6),
    "wrong number of spectra");
  fitshandle inp;
  inp.open(infile);
  inp.goto_hdu(2);

  arr<double> tt(lmax+1,0),gg(lmax+1,0),cc(lmax+1,0),tg(lmax+1,0),
              tc(lmax+1,0),gc(lmax+1,0);

  int lmax_file = safe_cast<int>(inp.nelems(1)-1);
  if (lmax_file<lmax)
    cerr << "warning: lmax in file smaller than expected; padding with 0."
          << endl;
  int lmax_read = min (lmax,lmax_file);
  inp.read_column_raw (1,&tt[0],lmax_read+1);
  if (nspecs>=4)
    {
    inp.read_column_raw (2,&gg[0],lmax_read+1);
    inp.read_column_raw (3,&cc[0],lmax_read+1);
    inp.read_column_raw (4,&tg[0],lmax_read+1);
    }
  if (nspecs==6)
    {
    inp.read_column_raw (5,&tc[0],lmax_read+1);
    inp.read_column_raw (6,&gc[0],lmax_read+1);
    }

  if (nspecs==1) powspec.Set(tt);
  if (nspecs==4) powspec.Set(tt,gg,cc,tg);
  if (nspecs==6) powspec.Set(tt,gg,cc,tg,tc,gc);
  }

void write_powspec_to_fits (fitshandle &out,
  const PowSpec &powspec, int nspecs)
  {
  planck_assert ((nspecs==1)||(nspecs==4)||(nspecs==6),
    "incorrect number of spectra");
  vector<fitscolumn> cols;
  cols.push_back(fitscolumn("Temperature C_l","unknown",1,PLANCK_FLOAT64));
  if (nspecs>1)
    {
    cols.push_back(fitscolumn("E-mode C_l","unknown",1,PLANCK_FLOAT64));
    cols.push_back(fitscolumn("B-mode C_l","unknown",1,PLANCK_FLOAT64));
    cols.push_back(fitscolumn("T-E cross-corr.","unknown",1,
      PLANCK_FLOAT64));
    }
  if (nspecs>4)
    {
    cols.push_back(fitscolumn("T-B cross-corr.","unknown",1,PLANCK_FLOAT64));
    cols.push_back(fitscolumn("E-B cross-corr.","unknown",1,PLANCK_FLOAT64));
    }
  out.insert_bintab(cols);
  out.write_column(1,powspec.tt());
  if (nspecs>1)
    {
    out.write_column(2,powspec.gg());
    out.write_column(3,powspec.cc());
    out.write_column(4,powspec.tg());
    }
  if (nspecs>4)
    {
    out.write_column(5,powspec.tc());
    out.write_column(6,powspec.gc());
    }
  }

void write_powspec_to_fits (const string &outfile,
  const PowSpec &powspec, int nspecs)
  {
  fitshandle out;
  out.create(outfile);
  write_powspec_to_fits(out,powspec,nspecs);
  }

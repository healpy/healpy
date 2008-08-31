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
 *  Copyright (C) 2003, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_data_io.h"
#include "arr.h"
#include "fitshandle.h"
#include "paramfile.h"
#include "simparams.h"

using namespace std;

void read_weight_ring (const string &dir, int nside, arr<double> &weight)
  {
  fitshandle inp;
  inp.open(dir+"/weight_ring_n"+intToString(nside,5)+".fits");
  inp.goto_hdu(2);
  weight.alloc (2*nside);
  inp.read_column(1,weight);
  }

void get_ring_weights (paramfile &params, simparams &par, int nside,
  arr<double> &weight)
  {
  bool weighted = params.find<bool>("weighted",false);
  par.add ("weighted","WEIGHTED",weighted,"ring weights used?");
  weight.alloc (2*nside);
  if (weighted)
    {
    string datadir = params.find<string>("healpix_data");
    read_weight_ring (datadir, nside, weight);
    for (int m=0; m<weight.size(); ++m) weight[m]+=1;
    }
  else
    weight.fill(1);
  }

void read_pixwin (const string &dir, int nside, arr<double> &temp)
  {
  fitshandle inp;
  inp.open(dir+"/pixel_window_n"+intToString(nside,4)+".fits");
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  }
void read_pixwin (const string &dir, int nside, arr<double> &temp,
  arr<double> &pol)
  {
  fitshandle inp;
  inp.open(dir+"/pixel_window_n"+intToString(nside,4)+".fits");
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  if (pol.size()==0)
    inp.read_entire_column(2,pol);
  else
    inp.read_column(2,pol);
  }

void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin)
  {
  bool do_pixwin = params.find<bool>("pixel_window",false);
  par.add("pixel_window","PIXWIN",do_pixwin,"pixel window used?");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  if (do_pixwin)
    {
    string datadir = params.find<string>("healpix_data");
    read_pixwin (datadir,nside,pixwin);
    }
  }
void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin, arr<double> &pixwin_pol)
  {
  bool do_pixwin = params.find<bool>("pixel_window",false);
  par.add("pixel_window","PIXWIN",do_pixwin,"pixel window used?");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  pixwin_pol.alloc(lmax+1);
  pixwin_pol.fill(1);
  if (do_pixwin)
    {
    string datadir = params.find<string>("healpix_data");
    read_pixwin (datadir,nside,pixwin,pixwin_pol);
    }
  }

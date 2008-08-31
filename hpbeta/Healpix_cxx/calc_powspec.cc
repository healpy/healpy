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
 *  Copyright (C) 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "fitshandle.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "alm_powspec_tools.h"

using namespace std;

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  announce ("calc_powspec");
  planck_assert (argc==3||argc==4,
    "usage: calc_powspec <almfile1> [<almfile2>] <powspec_file>");

  if (argc==3)
    {
    fitshandle inp;
    inp.open(argv[1]);
    inp.goto_hdu(2);
    int lmax,mmax;
    get_almsize(inp,lmax,mmax);
    Alm<xcomplex<float> > alm;
    read_Alm_from_fits (inp,alm,lmax,mmax);
    PowSpec powspec;
    extract_powspec (alm,powspec);
    fitshandle out;
    out.create (argv[2]);
    write_powspec_to_fits (out,powspec,1);
    }
  else
    {
    fitshandle inp;
    inp.open(argv[1]);
    inp.goto_hdu(2);
    int lmax,mmax;
    get_almsize(inp,lmax,mmax);
    Alm<xcomplex<float> > alm1;
    read_Alm_from_fits (inp,alm1,lmax,mmax);
    inp.close();
    inp.open(argv[2]);
    inp.goto_hdu(2);
    get_almsize(inp,lmax,mmax);
    Alm<xcomplex<float> > alm2;
    read_Alm_from_fits (inp,alm2,lmax,mmax);
    inp.close();
    PowSpec powspec;
    extract_crosspowspec (alm1,alm2,powspec);
    fitshandle out;
    out.create (argv[3]);
    write_powspec_to_fits (out,powspec,1);
    }
PLANCK_DIAGNOSIS_END
  }

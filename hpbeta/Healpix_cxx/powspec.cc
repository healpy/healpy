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

#include "powspec.h"
#include "lsconstants.h"

using namespace std;

void PowSpec::dealloc()
  {
  tt_.dealloc();
  gg_.dealloc();
  cc_.dealloc();
  tg_.dealloc();
  tc_.dealloc();
  gc_.dealloc();
  }

PowSpec::PowSpec(int nspec, int lmax)
  { Set(nspec,lmax); }

PowSpec::~PowSpec()
  {}

void PowSpec::assertArraySizes() const
  {
  planck_assert((num_specs==1) || (num_specs==4) || (num_specs==6),
    "incorrect number of spectral components");
  if (num_specs==1)
    planck_assert(multiequal(tsize(0),gg_.size(),cc_.size(),tg_.size(),
      tc_.size(),gc_.size()), "incorrect array sizes");
  if (num_specs==4)
    {
    planck_assert(multiequal(tt_.size(),gg_.size(),cc_.size(),tg_.size()),
      "incorrect array sizes");
    planck_assert(multiequal(tsize(0),tc_.size(),gc_.size()),
      "incorrect array sizes");
    }
  if (num_specs==6)
    planck_assert(multiequal(tt_.size(),gg_.size(),cc_.size(),tg_.size(),
      tc_.size(),gc_.size()), "incorrect array sizes");
  }

bool PowSpec::consistentAutoPowspec() const
  {
  for (tsize l=0; l<tt_.size(); ++l)
    if (tt_[l]<0) return false;
  if (num_specs>=4)
    for (tsize l=0; l<tt_.size(); ++l)
      {
      if (gg_[l]<0) return false;
      if (cc_[l]<0) return false;
      if (abs(tg_[l])>sqrt(tt_[l]*gg_[l])) return false;
      }
  if (num_specs==6)
    for (tsize l=0; l<tt_.size(); ++l)
      {
      if (abs(tc_[l])>sqrt(tt_[l]*cc_[l])) return false;
      if (abs(gc_[l])>sqrt(gg_[l]*cc_[l])) return false;
      }
  return true;
  }

void PowSpec::Set(int nspec, int lmax)
  {
  num_specs=nspec;
  planck_assert ((num_specs==1) || (num_specs==4) || (num_specs==6),
    "wrong number of spectrums");
  tt_.alloc(lmax+1);
  if (num_specs>1)
    {
    gg_.alloc(lmax+1);
    cc_.alloc(lmax+1);
    tg_.alloc(lmax+1);
    }
  if (num_specs>4)
    {
    tc_.alloc(lmax+1);
    gc_.alloc(lmax+1);
    }
  }

void PowSpec::Set(arr<double> &tt_new)
  {
  dealloc();
  num_specs = 1;
  tt_.transfer(tt_new);
//FIXME: temporarily relaxed to allow cross-spectra
  if (!consistentAutoPowspec())
    cerr << "Warning: negative values in TT spectrum" << endl;
  }

void PowSpec::Set(arr<double> &tt_new, arr<double> &gg_new,
  arr<double> &cc_new, arr<double> &tg_new)
  {
  dealloc();
  num_specs = 4;
  tt_.transfer(tt_new); gg_.transfer(gg_new);
  cc_.transfer(cc_new); tg_.transfer(tg_new);
  assertArraySizes();
  }

void PowSpec::Set(arr<double> &tt_new, arr<double> &gg_new, arr<double> &cc_new,
  arr<double> &tg_new, arr<double> &tc_new, arr<double> &gc_new)
  {
  Set (tt_new, gg_new, cc_new, tg_new);
  num_specs = 6;
  tc_.transfer(tc_new); gc_.transfer(gc_new);
  assertArraySizes();
  }

void PowSpec::smoothWithGauss (double fwhm)
  {
  planck_assert (num_specs<=4, "not yet implemented for num_specs>4");
  double sigma = fwhm*fwhm2sigma;
  double fact_pol = exp(2*sigma*sigma);
  for (tsize l=0; l<tt_.size(); ++l)
    {
    double f1 = exp(-.5*l*(l+1)*sigma*sigma);
    double f2 = f1*fact_pol;
    tt_[l] *= f1*f1;
    if (num_specs>1)
      {
      gg_[l] *= f2*f2;
      cc_[l] *= f2*f2;
      tg_[l] *= f1*f2;
      }
    }
  }

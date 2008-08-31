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
 *  Copyright (C) 2003,2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "powspec.h"

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

void PowSpec::Set(arr<double> &tt_new)
  {
  dealloc();
  num_specs = 1;
  tt_.transfer(tt_new);
//FIXME: temporarily relaxed to allow cross-spectra
  for (int l=0; l<tt_.size(); ++l)
    if (tt_[l]<0)
      {
      cerr << "Warning: negative values in TT spectrum" << endl;
      break;
      }
  }

void PowSpec::Set(arr<double> &tt_new, arr<double> &gg_new,
  arr<double> &cc_new, arr<double> &tg_new)
  {
  dealloc();
  num_specs = 4;
  tt_.transfer(tt_new);
  gg_.transfer(gg_new);
  cc_.transfer(cc_new);
  tg_.transfer(tg_new);
  planck_assert((tt_.size()==gg_.size()) && (tt_.size()==cc_.size())
    && (tt_.size()==tg_.size()), "PowSpec::Set: size mismatch");
  for (int l=0; l<tt_.size(); ++l)
    {
    planck_assert (tt_[l]>=0, "negative TT spectrum at l="+dataToString(l));
    planck_assert (gg_[l]>=0, "negative GG spectrum at l="+dataToString(l));
    planck_assert (cc_[l]>=0, "negative CC spectrum at l="+dataToString(l));
    planck_assert (abs(tg_[l]<=sqrt(tt_[l]*gg_[l])),
      "Inconsistent T, E and TxE terms at l="+dataToString(l));
    }
  }

void PowSpec::Smooth_with_Gauss (double fwhm)
  {
  planck_assert (num_specs<=4, "not yet implemented for num_specs>4");
  double sigma = fwhm*fwhm2sigma;
  double fact_pol = exp(2*sigma*sigma);
  for (int l=0; l<tt_.size(); ++l)
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

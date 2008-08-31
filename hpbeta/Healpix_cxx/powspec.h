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

#ifndef POWSPEC_H
#define POWSPEC_H

#include "arr.h"

/*! Class for storing unpolarised and polarised power spectra. */
class PowSpec
  {
  private:
    arr<double> tt_, gg_, cc_, tg_, tc_, gc_;
    int num_specs;

    void dealloc();

  public:
    /*! */
    PowSpec() {}
    /*! Constructs a \a PowSpec with \a nspec components and a maximum
        multipole of \a lmax. \a nspec can be 1 (TT), 4 (TT,GG,CC,TG) or
        6 (TT,GG,CC,TG,TC,GC). */
    PowSpec(int nspec, int lmax)
      : num_specs(nspec)
      {
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

    /*! Returns the number of spectral components. */
    int Num_specs() const { return num_specs; }
    /*! Returns the maximum \a l. */
    int Lmax() const { return tt_.size()-1; }
    /*! Returns the TT array (read-only). */
    const arr<double> &tt() const { return tt_; }
    /*! Returns the GG array (read-only). */
    const arr<double> &gg() const { return gg_; }
    /*! Returns the CC array (read-only). */
    const arr<double> &cc() const { return cc_; }
    /*! Returns the TG array (read-only). */
    const arr<double> &tg() const { return tg_; }
    /*! Returns the TC array (read-only). */
    const arr<double> &tc() const { return tc_; }
    /*! Returns the GC array (read-only). */
    const arr<double> &gc() const { return gc_; }
    /*! Returns TT(l) (read-write). */
    double &tt (int l) { return tt_[l]; }
    /*! Returns GG(l) (read-write). */
    double &gg (int l) { return gg_[l]; }
    /*! Returns CC(l) (read-write). */
    double &cc (int l) { return cc_[l]; }
    /*! Returns TG(l) (read-write). */
    double &tg (int l) { return tg_[l]; }
    /*! Returns TC(l) (read-write). */
    double &tc (int l) { return tc_[l]; }
    /*! Returns GC(l) (read-write). */
    double &gc (int l) { return gc_[l]; }
    /*! Returns TT(l) (read-only). */
    const double &tt (int l) const { return tt_[l]; }
    /*! Returns GG(l) (read-only). */
    const double &gg (int l) const { return gg_[l]; }
    /*! Returns CC(l) (read-only). */
    const double &cc (int l) const { return cc_[l]; }
    /*! Returns TG(l) (read-only). */
    const double &tg (int l) const { return tg_[l]; }
    /*! Returns TC(l) (read-only). */
    const double &tc (int l) const { return tc_[l]; }
    /*! Returns GC(l) (read-only). */
    const double &gc (int l) const { return gc_[l]; }

    /*! Sets the whole TT array. */
    void Set(arr<double> &tt_new);
    /*! Sets all components. */
    void Set(arr<double> &tt_new, arr<double> &gg_new,
             arr<double> &cc_new, arr<double> &tg_new);
    /* Smooths the spectrum with a Gaussian beam.
       \a fwhm is given in radian.
       \note This is only implememted for 1 and 4 spectra so far. */
    void Smooth_with_Gauss (double fwhm);
  };

#endif

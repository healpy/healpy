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

/*! \file healpix_map.h
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_MAP_H
#define HEALPIX_MAP_H

#include "arr.h"
#include "healpix_base.h"

//! Healpix value representing "undefined"
const double Healpix_undef=-1.6375e30;

/*! A HEALPix map of a given datatype */
template<typename T> class Healpix_Map: public Healpix_Base
  {
  private:
    arr<T> map;

  public:
    /*! Constructs an unallocated map. */
    Healpix_Map () {}
    /*! Constructs a map with a given \a order and the ordering
        scheme \a scheme. */
    Healpix_Map (int order, Healpix_Ordering_Scheme scheme)
      : Healpix_Base (order, scheme), map(npix_) {}
    /*! Constructs a map with a given \a nside and the ordering
        scheme \a scheme. */
    Healpix_Map (int nside, Healpix_Ordering_Scheme scheme, const nside_dummy)
      : Healpix_Base (nside, scheme, SET_NSIDE), map(npix_) {}
    /*! Constructs a map from the contents of \a data and sets the ordering
        scheme to \a Scheme. The size of \a data must be a valid HEALPix
        map size. */
    Healpix_Map (const arr<T> &data, Healpix_Ordering_Scheme scheme)
      : Healpix_Base (npix2nside(data.size()), scheme, SET_NSIDE), map(data) {}

    /*! Deletes the old map, creates a map from the contents of \a data and
        sets the ordering scheme to \a scheme. The size of \a data must be a
        valid HEALPix map size.
        \note On exit, \a data is zero-sized! */
    void Set (arr<T> &data, Healpix_Ordering_Scheme scheme)
      {
      Healpix_Base::SetNside(npix2nside (data.size()), scheme);
      map.transfer(data);
      }

    /*! Deletes the old map and creates a new map  with a given \a order
        and the ordering scheme \a scheme. */
    void Set (int order, Healpix_Ordering_Scheme scheme)
      {
      Healpix_Base::Set(order, scheme);
      map.alloc(npix_);
      }
    /*! Deletes the old map and creates a new map  with a given \a nside
        and the ordering scheme \a scheme. */
    void SetNside (int nside, Healpix_Ordering_Scheme scheme)
      {
      Healpix_Base::SetNside(nside, scheme);
      map.alloc(npix_);
      }

    /*! Fills the map with \a val. */
    void fill (const T &val)
      { map.fill(val); }

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme. \a orig must have the same resolution as the
        current map. */
    void Import_nograde (const Healpix_Map<T> &orig)
      {
      planck_assert (nside_==orig.nside_,
        "Import_nograde: maps have different nside");
      if (orig.scheme_ == scheme_)
        for (int m=0; m<npix_; ++m) map[m] = orig.map[m];
      else
        {
        swapfunc swapper = (scheme_ == NEST) ?
          &Healpix_Base::ring2nest : &Healpix_Base::nest2ring;
#pragma omp parallel
{
        int m;
#pragma omp for schedule (dynamic,5000)
        for (m=0; m<npix_; ++m) map[(this->*swapper)(m)] = orig.map[m];
}
        }
      }

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme and the map resolution. \a orig must have higher
        resolution than the current map. */
    void Import_upgrade (const Healpix_Map<T> &orig)
      {
      planck_assert(nside_>orig.nside_,"Import_upgrade: this is no upgrade");
      int fact = nside_/orig.nside_;
      planck_assert (nside_==orig.nside_*fact,
        "the larger Nside must be a multiple of the smaller one");

#pragma omp parallel
{
      int m;
#pragma omp for schedule (dynamic,5000)
      for (m=0; m<orig.npix_; ++m)
        {
        int x,y,f;
        orig.pix2xyf(m,x,y,f);
        for (int j=fact*y; j<fact*(y+1); ++j)
          for (int i=fact*x; i<fact*(x+1); ++i)
            {
            int mypix = xyf2pix(i,j,f);
            map[mypix] = orig.map[m];
            }
        }
}
      }

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme and the map resolution. \a orig must have higher
        resolution than the current map.
        \a pessimistic determines whether or not
        pixels are set to \a Healpix_undef when not all of the corresponding
        high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
    void Import_degrade (const Healpix_Map<T> &orig, bool pessimistic=false);

    /*! Imports the map \a orig into the current map, adjusting the
        ordering scheme and the map resolution if necessary.
        When downgrading, \a pessimistic determines whether or not
        pixels are set to \a Healpix_undef when not all of the corresponding
        high-resolution pixels are defined.

        This method is instantiated for \a float and \a double only. */
    void Import (const Healpix_Map<T> &orig, bool pessimistic=false)
      {
      if (orig.nside_ == nside_) // no up/degrading
        Import_nograde(orig);
      else if (orig.nside_ < nside_) // upgrading
        Import_upgrade(orig);
      else
        Import_degrade(orig,pessimistic);
      }

    /*! Returns a constant reference to the pixel with the number \a pix. */
    const T &operator[] (int pix) const { return map[pix]; }
    /*! Returns a reference to the pixel with the number \a pix. */
    T &operator[] (int pix) { return map[pix]; }

    /*! Swaps the map ordering from RING to NEST and vice versa.
        This is done in-place (i.e. with negligible space overhead). */
    void swap_scheme()
      {
      swapfunc swapper = (scheme_ == NEST) ?
        &Healpix_Base::ring2nest : &Healpix_Base::nest2ring;

      arr<int> cycle=swap_cycles();

      for (tsize m=0; m<cycle.size(); ++m)
        {
        int istart = cycle[m];

        T pixbuf = map[istart];
        int iold = istart, inew = (this->*swapper)(istart);
        while (inew != istart)
          {
          map[iold] = map[inew];
          iold = inew;
          inew = (this->*swapper)(inew);
          }
        map[iold] = pixbuf;
        }
      scheme_ = (scheme_==RING) ? NEST : RING;
      }

    /*! performs the actual interpolation using \a pix and \a wgt. */
    T interpolation (const fix_arr<int,4> &pix,
      const fix_arr<double,4> &wgt) const
      {
      double wtot=0;
      T res=T(0);
      for (tsize i=0; i<4; ++i)
        {
        T val=map[pix[i]];
        if (!approx<double>(val,Healpix_undef))
          { res+=T(val*wgt[i]); wtot+=wgt[i]; }
        }
      return (wtot==0.) ? T(Healpix_undef) : T(res/wtot);
      }
    /*! Returns the interpolated map value at \a ptg */
    T interpolated_value (const pointing &ptg) const
      {
      fix_arr<int,4> pix;
      fix_arr<double,4> wgt;
      get_interpol (ptg, pix, wgt);
      return interpolation (pix, wgt);
      }

    /*! Returns a constant reference to the map data. */
    const arr<T> &Map() const { return map; }

    /*! Returns the minimum and maximum value of the map in
        \a Min and \a Max.

        This method is instantiated for \a float and \a double only. */
    void minmax (T &Min, T &Max) const;

    /*! Swaps the contents of two Healpix_Map objects. */
    void swap (Healpix_Map &other)
      {
      Healpix_Base::swap(other);
      map.swap(other.map);
      }

    /*! Returns the average of all defined map pixels. */
    double average() const
      {
      double avg=0;
      int pix=0;
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { ++pix; avg+=map[m]; }
      return (pix>0) ? avg/pix : Healpix_undef;
      }

    /*! Adds \a val to all defined map pixels. */
    void Add (T val)
      {
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { map[m]+=val; }
      }

    /*! Multiplies all defined map pixels by \a val. */
    void Scale (T val)
      {
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { map[m]*=val; }
      }

    /*! Returns the root mean square of the map, not counting undefined
        pixels. */
    double rms() const
      {
      using namespace std;

      double result=0;
      int pix=0;
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { ++pix; result+=map[m]*map[m]; }
      return (pix>0) ? sqrt(result/pix) : Healpix_undef;
      }
    /*! Returns the maximum absolute value in the map, ignoring undefined
        pixels. */
    T absmax() const
      {
      using namespace std;

      T result=0;
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { result = max(result,abs(map[m])); }
      return result;
      }
    /*! Returns \a true, if no pixel has the value \a Healpix_undef,
        else \a false. */
    bool fullyDefined() const
      {
      for (int m=0; m<npix_; ++m)
        if (approx<double>(map[m],Healpix_undef))
          return false;
      return true;
      }
    /*! Sets all pixels with the value \a Healpix_undef to 0, and returns
        the number of modified pixels. */
    tsize replaceUndefWith0()
      {
      tsize res=0;
      for (int m=0; m<npix_; ++m)
        if (approx<double>(map[m],Healpix_undef))
          { map[m]=0.; ++res; }
      return res;
      }
  };

#endif

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
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_MAP_H
#define HEALPIX_MAP_H

#include "arr.h"
#include "healpix_base.h"

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
        valid HEALPix map size. */
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
      pix2xyf to_xyf = (orig.scheme_==RING) ?
        &Healpix_Map::ring2xyf : &Healpix_Map::nest2xyf;
      xyf2pix from_xyf = (scheme_==RING) ?
        &Healpix_Map::xyf2ring : &Healpix_Map::xyf2nest;

#pragma omp parallel
{
      int m;
#pragma omp for schedule (dynamic,5000)
      for (m=0; m<orig.npix_; ++m)
        {
        int x,y,f;
        (orig.*to_xyf)(m,x,y,f);
        for (int j=fact*y; j<fact*(y+1); ++j)
          for (int i=fact*x; i<fact*(x+1); ++i)
            {
            int mypix = (this->*from_xyf)(i,j,f);
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
      static const int clen[] = { 0,7,5,4,12,10,13,18,14,19,18,17,27,21 };
      static const int cycle[][30] = {
        { },
        { 0,1,8,12,16,21,40 },
        { 0,1,2,40,114 },
        { 0,4,160,263 },
        { 0,4,30,49,51,87,526,1027,1105,1387,1807,2637 },
        { 0,8,10,18,39,74,146,307,452,4737 },
        { 0,1,2,7,9,17,80,410,1526,1921,32859,33566,38931 },
        { 0,5,6,10,12,24,27,95,372,494,924,1409,3492,4248,9137,66043,103369,
          156899 },
        { 0,1,2,3,4,45,125,351,697,24337,102940,266194,341855,419857 },
        { 0,1,2,3,9,16,1705,2082,2126,8177,12753,15410,52642,80493,83235,
          88387,99444,1675361,2495125 },
        { 0,2,6,8,9,11,20,50,93,152,183,2137,13671,44794,486954,741908,
          4803258,5692573 },
        { 0,1,5,6,44,53,470,2847,3433,4906,13654,14710,400447,1797382,
          2744492,18775974,23541521 },
        { 0,4,9,10,16,33,83,117,318,451,5759,10015,128975,171834,211256,
          347608,1278690,2154097,2590798,3427694,5581717,21012301,27023976,
          72522811,95032729,139166747,171822389 },
        { 0,5,10,267,344,363,2968,3159,9083,18437,76602,147614,1246902,
          1593138,2035574,6529391,9511830,11340287,29565945,281666026,
          677946848 } };

      swapfunc swapper = (scheme_ == NEST) ?
        &Healpix_Base::ring2nest : &Healpix_Base::nest2ring;

      planck_assert (order_>=0, "swap_scheme(): need hierarchical map");

      for (int m=0; m<clen[order_]; ++m)
        {
        int istart = cycle[order_][m];

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
      return map[pix[0]]*wgt[0] + map[pix[1]]*wgt[1]
           + map[pix[2]]*wgt[2] + map[pix[3]]*wgt[3];
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
      return avg/pix;
      }

    /*! Adds \a val to all defined map pixels. */
    void add (T val)
      {
      for (int m=0; m<npix_; ++m)
        if (!approx<double>(map[m],Healpix_undef))
          { map[m]+=val; }
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
      return sqrt(result/pix);
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
  };

#endif

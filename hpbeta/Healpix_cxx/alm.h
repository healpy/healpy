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

/*! \file alm.h
 *  Class for storing spherical harmonic coefficients.
 *
 *  Copyright (C) 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALM_H
#define PLANCK_ALM_H

#include "arr.h"

/*! Class for storing spherical harmonic coefficients. */
template<typename T> class Alm
  {
  private:
    int lmax, mmax, tval;
    arr<T> alm;

  public:
    /*! Returns the number of coefficients in an Alm object with maximum
        quantum numbers \a l and \a m. */
    static long Num_Alms (int l, int m)
      {
      planck_assert(m<=l,"mmax must not be larger than lmax");
      return ((m+1)*(m+2))/2 + (m+1)*(l-m);
      }

    /*! Constructs an Alm object with given \a lmax and \a mmax. */
    Alm (int lmax_=0, int mmax_=0)
      : lmax(lmax_), mmax(mmax_), tval(2*lmax+1),
        alm (Num_Alms(lmax,mmax))
      {}

    /*! Deletes the old coefficients and allocates storage according to
        \a lmax and \a mmax. */
    void Set (int lmax_, int mmax_)
      {
      lmax=lmax_;
      mmax=mmax_;
      tval=2*lmax+1;
      alm.alloc(Num_Alms(lmax,mmax));
      }

    /*! Deallocates the old coefficients and uses the content of \a data
        for storage. \a data is deallocated during the call. */
    void Set (arr<T> &data, int lmax_, int mmax_)
      {
      planck_assert (Num_Alms(lmax_,mmax_)==data.size(),"wrong array size");
      lmax=lmax_;
      mmax=mmax_;
      tval=2*lmax+1;
      alm.transfer(data);
      }

    /*! Sets all coefficients to zero. */
    void SetToZero ()
      { alm.fill (0); }

    /*! Multiplies all coefficients by \a factor. */
    template<typename T2> void Scale (const T2 &factor)
      { for (int m=0; m<alm.size(); ++m) alm[m]*=factor; }
    /*! \a a(l,m) *= \a factor[l] for all \a l,m. */
    template<typename T2> void ScaleL (const arr<T2> &factor)
      {
      planck_assert(factor.size()>lmax, "alm.ScaleL: factor array too short");
      for (int m=0; m<=mmax; ++m)
        for (int l=m; l<=lmax; ++l)
          operator()(l,m)*=factor[l];
      }
    /*! Adds \a num to a_00. */
    template<typename T2> void Add (const T2 &num)
      { alm[0]+=num; }

    /*! Returns a reference to the specified coefficient. */
    T &operator() (int l, int m)
      { return alm[((m*(tval-m))>>1) + l]; }
    /*! Returns a constant reference to the specified coefficient. */
    const T &operator() (int l, int m) const
      { return alm[((m*(tval-m))>>1) + l]; }

    /*! Returns a pointer for a given m, from which the address of a_lm
        can be obtained by adding l. */
    T *mstart (int m)
      { return &alm[(m*(tval-m))>>1]; }
    /*! Returns a pointer for a given m, from which the address of a_lm
        can be obtained by adding l. */
    const T *mstart (int m) const
      { return &alm[(m*(tval-m))>>1]; }

    /*! Returns the maximum \a l */
    int Lmax() const { return lmax; }
    /*! Returns the maximum \a m */
    int Mmax() const { return mmax; }

    /*! Returns a constant reference to the a_lm data. */
    const arr<T> &Alms () const { return alm; }

    /*! Swaps the contents of two A_lm objects. */
    void swap (Alm &other)
      {
      std::swap(lmax, other.lmax);
      std::swap(mmax, other.mmax);
      std::swap(tval, other.tval);
      alm.swap(other.alm);
      }

    /*! Returns \a true, if both objects have the same \a lmax and \a mmax,
        else  \a false. */
    bool conformable (const Alm &other) const
      { return ((lmax==other.lmax) && (mmax==other.mmax)); }

    /*! Adds all coefficients from \a other to the own coefficients. */
    void Add (const Alm &other)
      {
      planck_assert (conformable(other), "A_lm are not conformable");
      for (int m=0; m<alm.size(); ++m)
        alm[m] += other.alm[m];
      }
  };

#endif

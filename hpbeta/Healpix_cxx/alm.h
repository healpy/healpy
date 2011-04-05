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
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALM_H
#define PLANCK_ALM_H

#include "arr.h"

/*! Base class for calculating the storage layout of spherical harmonic
    coefficients. */
class Alm_Base
  {
  protected:
    int lmax, mmax, tval;

  public:
    /*! Returns the total number of coefficients for maximum quantum numbers
        \a l and \a m. */
    static tsize Num_Alms (int l, int m);

    /*! Constructs an Alm_Base object with given \a lmax and \a mmax. */
    Alm_Base (int lmax_=0, int mmax_=0)
      : lmax(lmax_), mmax(mmax_), tval(2*lmax+1) {}

    /*! Changes the object's maximum quantum numbers to \a lmax and \a mmax. */
    void Set (int lmax_, int mmax_)
      {
      lmax=lmax_;
      mmax=mmax_;
      tval=2*lmax+1;
      }

    /*! Returns the maximum \a l */
    int Lmax() const { return lmax; }
    /*! Returns the maximum \a m */
    int Mmax() const { return mmax; }

    /*! Returns an array index for a given m, from which the index of a_lm
        can be obtained by adding l. */
    int index_l0 (int m) const
      { return ((m*(tval-m))>>1); }

    /*! Returns the array index of the specified coefficient. */
    int index (int l, int m) const
      { return index_l0(m) + l; }

    /*! Returns \a true, if both objects have the same \a lmax and \a mmax,
        else  \a false. */
    bool conformable (const Alm_Base &other) const
      { return ((lmax==other.lmax) && (mmax==other.mmax)); }

    /*! Swaps the contents of two Alm_Base objects. */
    void swap (Alm_Base &other);
  };

/*! Class for storing spherical harmonic coefficients. */
template<typename T> class Alm: public Alm_Base
  {
  private:
    arr<T> alm;

  public:
    /*! Constructs an Alm object with given \a lmax and \a mmax. */
    Alm (int lmax_=0, int mmax_=0)
      : Alm_Base(lmax_,mmax_), alm (Num_Alms(lmax,mmax)) {}

    /*! Deletes the old coefficients and allocates storage according to
        \a lmax and \a mmax. */
    void Set (int lmax_, int mmax_)
      {
      Alm_Base::Set(lmax_, mmax_);
      alm.alloc(Num_Alms(lmax,mmax));
      }

    /*! Deallocates the old coefficients and uses the content of \a data
        for storage. \a data is deallocated during the call. */
    void Set (arr<T> &data, int lmax_, int mmax_)
      {
      planck_assert (Num_Alms(lmax_,mmax_)==data.size(),"wrong array size");
      Alm_Base::Set(lmax_, mmax_);
      alm.transfer(data);
      }

    /*! Sets all coefficients to zero. */
    void SetToZero ()
      { alm.fill (0); }

    /*! Multiplies all coefficients by \a factor. */
    template<typename T2> void Scale (const T2 &factor)
      { for (tsize m=0; m<alm.size(); ++m) alm[m]*=factor; }
    /*! \a a(l,m) *= \a factor[l] for all \a l,m. */
    template<typename T2> void ScaleL (const arr<T2> &factor)
      {
      planck_assert(factor.size()>tsize(lmax),
        "alm.ScaleL: factor array too short");
      for (int m=0; m<=mmax; ++m)
        for (int l=m; l<=lmax; ++l)
          operator()(l,m)*=factor[l];
      }
    /*! Adds \a num to a_00. */
    template<typename T2> void Add (const T2 &num)
      { alm[0]+=num; }

    /*! Returns a reference to the specified coefficient. */
    T &operator() (int l, int m)
      { return alm[index(l,m)]; }
    /*! Returns a constant reference to the specified coefficient. */
    const T &operator() (int l, int m) const
      { return alm[index(l,m)]; }

    /*! Returns a pointer for a given m, from which the address of a_lm
        can be obtained by adding l. */
    T *mstart (int m)
      { return &alm[index_l0(m)]; }
    /*! Returns a pointer for a given m, from which the address of a_lm
        can be obtained by adding l. */
    const T *mstart (int m) const
      { return &alm[index_l0(m)]; }

    /*! Returns a constant reference to the a_lm data. */
    const arr<T> &Alms () const { return alm; }

    /*! Swaps the contents of two Alm objects. */
    void swap (Alm &other)
      {
      Alm_Base::swap(other);
      alm.swap(other.alm);
      }

    /*! Adds all coefficients from \a other to the own coefficients. */
    void Add (const Alm &other)
      {
      planck_assert (conformable(other), "A_lm are not conformable");
      for (tsize m=0; m<alm.size(); ++m)
        alm[m] += other.alm[m];
      }
  };

#endif

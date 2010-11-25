/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file planck_rng.h
 *  This file contains the random number generator
 *  used by the Planck LevelS package.
 *  The generator is a C++ port of the xorshift generator xor128() described
 *  in Marsaglia, Journal of Statistical Software 2003, vol 8.
 *  It has a period of 2^128 - 1.
 *
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_RNG_H
#define PLANCK_RNG_H

#include <cmath>
#include "error_handling.h"

/*! C++ port of the xorshift generator xor128() described in Marsaglia,
    Journal of Statistical Software 2003, vol 8.
    It has a period of 2^128 - 1. */
class planck_rng
  {
  private:
    unsigned int x,y,z,w;
    double small, gset;
    bool empty;

    void twiddle (unsigned int &v)
      {
      for (int i=0; i<9; ++i)
        {
        v ^= v<<13;
        v ^= v>>17;
        v ^= v<<5;
        }
      }

    void init_rng ()
      {
      // avoid zero seeds
      if (x==0) x = 123456789;
      if (y==0) y = 362436069;
      if (z==0) z = 521288629;
      if (w==0) w = 88675123;

      // shuffle the bits of the seeds
      twiddle(x); twiddle(y); twiddle(z); twiddle(w);

      // burn in the RNG
      for (int i=0; i<16; ++i)
        int_rand_uni();
      }

  public:
    /*! Initializes the generator with 0 to 4 seed values. */
    planck_rng (unsigned int x1=123456789, unsigned int y1=362436069,
                unsigned int z1=521288629, unsigned int w1=88675123)
      : x(x1), y(y1), z(z1), w(w1),
        small(1./(1.+double(0xFFFFFFFF))), gset(0.), empty(true)
      {
      planck_assert (sizeof(unsigned int)==4, "wrong integer size for RNG");
      init_rng();
      }

    /*! Re-initializes the generator with 0 to 4 seed values. */
    void seed (unsigned int x1=123456789, unsigned int y1=362436069,
      unsigned int z1=521288629, unsigned int w1=88675123)
      {
      x = x1; y = y1; z = z1; w = w1;
      empty = true;
      init_rng();
      }

    /*! Returns uniformly distributed random integer numbers from the
        interval [0;0xFFFFFFFF]. */
    unsigned int int_rand_uni()
      {
      unsigned int t = x^(x<<11);
      x = y;
      y = z;
      z = w;

      return w=(w^(w>>19))^(t^(t>>8));
      }

    //! Returns uniformly distributed random numbers from the interval [0;1[.
    double rand_uni()
      {
      return small*int_rand_uni();
      }

    //! Returns random numbers with Gaussian distribution (mean=0, sigma=1).
    /*! Uses rand_uni() internally. */
    double rand_gauss()
      {
      using namespace std;
      if (empty)
        {
        double v1,v2,rsq;
        do
          {
          v1=2*rand_uni()-1.;
          v2=2*rand_uni()-1.;
          rsq=v1*v1+v2*v2;
          }
        while ((rsq>=1) || (rsq==0));
        double fac=sqrt(-2*log(rsq)/rsq);
        gset=v1*fac;
        empty=false;
        return v2*fac;
        }
      else
        {
        empty=true;
        return gset;
        }
      }

    //! Returns exponentially distributed random numbers (mean=1, nonnegative)
    /*! Uses rand_uni() internally. */
    double rand_exp()
      {
      using namespace std;
      double val=rand_uni();
      if (val==0.) val=1.;
      return -log(val);
      }
  };

#endif

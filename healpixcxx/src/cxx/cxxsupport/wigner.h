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

/*! \file wigner.h
 *  Several C++ classes for calculating Wigner matrices
 *
 *  Copyright (C) 2009-2011 Max-Planck-Society
 *  \author Martin Reinecke and others (see individual classes)
 */

#ifndef PLANCK_WIGNER_H
#define PLANCK_WIGNER_H

#include <cmath>
#include "arr.h"

#include "sse_utils_cxx.h"

/*! Class for calculation of the Wigner matrix at pi/2, using Risbo recursion
    in a way that cannot easily be parallelised, but is fairly efficient on
    scalar machines. */
class wigner_d_halfpi_risbo_scalar
  {
  private:
    double pq;
    arr<double> sqt;
    arr2<double> d;
    int n;

    void do_line0 (double *l1, int j);
    void do_line (const double *l1, double *l2, int j, int k);

  public:
    wigner_d_halfpi_risbo_scalar(int lmax);

    const arr2<double> &recurse ();
  };

/*! Class for calculation of the Wigner matrix at arbitrary angles, using Risbo
    recursion in a way that cannot easily be parallelised, but is fairly
    efficient on scalar machines. */
class wigner_d_risbo_scalar
  {
  private:
    double p,q;
    arr<double> sqt;
    arr2<double> d;
    int n;

    void do_line0 (double *l1, int j);
    void do_line (const double *l1, double *l2, int j, int k);

  public:
    wigner_d_risbo_scalar(int lmax, double ang);

    const arr2<double> &recurse ();
  };

/*! Class for calculation of the Wigner matrix at pi/2, using Risbo recursion
    in a way that can be OpenMP-parallelised. This approach uses more memory
    and is slightly slower than wigner_d_halfpi_risbo_scalar. */
class wigner_d_halfpi_risbo_openmp
  {
  private:
    double pq;
    arr<double> sqt;
    arr2<double> d,dd;
    int n;

  public:
    wigner_d_halfpi_risbo_openmp(int lmax);

    const arr2<double> &recurse ();
  };

/*! Class for calculation of the Wigner matrix at arbitrary angles, using Risbo
    recursion in a way that can be OpenMP-parallelised. This approach uses more
    memory and is slightly slower than wigner_d_risbo_scalar. */
class wigner_d_risbo_openmp
  {
  private:
    double p,q;
    arr<double> sqt;
    arr2<double> d, dd;
    int n;

  public:
    wigner_d_risbo_openmp(int lmax, double ang);

    const arr2<double> &recurse ();
  };

/*! Class for calculation of the Wigner matrix elements by l-recursion.
    For details, see Prezeau & Reinecke 2010, http://arxiv.org/pdf/1002.1050 */
class wignergen_scalar
  {
  protected:
    typedef double dbl3[3];

    // members set in the constructor and kept fixed afterwards
    double fsmall, fbig, eps;
    int lmax;
    arr<long double> logsum, lc05, ls05;
    arr<double> flm1, flm2, cf, costh, xl;
    arr<bool> thetaflip;

    // members depending on m and m'
    int m1, m2, am1, am2, mlo, mhi, cosPow, sinPow;
    long double prefactor;
    arr<dbl3> fx;
    bool preMinus;

    // members depending on theta
    arr<double> result;

    enum { large_exponent2=90, minscale=-4, maxscale=14 };

  public:
    /*! Constructs an object that can compute Wigner matrix elements up
        to a maximum \a l value of \a lmax_, at the colatitudes provided
        in \a thetas. The generator will be allowed to regard values with
        absolute magnitudes smaller than \a epsilon as zero; a typical value
        is 1e-30. */
    wignergen_scalar (int lmax_, const arr<double> &thetas, double epsilon);

    /*! Prepares the object to produce Wigner matrix elements with \a m=m1_
        and \a m'=m2_ in subsequent calls to calc(). This operation is not cheap
        so it is recommended to use calc() for many different colatitudes after
        every call to prepare(), if possible. */
    void prepare (int m1_, int m2_);

    /*! Computes the Wigner matrix elements for the values of \a m and \a m'
        set by the preceding call to prepare(), for all \a l up to \a lmax
        (set in the constructor), and for the \a nth colatitude passed to the
        constructor. On return, \a firstl contains the index of the first
        matrix element larger than \a epsilon; all values with smaller indices
        in the result array are undefined. */
    const arr<double> &calc (int nth, int &firstl);
    void calc (int nth, int &firstl, arr<double> &resx) const;
  };

class wignergen: public wignergen_scalar
  {
#ifdef __SSE2__
  private:
    arr_align<V2df,16> result2;

  public:
    wignergen (int lmax_, const arr<double> &thetas, double epsilon)
      : wignergen_scalar (lmax_,thetas,epsilon), result2(lmax_+1) {}

    using wignergen_scalar::calc;
    const arr_align<V2df,16> &calc (int nth1, int nth2, int &firstl);
    void calc (int nth1, int nth2, int &firstl, arr_align<V2df,16> &resx) const;
#else
  public:
    wignergen (int lmax_, const arr<double> &thetas, double epsilon)
      : wignergen_scalar (lmax_,thetas,epsilon) {}
#endif
  };

class wigner_estimator
  {
  private:
    int lmax, m1, m2, mbig;
    double xlmax, epsPow, cosm1m2;

  public:
    wigner_estimator (int lmax_, double epsPow_);

    void prepare_m (int m1_, int m2_);
    bool canSkip (double theta) const;
  };

#endif

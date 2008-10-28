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
 *  Code for efficient calculation of Y_lm(theta,phi=0)
 *
 *  Copyright (C) 2005, 2006 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_YLMGEN_H
#define PLANCK_YLMGEN_H

#include <cmath>
#include "arr.h"
#include "lsconstants.h"

/*! Class for efficient calculation of Y_lm(theta,phi=0) */
class Ylmgen
  {
  private:
    double fsmall, fbig, eps, cth_crit;
    int lmax, mmax, m_last, m_crit;
    arr<double> cf;
    arr<double[2]> recfac;
    arr<double> mfac;
    arr<double> t1fac, t2fac;

    enum { large_exponent2 = 90, minscale=-4 };

    void recalc_recfac (int m)
      {
      using namespace std;

      if (m_last==m) return;

      double f_old=1;
      for (int l=m; l<recfac.size(); ++l)
        {
        recfac[l][0] = t1fac[l]*t2fac[l+m]*t2fac[l-m];
        recfac[l][1] = recfac[l][0]/f_old;
        f_old = recfac[l][0];
        }

      m_last=m;
      }

  public:
    /*! Creates a generator which will calculate Y_lm(theta,phi=0)
        up to \a l=l_max and \a m=m_max. It may regard Y_lm whose absolute
        magnitude is smaller than \a epsilon as zero. */
    Ylmgen (int l_max, int m_max, double epsilon=1e-30)
      : eps(epsilon), cth_crit(2.), lmax(l_max), mmax(m_max), m_last(-1),
        m_crit(mmax+1), cf(-minscale+11), recfac(lmax+1), mfac(mmax+1),
        t1fac(lmax+1), t2fac(lmax+mmax+1)
      {
      using namespace std;

      fsmall = ldexp(1.,-large_exponent2);
      fbig   = ldexp(1., large_exponent2);
      for (int m=0; m<cf.size(); ++m)
        cf[m] = ldexp(1.,(m+minscale)*large_exponent2);

      mfac[0] = 1;
      for (int m=1; m<mfac.size(); ++m)
        mfac[m] = mfac[m-1]*sqrt((2*m+1.)/(2*m));
      for (int m=0; m<mfac.size(); ++m)
        mfac[m] = inv_ln2*log(inv_sqrt4pi*mfac[m]);
      for (int l=0; l<t1fac.size(); ++l)
        t1fac[l] = sqrt(4.*(l+1)*(l+1)-1.);
      for (int i=0; i<t2fac.size(); ++i)
        t2fac[i] = 1./sqrt(i+1.);
      }

    /*! For a colatitude given by \a cth and \a sth (representing cos(theta)
        and sin(theta)) and a multipole moment \a m, calculate the
        Y_lm(theta,phi=0) for \a m<=l<=lmax and return in it \a result[l].
        On exit, \a firstl is the \a l index of the first Y_lm with an
        absolute magnitude larger than \a epsilon. If \a firstl>lmax, all
        absolute values are smaller than \a epsilon.
        \a result[l] is undefined for all \a l<firstl. */
    void get_Ylm (double cth, double sth, int m, arr<double> &result,
      int &firstl)
      {
      using namespace std;

      planck_assert (m<=mmax, "get_Ylm: m larger than mmax");

      if (((m>=m_crit)&&(abs(cth)>=cth_crit)) || ((m>0)&&(sth==0)))
        { firstl=lmax+1; return; }

      recalc_recfac(m);
      result.alloc(lmax+1);

      double logval = mfac[m];
      if (m>0) logval += m*inv_ln2*log(sth);
      int scale = int (logval/large_exponent2)-minscale;
      double corfac = (scale<0) ? 0. : cf[scale];
      double lam_1 = 0;
      double lam_2 = exp(ln2*(logval-(scale+minscale)*large_exponent2));
      if (m&1) lam_2 = -lam_2;

      int l=m;
      while (true)
        {
        if (abs(lam_2*corfac)>eps) break;
        if (++l>lmax) break;
        double lam_0 = cth*lam_2*recfac[l-1][0] - lam_1*recfac[l-1][1];
        if (abs(lam_0*corfac)>eps) { lam_1=lam_2; lam_2=lam_0; break; }
        if (++l>lmax) break;
        lam_1 = cth*lam_0*recfac[l-1][0] - lam_2*recfac[l-1][1];
        if (abs(lam_1*corfac)>eps) { lam_2=lam_1; lam_1=lam_0; break; }
        if (++l>lmax) break;
        lam_2 = cth*lam_1*recfac[l-1][0] - lam_0*recfac[l-1][1];

        while (abs(lam_2)>fbig)
          {
          lam_1 *= fsmall;
          lam_2 *= fsmall;
          ++scale;
          corfac = (scale<0) ? 0. : cf[scale];
          }
        }

      firstl=l;
      if (l>lmax)
        { m_crit=m; cth_crit=abs(cth); return; }

      lam_1*=corfac;
      lam_2*=corfac;

      for (;l<lmax-2;l+=3)
        {
        result[l]=lam_2;
        double lam_0 = cth*lam_2*recfac[l][0] - lam_1*recfac[l][1];
        result[l+1] = lam_0;
        lam_1 = cth*lam_0*recfac[l+1][0] - lam_2*recfac[l+1][1];
        result[l+2] = lam_1;
        lam_2 = cth*lam_1*recfac[l+2][0] - lam_0*recfac[l+2][1];
        }
      while (true)
        {
        result[l]=lam_2;
        if (++l>lmax) break;
        double lam_0 = cth*lam_2*recfac[l-1][0] - lam_1*recfac[l-1][1];
        result[l] = lam_0;
        if (++l>lmax) break;
        lam_1 = cth*lam_0*recfac[l-1][0] - lam_2*recfac[l-1][1];
        result[l] = lam_1;
        if (++l>lmax) break;
        lam_2 = cth*lam_1*recfac[l-1][0] - lam_0*recfac[l-1][1];
        }
      }

  };

#endif

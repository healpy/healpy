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

/*! \file alm_map_tools.h
 *  Copyright (C) 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALM_MAP_TOOLS_H
#define PLANCK_ALM_MAP_TOOLS_H

#include <vector>
#include "cxxutils.h"
#include "alm.h"
#include "xcomplex.h"

/*! A class holding information about a ring of pixels in a spherical map. */
class ringinfo
  {
  public:
    double theta, phi0, weight, cth, sth;
    int nph, ofs;

    ringinfo()
      : nph(0) {}
    /*! Constructs a \a ringinfo object.
        \param theta_ colatitude of the ring (in radian)
        \param phi0_ longitude of the first pixel in the ring (in radian)
        \param weight_ weighting factor for all pixels in the ring. This
               is typically the surface of a pixel in sterad.
        \note \a weight_ is only needed for map analysis, not synthesis.
        \param nph_ number of pixels in the ring
        \param ofs_ index of the first ring pixel in the total map array
               (counting from zero) */
    ringinfo (double theta_, double phi0_, double weight_, int nph_, int ofs_)
      : theta(theta_), phi0(phi0_), weight(weight_),
        cth(cos(theta)), sth(sin(theta)), nph(nph_), ofs(ofs_)
      {}
  };

/*! A class holding information about a ring pair in a spherical map. */
class ringpair
  {
  public:
    ringinfo r1, r2;

    /*! Initialize the object with the ring described by \a info.
        The second ring is left empty. */
    ringpair (const ringinfo &info)
      : r1(info) {}
    /*! Initialize the object with the rings described by \a info1
        and \a info2.
        \note The colatitude of \a info2 must be \f$\pi\f$ minus the colatitude
        of \a info1. */
    ringpair (const ringinfo &info1,const ringinfo &info2)
      : r1(info1), r2(info2)
      {
      planck_assert(approx(r1.theta,pi-r2.theta,1e-10), "invalid ringpair");
      }
  };

void info2pair(const std::vector<ringinfo> &info, std::vector<ringpair> &pair);

template<typename T> void map2alm (const std::vector<ringpair> &pair,
  const T *map, Alm<xcomplex<T> > &alm, bool add_alm);

template<typename T> void map2alm_pol
  (const std::vector<ringpair> &pair, const T *mapT, const T *mapQ,
   const T *mapU, Alm<xcomplex<T> > &almT, Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC, bool add_alm);

template<typename T> void map2alm_pol
  (const std::vector<ringpair> &pairT, const T *mapT,
   const std::vector<ringpair> &pairQ, const T *mapQ,
   const std::vector<ringpair> &pairU, const T *mapU,
   Alm<xcomplex<T> > &almT, Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC, bool add_alm);

template<typename T> void alm2map (const Alm<xcomplex<T> > &alm,
  const std::vector<ringpair> &pair, T *map);

template<typename T> void alm2map_pol
  (const Alm<xcomplex<T> > &almT, const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC, const std::vector<ringpair> &pair,
   T *mapT, T *mapQ, T *mapU);

template<typename T> void alm2map_der1 (const Alm<xcomplex<T> > &alm,
   const std::vector<ringpair> &pair, T *map, T *dth, T *dph);

#endif

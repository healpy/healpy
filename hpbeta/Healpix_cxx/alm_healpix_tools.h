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

/*! \file alm_healpix_tools.h
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_ALM_HEALPIX_TOOLS_H
#define HEALPIX_ALM_HEALPIX_TOOLS_H

#include "xcomplex.h"
#include "arr.h"

template<typename T> class Alm;
template<typename T> class Healpix_Map;

/*! \defgroup alm_healpix_group Conversions between a_lm and HEALPix maps */
/*! \{ */

/*! Converts a Healpix map to a set of a_lms.
    \param map the input map, which must have RING ordering
    \param alm the output a_lms. l_max and m_max of the conversion are
           determined from this object.
    \param weight array containing the weights for the individual rings of
           the map. It must have at least 2*\a map.Nside() entries.
    \param add_alm If this is \a true, then the computed a_lm are added
           to the values already residing in \a alm. */
template<typename T> void map2alm (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, const arr<double> &weight,
    bool add_alm=false);

/*! Converts a Healpix map to a set of a_lms, using an iterative scheme
    which is more accurate than plain map2alm().
    \param map the input map, which must have RING ordering.
    \param alm the output a_lms. l_max and m_max of the conversion are
           determined from this object.
    \param num_iter the number of iterations (0 is identical to map2alm()).
    \param weight array containing the weights for the individual rings of
           the map. It must have at least 2*\a map.Nside() entries. */
template<typename T> void map2alm_iter (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, int num_iter, const arr<double> &weight);

template<typename T> void map2alm_iter (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, int num_iter)
  {
  arr<double> wgt(2*map.Nside());
  wgt.fill(1);
  map2alm_iter(map,alm,num_iter,wgt);
  }

template<typename T> void map2alm_iter2 (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, double err_abs, double err_rel);

/*! Converts Healpix maps containing the I, Q and U Stokes parameters
    to sets of a_lms.
    \param mapT the I-Stokes parameter input map
    \param mapQ the Q-Stokes parameter input map
    \param mapU the U-Stokes parameter input map
    \note All maps must have the same nside, and must be in RING scheme.
    \param almT the output temperature a_lms
    \param almG the output gradient a_lms
    \param almC the output curl a_lms
    \note all a_lm sets must have the the same lmax and mmax.
    \param weight ring weights for the maps.
    \param add_alm If this is \a true, then the computed a_lm are added
           to the values already residing in \a alm.
    \note The weight array must have at least 2*\a mapT.Nside() entries. */
template<typename T> void map2alm_pol
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   const arr<double> &weight,
   bool add_alm=false);
/*! Converts Healpix maps containing the I, Q and U Stokes parameters
    to sets of a_lms, using an iterative scheme which is more accurate than
    plain map2alm_pol().
    \param mapT the I-Stokes parameter input map
    \param mapQ the Q-Stokes parameter input map
    \param mapU the U-Stokes parameter input map
    \note All maps must have the same nside, and must be in RING scheme.
    \param almT the output temperature a_lms
    \param almG the output gradient a_lms
    \param almC the output curl a_lms
    \note all a_lm sets must have the the same lmax and mmax.
    \param num_iter the number of iterations (0 is identical to map2alm_pol()).
    \param weight ring weights for the maps.
    \note The weight array must have at least 2*\a mapT.Nside() entries. */
template<typename T> void map2alm_pol_iter
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   int num_iter,
   const arr<double> &weight);

template<typename T> void map2alm_pol_iter
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   int num_iter)
  {
  arr<double> wgt(2*mapT.Nside());
  wgt.fill(1);
  map2alm_pol_iter(mapT,mapQ,mapU,almT,almG,almC,num_iter,wgt);
  }

template<typename T> void map2alm_pol_iter2
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   double err_abs, double err_rel);

/*! Converts a a set of a_lm to a HEALPix map.
    \param alm the input a_lms. l_max and m_max of the conversion are
           determined from this object.
    \param map the output map, which must have RING ordering. */
template<typename T> void alm2map (const Alm<xcomplex<T> > &alm,
  Healpix_Map<T> &map);

/*! Converts a a set of polarised a_lm to a HEALPix map.
    \param almT the input temperature a_lms
    \param almG the input gradient a_lms
    \param almC the input curl a_lms
    \param mapT the I-Stokes parameter output map
    \param mapQ the Q-Stokes parameter output map
    \param mapU the U-Stokes parameter output map */
template<typename T> void alm2map_pol
  (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC,
   Healpix_Map<T> &mapT,
   Healpix_Map<T> &mapQ,
   Healpix_Map<T> &mapU);

/*! Converts a a set of a_lm to a HEALPix map and its first derivatives.
    \param alm the input a_lms. l_max and m_max of the conversion are
           determined from this object.
    \param map the output map, which must have RING ordering.
    \param mapdth an output map containing \f$d (\mbox{map})/d\vartheta\f$,
           which must have RING ordering.
    \param mapdph an output map containing
           \f$(\sin\vartheta)^{-1}d(\mbox{map})/d\varphi\f$,
           which must have RING ordering. */
template<typename T> void alm2map_der1
  (const Alm<xcomplex<T> > &alm,
   Healpix_Map<T> &map,
   Healpix_Map<T> &mapdth,
   Healpix_Map<T> &mapdph);

/*! \} */

#endif

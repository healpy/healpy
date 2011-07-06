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

/*! \file trafos.h
 *  Celestial coordinate transformations.
 *
 *  Copyright (C) 2005-2011 Max-Planck-Society
 * \author Martin Reinecke
 */
#ifndef PLANCK_TRAFOS_H
#define PLANCK_TRAFOS_H

#include "vec3.h"
#include "pointing.h"
#include "rotmatrix.h"

enum coordsys { Ecliptic, Equatorial, Galactic };

/*! Class for celestial coordinate transformations. */
class Trafo
  {
  private:
    rotmatrix mat;

    static vec3 xcc_dp_precess (const vec3 &iv, double iepoch, double oepoch);
    static double get_epsilon (double epoch);
    static vec3 xcc_dp_e_to_q (const vec3 &iv, double epoch);
    static vec3 xcc_dp_q_to_e (const vec3 &iv, double epoch);
    static vec3 xcc_dp_g_to_e (const vec3 &iv, double epoch);
    static vec3 xcc_dp_e_to_g (const vec3 &iv, double epoch);
    static vec3 xcc_v_convert(const vec3 &iv, double iepoch, double oepoch,
      coordsys isys,coordsys osys);
    static void coordsys2matrix (double iepoch, double oepoch, coordsys isys,
      coordsys osys, rotmatrix &matrix);

  public:
    /*! Creates a \a Trafo for transformation from \a iepoch and \a isys
        to \a oepoch and \a osys. */
    Trafo (double iepoch, double oepoch, coordsys isys, coordsys osys);

    /*! Transforms the vector \a vec and returns the result. */
    vec3 operator() (const vec3 &vec) const
      { return mat.Transform(vec); }

    /*! Transforms the pointing \a ptg and returns the result. */
    pointing operator() (const pointing &ptg) const;

    /*! Transforms the pointing \a ptg and returns it in \a newptg.
        On exit, \a delta_psi holds the change in orientation. */
    void rotatefull (const pointing &ptg, pointing &newptg,
      double &delta_psi) const;

    /*! Transforms the pointing \a ptg and adjusts \a psi accordingly. */
    void rotatefull (pointing &ptg, double &psi) const;

    /*! Transforms the vector \a vec and returns it in \a newvec.
        On exit, \a delta_psi holds the change in orientation. */
    void rotatefull (const vec3 &vec, vec3 &newvec, double &delta_psi) const;

    /*! Transforms the vector \a vec and adjusts \a psi accordingly. */
    void rotatefull (vec3 &vec, double &psi) const;

    /*! Returns the internally used rotation matrix. */
    const rotmatrix &Matrix() const
      { return mat; }
  };

#endif

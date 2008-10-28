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

/*! \file trafos.cc
 *  Celestial coordinate transformations.
 *
 *  Copyright (C) 2005 Max-Planck-Society
 * \author Martin Reinecke
 */

#include "trafos.h"
#include "lsconstants.h"

using namespace std;

vec3 Trafo::xcc_dp_precess (const vec3 &iv, double iepoch,
  double oepoch)
  {
  // Z-axis rotation by OBL_LONG:
  double Tm = ((oepoch+iepoch)*0.5 - 1900.) *0.01;
  double gp_long  = (oepoch-iepoch) * (50.2564+0.0222*Tm) / 3600.;
  double obl_long = 180. - (173. + (57.06+54.77*Tm) / 60.) + gp_long*0.5;
  double dco = cos(obl_long*degr2rad), dso = sin(obl_long*degr2rad);
  vec3 ov (iv.x*dco-iv.y*dso, iv.x*dso+iv.y*dco, iv.z);

  // X-axis rotation by dE:
  double dE = (oepoch-iepoch) * (0.4711-0.0007*Tm) / 3600.;
  double dce = cos(dE*degr2rad), dse = sin(dE*degr2rad);
  double temp = ov.y*dce - ov.z*dse;
  ov.z = ov.y*dse + ov.z*dce;
  ov.y = temp;

  // Z-axis rotation by GP_LONG - OBL_LONG:
  double dL = gp_long - obl_long;
  double dcl = cos(dL*degr2rad), dsl = sin(dL*degr2rad);
  temp = ov.x*dcl - ov.y*dsl;
  ov.y = ov.x*dsl + ov.y*dcl;
  ov.x = temp;

  return ov;
  }

double Trafo::get_epsilon (double epoch)
  {
  double T = (epoch - 1900.) * 0.01; 
  double epsilon = 23.452294 - 0.0130125*T - 1.63889e-6*T*T + 5.02778e-7*T*T*T;
  return epsilon*degr2rad;
  }

/*! Routine to convert from ecliptic coordinates to equatorial (celestial)
    coordinates at the given epoch.  Adapted from the COBLIB routine by Dr.
    E. Wright. */
vec3 Trafo::xcc_dp_e_to_q (const vec3 &iv, double epoch)
  {
  double epsilon=get_epsilon(epoch);
  double dc = cos(epsilon), ds = sin(epsilon);
  return vec3 (iv.x, dc*iv.y-ds*iv.z, dc*iv.z+ds*iv.y);
  }

vec3 Trafo::xcc_dp_q_to_e (const vec3 &iv, double epoch)
  {
  double epsilon=-get_epsilon(epoch);
  double dc = cos(epsilon), ds = sin(epsilon);
  return vec3 (iv.x, dc*iv.y-ds*iv.z, dc*iv.z+ds*iv.y);
  }

/*! Routine to convert galactic coordinates to ecliptic (celestial) 
    coordinates at the given epoch.  First the conversion to ecliptic 
    2000 is done, then if necessary the results are precessed. */
vec3 Trafo::xcc_dp_g_to_e (const vec3 &iv, double epoch)
  {
  static rotmatrix T (-0.054882486,  0.494116468, -0.867661702,
                      -0.993821033, -0.110993846, -0.000346354,
                      -0.096476249,  0.862281440,  0.497154957);
  vec3 hv=T.Transform(iv);

  if (!approx(epoch,2000.))
    hv=xcc_dp_precess(hv,2000.,epoch);

  return hv;
  }

/*! Routine to convert ecliptic (celestial) coordinates at the given
    epoch to galactic coordinates.  The ecliptic coordinates are first
    precessed to 2000.0, then converted. */
vec3 Trafo::xcc_dp_e_to_g (const vec3 &iv, double epoch)
  {
  static rotmatrix T (-0.054882486, -0.993821033, -0.096476249,
                       0.494116468, -0.110993846,  0.862281440,
                      -0.867661702, -0.000346354,  0.497154957);
  vec3 hv=iv;
  if (!approx(epoch,2000.))
    hv=xcc_dp_precess(hv,epoch,2000.);

  return T.Transform(hv);
  }

/*! Function to convert between standard coordinate systems, including
    precession. */
vec3 Trafo::xcc_v_convert(const vec3 &iv, double iepoch, double oepoch,
  coordsys isys,coordsys osys)
  {
  vec3 xv;
  if (isys == Ecliptic)
    xv=iv;
  else if (isys == Equatorial)
    xv = xcc_dp_q_to_e(iv,iepoch);
  else if (isys == Galactic)
    xv = xcc_dp_g_to_e(iv,iepoch);
  else
    throw Message_error("Unknown input coordinate system");

  vec3 yv = approx(iepoch,oepoch) ? xv : xcc_dp_precess(xv,iepoch,oepoch);

  vec3 ov;
  if (osys == Ecliptic)
    ov = yv;
  else if (osys == Equatorial)
    ov = xcc_dp_e_to_q(yv,oepoch);
  else if (osys == Galactic)
    ov = xcc_dp_e_to_g(yv,oepoch);
  else
    throw Message_error("Unknown output coordinate system");

  return ov;
  }

void Trafo::coordsys2matrix (double iepoch, double oepoch,
  coordsys isys, coordsys osys, rotmatrix &matrix)
  {
  vec3 v1p = xcc_v_convert(vec3(1,0,0),iepoch,oepoch,isys,osys),
       v2p = xcc_v_convert(vec3(0,1,0),iepoch,oepoch,isys,osys),
       v3p = xcc_v_convert(vec3(0,0,1),iepoch,oepoch,isys,osys);
  v1p.Normalize(); v2p.Normalize(); v3p.Normalize();
  matrix=rotmatrix(v1p,v2p,v3p);
  }

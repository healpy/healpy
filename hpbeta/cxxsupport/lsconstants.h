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

/*! \file lsconstants.h
 *  Mathematical, physical and technical constants for LevelS.
 */

#ifndef PLANCK_CONSTANTS_H
#define PLANCK_CONSTANTS_H

/*! \defgroup mathconstgroup Mathematical constants */
/*! \{ */

const double pi=3.141592653589793238462643383279502884197;
const double twopi=6.283185307179586476925286766559005768394;
const double inv_twopi=1.0/twopi;
const double fourpi=12.56637061435917295385057353311801153679;
const double halfpi=1.570796326794896619231321691639751442099;
const double inv_halfpi=0.6366197723675813430755350534900574;
const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;

const double ln2 = 0.6931471805599453094172321214581766;
const double inv_ln2 = 1.4426950408889634073599246810018921;
const double ln10 = 2.3025850929940456840179914546843642;

const double onethird=1.0/3.0;
const double twothird=2.0/3.0;
const double fourthird=4.0/3.0;

const double degr2rad=pi/180.0;
const double arcmin2rad=degr2rad/60;
const double rad2degr=180.0/pi;

//! Ratio between FWHM and sigma of a Gauss curve (\f$\sqrt{8\ln2}\f$).
const double sigma2fwhm=2.3548200450309493; // sqrt(8*log(2.))
const double fwhm2sigma=1/sigma2fwhm;

/*! \} */

/*! \defgroup physconstgroup Physical constants */
/*! \{ */

const double Jansky2SI=1.0e-26;
const double SI2Jansky=1.0e+26;

//! Light speed in m/s (CODATA 2006)
const double speedOfLight=2.99792458e8;

//! Boltzmann's constant in J/K (CODATA 2006)
const double kBoltzmann=1.3806504e-23;

//! Stefan-Boltzmann constant in W/m^2/K^4 (CODATA 2006)
const double sigmaStefanBoltzmann=5.6704e-8;

//! Planck's constant in J s (CODATA 2006)
const double hPlanck=6.62606896e-34;

//! Astronomical unit in m
const double astronomicalUnit=1.49597870691e11;

//! Solar constant in W/m^2
const double solarConstant=1368.0;

//! Average CMB temperature in K (Mather et al. 1999, ApJ 512, 511)
const double tcmb = 2.725;

//! offset (in seconds) between Jan 1, 1958 and Jan 1, 1970
const double sec_58_70 = 378691200.;

/*! \} */

#endif

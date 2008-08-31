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

/*! \file alm_powspec_tools.h
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALM_POWSPEC_TOOLS_H
#define PLANCK_ALM_POWSPEC_TOOLS_H

#include "xcomplex.h"

template<typename T> class Alm;
class PowSpec;
class planck_rng;
class rotmatrix;

/*! \defgroup alm_ps_group Conversions between a_lm and power spectra */
/*! \{ */

/*! Creates a Gaussian realisation of the input power spectrum \a powspec,
    using the random number generator \a rng, and returns the result
    in \a alm. */
template<typename T> void create_alm (const PowSpec &powspec,
  Alm<xcomplex<T> > &alm, planck_rng &rng);

/*! Creates a Gaussian realisation of the polarised input power spectrum
    \a powspec, using the random number generator \a rng, and returns the
    result in \a almT, \a almG and \a almC. */
template<typename T> void create_alm_pol
  (const PowSpec &powspec,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   planck_rng &rng);

/*! Returns the unpolarised power spectrum of \a alm in \a powspec. */
template<typename T> void extract_powspec
  (const Alm<xcomplex<T> > &alm, PowSpec &powspec);
/*! Returns the cross power spectrum of \a alm1 and \a alm2 in \a powspec. */
template<typename T> void extract_crosspowspec
  (const Alm<xcomplex<T> > &alm1,
   const Alm<xcomplex<T> > &alm2, PowSpec &powspec);
/*! Returns the polarised power spectrum of \a almT, \a almG and \a almC
    in \a powspec. */
template<typename T> void extract_powspec
  (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC,
   PowSpec &powspec);

/*! \} */

/*! Applies a convolution with a Gaussian beam with an FWHM of
    \a fwhm_arcmin arcmin to \a alm.
    \note If \a fwhm_arcmin<0, a deconvolution with \a -fwhm_arcmin
      is performed.
    \relates Alm */
template<typename T> void smooth_with_Gauss
  (Alm<xcomplex<T> > &alm, double fwhm_arcmin);
/*! Applies a convolution with a Gaussian beam with an FWHM of
    \a fwhm_arcmin arcmin to \a almT, \a almG and \a almC.
    \note If \a fwhm_arcmin<0, a deconvolution with \a -fwhm_arcmin
      is performed.
    \relates Alm */
template<typename T> void smooth_with_Gauss
  (Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   double fwhm_arcmin);

/*! Rotates \a alm through the Euler angles \a psi, \a theta and \a phi.
    The Euler angle convention  is right handed, rotations are active.
    - \a psi is the first rotation about the z-axis (vertical)
    - then \a theta about the ORIGINAL (unrotated) y-axis
    - then \a phi  about the ORIGINAL (unrotated) z-axis (vertical)
    \relates Alm */
template<typename T> void rotate_alm (Alm<xcomplex<T> > &alm,
  double psi, double theta, double phi);

/*! Rotates \a almT, \a almG and \a almC through the Euler angles
    \a psi, \a theta and \a phi.
    The Euler angle convention  is right handed, rotations are active.
    - \a psi is the first rotation about the z-axis (vertical)
    - then \a theta about the ORIGINAL (unrotated) y-axis
    - then \a phi  about the ORIGINAL (unrotated) z-axis (vertical)
    \relates Alm */
template<typename T> void rotate_alm (Alm<xcomplex<T> > &almT,
  Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC,
  double psi, double theta, double phi);

/*! Rotates \a alm through the rotation matrix \a mat.
    \relates Alm */
template<typename T> void rotate_alm (Alm<xcomplex<T> > &alm,
  const rotmatrix &mat);

/*! Rotates \a almT, \a almG and \a almC through the rotation matrix \a mat.
    \relates Alm */
template<typename T> void rotate_alm (Alm<xcomplex<T> > &almT,
  Alm<xcomplex<T> > &almG, Alm<xcomplex<T> > &almC,
  const rotmatrix &mat);

#endif

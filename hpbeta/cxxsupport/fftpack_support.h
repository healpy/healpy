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
 *  Copyright (C) 2004 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_FFTPACK_SUPPORT_H
#define PLANCK_FFTPACK_SUPPORT_H

#include "ls_fft.h"
#include "arr.h"
#include "xcomplex.h"

class cfft
  {
  private:
    int n;
    complex_plan plan;

  public:
    cfft () : n(-1), plan(0) {}
    cfft (int size_)
      : n(size_), plan(make_complex_plan(size_)) {}
    ~cfft ()
      { if (plan!=0) kill_complex_plan (plan); }
    void Set (int size_)
      {
      if (plan!=0) kill_complex_plan (plan);
      n=size_;
      plan=make_complex_plan(size_);
      }

    int size() const
      { return n; }

    void forward (double *data)
      { complex_plan_forward(plan,data); }
    void backward (double *data)
      { complex_plan_backward(plan,data); }
    void forward (arr<xcomplex<double> >&data)
      { forward(&(data[0].re)); }
    void backward (arr<xcomplex<double> >&data)
      { backward(&(data[0].re)); }
  };

class rfft
  {
  private:
    int n;
    real_plan plan;

  public:
    rfft () : n(-1), plan(0) {}
    rfft (int size_)
      : n(size_), plan(make_real_plan(size_)) {}
    ~rfft ()
      { if (plan!=0) kill_real_plan (plan); }
    void Set (int size_)
      {
      if (plan!=0) kill_real_plan (plan);
      n=size_;
      plan=make_real_plan(size_);
      }

    int size() const
      { return n; }

    void forward_fftpack (double *data)
      { real_plan_forward_fftpack(plan,data); }
    void backward_fftpack (double *data)
      { real_plan_backward_fftpack(plan,data); }
    void forward_fftpack (arr<double> &data)
      { forward_fftpack(&(data[0])); }
    void backward_fftpack (arr<double> &data)
      { backward_fftpack(&(data[0])); }
    void forward_c (arr<xcomplex<double> >&data)
      { real_plan_forward_c(plan,&(data[0].re)); }
    void backward_c (arr<xcomplex<double> >&data)
      { real_plan_backward_c(plan,&(data[0].re)); }
  };

#endif

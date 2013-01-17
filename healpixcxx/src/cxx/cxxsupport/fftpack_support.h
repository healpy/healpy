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

/*
 *  Copyright (C) 2004-2011 Max-Planck-Society
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
    tsize n;
    complex_plan plan;

  public:
    cfft () : n(0), plan(0) {}
    cfft (tsize size_)
      : n(size_), plan(make_complex_plan(size_)) {}
    cfft (const cfft &orig)
      : n(orig.n), plan(copy_complex_plan(orig.plan)) {}
    ~cfft ()
      { if (plan!=0) kill_complex_plan (plan); }
    cfft &operator=(const cfft &orig)
      {
      if (n!=orig.n)
        {
        if (plan!=0) kill_complex_plan (plan);
        n=orig.n;
        plan = copy_complex_plan(orig.plan);
        }
      return *this;
      }
    void Set (tsize size_)
      {
      if (n==size_) return;
      if (plan!=0) kill_complex_plan (plan);
      n=size_;
      plan=make_complex_plan(size_);
      }

    tsize size() const
      { return n; }

    void forward (double *data)
      { complex_plan_forward(plan,data); }
    void backward (double *data)
      { complex_plan_backward(plan,data); }
    void forward (xcomplex<double> *data)
      { complex_plan_forward(plan,&(data->re)); }
    void backward (xcomplex<double> *data)
      { complex_plan_backward(plan,&(data->re)); }
    void forward (arr<xcomplex<double> >&data)
      { forward(&(data[0].re)); }
    void backward (arr<xcomplex<double> >&data)
      { backward(&(data[0].re)); }
  };

class rfft
  {
  private:
    tsize n;
    real_plan plan;

  public:
    rfft () : n(0), plan(0) {}
    rfft (const rfft &orig)
      : n(orig.n), plan(copy_real_plan(orig.plan)) {}
    rfft (tsize size_)
      : n(size_), plan(make_real_plan(size_)) {}
    ~rfft ()
      { if (plan!=0) kill_real_plan (plan); }
    rfft &operator=(const rfft &orig)
      {
      if (n!=orig.n)
        {
        if (plan!=0) kill_real_plan (plan);
        n=orig.n;
        plan = copy_real_plan(orig.plan);
        }
      return *this;
      }
    void Set (tsize size_)
      {
      if (n==size_) return;
      if (plan!=0) kill_real_plan (plan);
      n=size_;
      plan=make_real_plan(size_);
      }

    tsize size() const
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

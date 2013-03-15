#ifndef ALICE_UTILS_H
#define ALICE_UTILS_H

#include "vec3.h"
#include "announce.h"
#include "rotmatrix.h"
#include "PolarizationHolder.h"
#include "TextureHolder.h"
#include "lsconstants.h"

/*! Returns vectors north and east, given a normalized vector location
  on the unit sphere */
void get_north_east(const vec3 &location, vec3 &north, vec3 &east)
  {
  if (fabs(location.x) + fabs(location.y) > 0.0)
    east = vec3(-location.y,location.x,0).Norm();
  else
    east.Set(1.0,0,0);
  north = crossprod(location, east);
}

/*! Returns a normalized direction parallel to and orthogonal to the
  the polarization given by q and u at a location on the unit sphere.
  Healpix conventions for q and u are used.  */
void get_qu_direction(const vec3 &location, double q, double u,
                      vec3 &direction, vec3 &orthogonal)
  {
  vec3 north, east;
  get_north_east(location, north, east);
  double angle = safe_atan2(u, q) / 2.0;
  direction = (north * -cos(angle)) + (east * sin(angle));
  orthogonal = crossprod(location, direction);
  }

/*! Returns a new_location, an angle theta away from the old location,
   in the direction of the polarization given by q and u, and in the
   approximate direction of a right-handed rotation around the
   old_axis.  The axis around which one rotates to get to the
   new_location is also returned.  */
void get_step(const vec3 &location, double q, double u, double theta,
              const vec3 &old_axis, vec3 &new_location, vec3 &new_axis)
  {
  vec3 dummy;
  rotmatrix rot;

  get_qu_direction(location, q, u, dummy, new_axis);
  if (dotprod(new_axis, old_axis) < 0.0) new_axis.Flip();

  rot.Make_Axis_Rotation_Transform(new_axis, theta);
  rot.Transform(location, new_location);
  }

/*! Performs one Runge-Kutta second order step.  Values of Q and U
  must be correct as input, and they are updated at the end of this
  function.  */
void runge_kutta_step(const vec3 &old_location, const PolarizationHolder &ph,
                      double &q, double &u, double theta, const vec3 &old_axis,
                      vec3 &new_location, vec3 &new_axis, pointing &p)
  {
  // Take a half-theta step and get new values of Q and U.
  get_step(old_location, q, u, theta/2.0, old_axis, new_location, new_axis);
  p = pointing(new_location);
  ph.getQU(p, q, u);

  // Then take a full step, with those values of Q and U, from the
  // original location.
  get_step(old_location, q, u, theta, old_axis, new_location, new_axis);
  p = pointing(new_location);
  ph.getQU(p, q, u);
  }

/*! Second order Runge-Kutta integration on the sphere.  Given a
  starting location, a qu map of the sky, and a step size theta, this
  subroutine returns an array of pointings extending in both
  directions from the starting location.  */
void runge_kutta_2(const vec3 &location, const PolarizationHolder &ph,
                   double theta, arr< pointing > &pointings)
  {
  double q, u;
  int i = pointings.size();
  pointing p(location);
  vec3 first_axis, old_axis, new_axis, new_location, dummy, old_location;

  ph.getQU(p, q, u);
  get_qu_direction(location, q, u, dummy, first_axis);
  old_axis = first_axis;
  old_location = location;

  pointings[pointings.size() / 2] = p;

  for(i = 1 + pointings.size() / 2; i < int(pointings.size()); i++)
    {
    runge_kutta_step(old_location, ph, q, u, theta, old_axis, new_location, new_axis, p);
    old_axis = new_axis;
    old_location = new_location;
    pointings[i] = p;
    }

  old_axis = -first_axis;
  old_location = location;
  for(i = -1 + pointings.size() / 2; i >= 0; i--)
    {
    runge_kutta_step(old_location, ph, q, u, theta, old_axis, new_location, new_axis, p);
    old_axis = new_axis;
    old_location = new_location;
    pointings[i] = p;
    }
  }

/*! Get an array of texture values from an array of pointings */
void pointings_to_textures(arr< pointing > &curve, const TextureHolder &th,
                          arr< double > &textures)
  {
  textures.alloc(curve.size());
  for(tsize i = 0; i < curve.size(); i++)
    textures[i] = th.getTextureDouble(curve[i]);
  }

/*! Create a sinusoidal kernel. */
void make_kernel(arr< double > &kernel)
  {
  for(tsize i = 0; i < kernel.size(); i++)
    {
    double sinx = sin(pi * (i + 1.0) / (kernel.size() + 1.0));
    kernel[i] = sinx * sinx;
    }
  }

/*! Convolve an array with a kernel. */
void convolve(const arr< double > &kernel, const arr< double > &raw, arr< double > &convolution)
  {
  convolution.alloc(raw.size() - kernel.size() + 1);
  for(tsize i = 0; i < convolution.size(); i++)
    {
    double total = 0;
    for (tsize j = 0; j < kernel.size(); j++)
      total += kernel[j] * raw[i+j];
    convolution[i] = total;
    }
  }

#endif // ALICE_UTILS_H

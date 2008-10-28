#ifndef SKYMAP
#define SKYMAP

#include "pointing.h"
#include "rotmatrix.h"

// This holds the viewed version of the sky, whether a projection onto
// a plane, or a Healpix format fits file.
class SkyMap
{
 public:
  rotmatrix d_mat; /* for rotating the projection */
  int d_rotate; /* flag for doing rotation */
  virtual ~SkyMap() 
    {
      d_rotate = 0;
    };
  virtual double get_pixel(int i) const = 0;
  virtual void set_pixel(int i, double val) = 0;
  virtual void add_to_pixel(int i, double val) = 0;
  virtual int is_valid_pixel(int i) const = 0;
  virtual int max_pixel() const = 0;
  virtual int get_next_pixel(int i) const = 0;
  virtual int project(pointing p) const = 0;
  virtual pointing deproject(int i) const = 0;
  virtual void minmax(float &min, float &max) const = 0;

  pointing rotate(pointing p) const
    {
      if (d_rotate)
	{
	  vec3 v = d_mat.Transform(p.to_vec3());
	  pointing p_out(v);
	  return p_out;
	}
      else
	return p;
    };
  
  pointing unrotate(pointing p) const
    {
      if (d_rotate)
	{
	  rotmatrix t = d_mat;
	  t.Transpose();
	  vec3 v = t.Transform(p.to_vec3());
	  pointing p_out(v);
	  return p_out;
	}
      else
	return p;
    };
};

#endif // SKYMAP

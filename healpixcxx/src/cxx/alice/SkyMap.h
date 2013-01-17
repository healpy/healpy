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
    bool d_rotate; /* flag for doing rotation */
    virtual ~SkyMap() {}
    virtual double get_pixel(int i) const = 0;
    virtual void set_pixel(int i, double val) = 0;
    virtual void add_to_pixel(int i, double val) = 0;
    virtual bool is_valid_pixel(int i) const = 0;
    virtual int max_pixel() const = 0;
    virtual int get_next_pixel(int i) const = 0;
    virtual int project(pointing p) const = 0;
    virtual pointing deproject(int i) const = 0;
    virtual void minmax(float &min, float &max) const = 0;

    pointing rotate (const pointing &p) const
      { return d_rotate ? pointing(d_mat.Transform(p.to_vec3())) : p; }

    pointing unrotate(const pointing &p) const
      {
      if (!d_rotate) return p;
      rotmatrix t = d_mat;
      t.Transpose();
      return pointing(t.Transform(p.to_vec3()));
      }
  };

#endif // SKYMAP

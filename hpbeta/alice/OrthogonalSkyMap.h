#ifndef ORTHOGONALSKYMAP
#define ORTHOGONALSKYMAP

#include <assert.h>
#include <math.h>
#include <iostream>
#include "RectSkyMap.h"

using namespace std;

class OrthogonalSkyMap : public RectSkyMap
{
 private:
  double d_xmin, d_xscale;
  double d_ymin, d_yscale;
  
 public:
  OrthogonalSkyMap() 
    {
      set_size(100, -1.0, 1.0, -1.0, 1.0);
    };
  
  OrthogonalSkyMap(int width, double xmin, double xmax, double ymin, double ymax)
    {
      set_size(width, xmin, xmax, ymin, ymax);
    };
  
  OrthogonalSkyMap(int width)
    {
      set_size(width, -1.0, 1.0, -1.0, 1.0);
    };
  
  void set_size(int width, double xmin, double xmax, double ymin, double ymax)
    {
      assert(width > 0);
      assert(xmax > xmin);
      assert(ymax > ymin); 
      int height = static_cast<int>(floor(0.5 + width * (ymax - ymin) / (xmax - xmin)));
      this->RectSkyMap::set_size(width, height);
      d_xmin = xmin;
      d_ymin = ymin;
      d_xscale = (xmax - xmin) / width;
      d_yscale = (ymax - ymin) / height;
/*       cout << "d_xmin = " << d_xmin << endl; */
/*       cout << "d_ymin = " << d_ymin << endl; */
/*       cout << "d_xscale = " << d_xscale << endl; */
/*       cout << "d_yscale = " << d_yscale << endl; */
    };
  
  int is_valid_pixel(int i) const
    {
      int x, y;
      double xp, yp;
      if (i >= 0)
    {
      i2xy(i, x, y);
      xy2xpyp(x, y, xp, yp);
      return (xp * xp + yp * yp < 1.0);
    }
      else
    return 0;
    };

  void xy2xpyp(int x, int y, double &xp, double &yp) const
    {
      xp = (x + 0.5) * d_xscale + d_xmin;
      yp = (d_y - y - 0.5) * d_yscale + d_ymin;
    };
  
  void xpyp2xy(double xp, double yp, int &x, int &y) const
    {
      x = static_cast<int>(floor((xp - d_xmin) / d_xscale));
      y = static_cast<int>(floor(d_y - (yp - d_ymin) / d_yscale));      
    };
  
  // Return 1 if the pointing is on the visible half of the sphere.
  int ang2xpyp(pointing p, double &xp, double &yp) const
    {
      vec3 v;
      // p.normalize();
      v = p.to_vec3();
      xp = v.y;
      yp = v.x;
      
      return (v.z > 0.0);
    };
  
  // Return 1 if result is in the projected circle. 
  int xpyp2ang(double xp, double yp, pointing &p) const
    {
      double r2 = xp * xp + yp * yp;
      double z;
      if (r2 <= 1.0)
    {
      z = sqrt(1.0 - r2);
      vec2pnt(vec3(yp, xp, z), p);
      return 1;
    }
      else 
    {
      p.theta = 0.0;
      p.phi = 0.0;
      return 0;
    }
    };
  
  int project(pointing p) const
    {
      int i, x, y;
      double xp, yp;

      p = rotate(p);
      if (ang2xpyp(p, xp, yp))
    {
      xpyp2xy(xp, yp, x, y);
      if (x >= 0 && x < d_x && y >= 0 && y < d_y)
        {
          xy2i(x, y, i);
          return i;
        }
      else
        return -1;
    }
      else
    return -1;      
    };
  
  pointing deproject(int i) const
    {
      int x, y;
      double xp, yp;
      pointing p;
      i2xy(i, x, y);
      xy2xpyp(x, y, xp, yp);
      xpyp2ang(xp, yp, p);
      p = unrotate(p);
      return p;    
    };
};

#endif // ORTHOGONALSKYMAP

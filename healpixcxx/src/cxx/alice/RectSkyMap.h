#ifndef RECTSKYMAP
#define RECTSKYMAP

#include <cassert>
#include "SkyMap.h"
#include "arr.h"

// This holds a planar projection of the sky.  d_x is the width of the
// rectangle, and d_y is the height.
class RectSkyMap : public SkyMap
  {
  public:
    arr2<float> d_array;
    int d_x, d_y;

    RectSkyMap() : d_x(0), d_y (0) {}

    void i2xy(int i, int &x, int &y) const
      {
      x = i % d_x;
      y = i / d_x;
      }

    void xy2i(int x, int y, int &i) const
      { i = x + y * d_x; }

    void set_size(int x, int y)
      {
      planck_assert((x>0) && (y>0),"bad map sizes");
      d_x = x;
      d_y = y;
      d_array.alloc(x, y);
      d_array.fill(0.0);
      }

    double get_pixel(int i) const
      {
      int x, y;
      i2xy(i, x, y);
      return d_array[x][y];
      }

    void set_pixel(int i, double pix)
      {
      int x, y;
      i2xy(i, x, y);
      d_array[x][y] = pix;
      }

    void add_to_pixel(int i, double pix)
      {
      int x, y;
      i2xy(i, x, y);
      d_array[x][y] += pix;
      }

    int max_pixel() const
      { return d_x * d_y - 1; }

    int get_next_pixel(int i) const
      {
      i++;
      while(i < max_pixel() && !is_valid_pixel(i))
        i++;
      return i;
      }

    void minmax(float &min, float &max) const
      { d_array.minmax(min, max); }
  };

#endif // RECTSKYMAP

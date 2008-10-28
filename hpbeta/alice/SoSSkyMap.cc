#include "SoSSkyMap.h"
#include <math.h>

SoSSkyMap::SoSSkyMap(int x)
{
  set_size(x);
}

void SoSSkyMap::set_size(int x) 
{
  int y;
  assert(x >= 2);
  y = x / 2;
  this->RectSkyMap::set_size(x, y);
  d_xscale = 2 * pi / x;
  d_yscale = pi / y;
}
     
int SoSSkyMap::project(pointing p) const
{
  int x, y, i;
  x = static_cast<int>(floor(p.phi / d_xscale));
  y = static_cast<int>(floor(p.theta / d_yscale));
  xy2i(x, y, i);
  return i;
}

pointing SoSSkyMap::deproject(int i) const
{
  pointing p;
  int x, y;
  i2xy(i, x, y);
  p.phi = (x + 0.5) * d_xscale;
  p.theta = (y + 0.5) * d_yscale;
  return p;
}

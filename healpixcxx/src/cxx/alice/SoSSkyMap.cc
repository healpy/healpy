#include "SoSSkyMap.h"
#include "lsconstants.h"
#include <cmath>

using namespace std;

void SoSSkyMap::set_size(int x)
  {
  planck_assert(x>=2,"bad map size");
  int y = x / 2;
  RectSkyMap::set_size(x, y);
  d_xscale = 2 * pi / x;
  d_yscale = pi / y;
  }

int SoSSkyMap::project(pointing p) const
  {
  int i;
  int x = int(floor(p.phi / d_xscale));
  int y = int(floor(p.theta / d_yscale));
  xy2i(x, y, i);
  return i;
  }

pointing SoSSkyMap::deproject(int i) const
  {
  int x, y;
  i2xy(i, x, y);
  return pointing ((y+0.5) * d_yscale, (x+0.5) * d_xscale);
  }

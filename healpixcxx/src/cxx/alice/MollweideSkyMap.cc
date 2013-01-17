#include "MollweideSkyMap.h"
#include "lsconstants.h"
using namespace std;

void MollweideSkyMap::set_size(int x)
  {
  planck_assert(x>=2,"bad map size");
  RectSkyMap::set_size(x, x/2);
  d_halfx = double(d_x) / 2.0;
  d_halfy = double(d_y) / 2.0;
  roottwo = sqrt(2.0);
  }

bool MollweideSkyMap::is_valid_pixel(int i) const
  {
  int x, y;
  double xp, yp;
  i2xy(i, x, y);
  xy2xpyp(x, y, xp, yp);
  return (xp * xp / 4.0 + yp * yp < 2.0);
  }

void MollweideSkyMap::xy2xpyp(int x, int y, double &xp, double &yp) const
  {
  xp = 2 * roottwo * (x + 0.5 - d_halfx) / d_halfx;
  yp = -roottwo * (y + 0.5 - d_halfy) / d_halfy;
  }

void MollweideSkyMap::xpyp2xy(double xp, double yp, int &x, int &y) const
  {
  x = int(floor(d_x * (xp + 2 * roottwo) / (4 * roottwo)));
  y = int(floor(d_y * (-yp + roottwo) / (2 * roottwo)));
  }

int MollweideSkyMap::project(pointing p) const
  {
  p.normalize();
  double lat = pi / 2 - p.theta;
  double lon = p.phi;
  if (lon > pi) lon -= twopi;
  double pisinlat = pi * sin(lat);

  double diff = 1.0;
  double psi = lat;
  // Newton's method: numerically solving for psi.
  while(abs(diff) > 1e-13)
    {
    diff = 2 * psi + sin(2 * psi) - pisinlat;
    double slope = 2 * (1.0 + cos(2 * psi));
    if (slope == 0.0)
      slope = 1.0e-5;
    psi -= diff / slope;
    }

  double xp = -2 * roottwo * cos(psi) * lon / pi;
  double yp = roottwo * sin(psi);

  int i, x, y;
  xpyp2xy(xp, yp, x, y);
  xy2i(x, y, i);
  return i;
  }

pointing MollweideSkyMap::deproject(int i) const
  {
  int x, y;
  i2xy(i, x, y);
  double x2, y2;
  xy2xpyp(x, y, x2, y2);
  double x3 = sqrt(2.0 - y2 * y2);
  double psi = (x2==0.0) ? halfpi : atan(y2 / x3);
  double foo = (2 * psi + sin(2 * psi)) / pi;
  double lat = asin(foo);
  double lon = -x2 * pi / (2 * roottwo * fabs(cos(psi)));

  return pointing(halfpi-lat,lon);
  }

#include <iostream>
#include <stdio.h>
#include "MollweideSkyMap.h"
using namespace std;


MollweideSkyMap::MollweideSkyMap(int x)
{
  set_size(x);
}

void MollweideSkyMap::set_size(int x)
{ 
  int y;
  assert(x >= 2);
  y = x / 2;
  this->RectSkyMap::set_size(x, y);
  d_halfx = static_cast<double>(d_x) / 2.0;
  d_halfy = static_cast<double>(d_y) / 2.0;
  roottwo = sqrt(2.0);
}

int MollweideSkyMap::is_valid_pixel(int i) const
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
  x = static_cast<int>(floor(d_x * (xp + 2 * roottwo) / (4 * roottwo)));
  y = static_cast<int>(floor(d_y * (-yp + roottwo) / (2 * roottwo)));
}

int MollweideSkyMap::project(pointing p) const
{
  int i, x, y;
  double xp, yp, psi, lat, lon;
  double slope, diff;
  double pisinlat;
  
  p.normalize();
  lat = pi / 2 - p.theta;
  lon = p.phi;
  if (lon > pi) lon -= 2 * pi;
  // cout << "lat, lon project = " << lat << ' ' << lon << endl;
  pisinlat = pi * sin(lat);

  diff = 1.0;
  psi = lat;
  // Newton's method: numerically solving for psi.
  while(fabs(diff) > 1.0e-13)
    {
      diff = 2 * psi + sin(2 * psi) - pisinlat;
      slope = 2 * (1.0 + cos(2 * psi));
      if (slope == 0.0)
	slope = 1.0e-5;
      psi -= diff / slope;
      // printf("psi, diff = %20.15f %20.15f\n", psi, diff);
      // cout << "psi, diff = " << psi << ' ' << diff << endl;
    }

  xp = -2 * roottwo * cos(psi) * lon / pi;
  yp = roottwo * sin(psi);

  xpyp2xy(xp, yp, x, y);
  // printf("xp, yp = %f %f\n", xp, yp);
  // cout << "xp, yp = " << xp << " " << yp << endl;
  xy2i(x, y, i);
  return i;
}

pointing MollweideSkyMap::deproject(int i) const
{
  pointing p;
  int x, y;
  double x2, y2, psi, lat, lon, foo;
  double x3;

  i2xy(i, x, y);
  xy2xpyp(x, y, x2, y2);
  // cout << "i, x, y = " << i << " " << x << " " << y << endl;
  // cout << "xp, yp = " << x2 << " " << y2 << endl;
  x3 = sqrt(2.0 - y2 * y2);
  if (x2 == 0.0)
    psi = pi / 2.0;
  else
    psi = atan(y2 / x3);
  // cout << "psi = " << psi << endl;
  foo = (2 * psi + sin(2 * psi)) / pi;
  // cout << "foo = " << foo << endl;
  lat = asin(foo);
  lon = -x2 * pi / (2 * roottwo * fabs(cos(psi)));
  // cout << "lat, lon deproject = " << lat << ' ' << lon << endl;
  p.phi = lon;
  p.theta = pi / 2 - lat;

  return p;
}

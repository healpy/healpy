#include <cstdio>
#include <iostream>
#include <assert.h>
#include <math.h>
#include "MollweideSkyMap.h"


using namespace std;

int main()
{
  cout << "hello" << endl;

  printf("pi = %25.20f, %25.20f\n", pi * 1.0, (float) pi);

  MollweideSkyMap m;
  m.set_size(2048);

  int i, j, x, y, x2, y2;
  double xp, yp;

  for(i = 0; i <= m.max_pixel(); i++)
    {
      m.i2xy(i, x, y);
      m.xy2xpyp(x, y, xp, yp);
      m.xpyp2xy(xp, yp, x2, y2);
      // cout << i << ' ' << x << ' ' << y << ' ' << xp << ' ' << yp << ' ' << x2 << ' ' << y2 << endl;
      assert(x2 == x);
      assert(y2 == y);
      m.xy2i(x2, y2, j);
      assert(i == j);
    }
  cout << "test passed a" << endl;

  pointing p;
  p.theta = pi / 2 + 0.1;
  p.phi = 0.0;
  m.project(p);
  
  for(i = 0; i < 100; i++)
    {
      p.theta = i / 100.0 * pi;
      m.project(p);
      // cout << endl;
    }
  
  p.theta = pi;
  m.project(p);
  
  for(i = 0; i <= m.max_pixel(); i++)
    if (m.is_valid_pixel(i))
      {
	p = m.deproject(i);
	// cout << "pointing = " << p;
	j = m.project(p);
	// cout << i << ' ' << j << endl;
	assert(i == j);
	// return 0;
      }
  cout << "test passed c" << endl;

  return 0;
}

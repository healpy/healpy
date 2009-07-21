#include <cstdio>
#include <iostream>
#include <assert.h>
#include <math.h>
#include "OrthogonalSkyMap.h"

using namespace std;

int main()
{
  cout << "hello" << endl;

  printf("pi = %25.20f, %25.20f\n", pi * 1.0, (float) pi);

  OrthogonalSkyMap orth(1024);

  int i, j, x, y, x2, y2;
  double xp, yp;

  for(i = 0; i <= orth.max_pixel(); i++)
    {
      orth.i2xy(i, x, y);
      orth.xy2xpyp(x, y, xp, yp);
      orth.xpyp2xy(xp, yp, x2, y2);
      //cout << i << ' ' << x << ' ' << y << ' ' << xp << ' ' << yp << ' ' << x2 << ' ' << y2 << endl;
      assert(x2 == x);
      assert(y2 == y);
      orth.xy2i(x2, y2, j);
      assert(i == j);
    }
  cout << "test passed a" << endl;

  pointing p; 
  for(i = 0; i <= orth.max_pixel(); i++)
    if (orth.is_valid_pixel(i))
      {
	p = orth.deproject(i);
	// cout << "pointing = " << p;
	j = orth.project(p);
	// cout << i << ' ' << j << endl;
	assert(i == j);
	// return 0;
      }
  cout << "test passed b" << endl;

  return 0;
}

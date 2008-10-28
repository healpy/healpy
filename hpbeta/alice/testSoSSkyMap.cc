#include <assert.h>
#include <iostream>
using namespace std;

#include "SoSSkyMap.h"

int main()
{
  cout << "Hello world" << endl;

  SoSSkyMap s;

  s.set_size(4096);
  int x, y, i, j;
  for (i = 0; i <= s.max_pixel(); i++)
    {
      s.i2xy(i, x, y);
      s.xy2i(x, y, j);
      assert(i == j);
    }

  double d, d2;
  for (i = 0; i <= s.max_pixel(); i++)
    {
      d = static_cast<double>(i);
      s.set_pixel(i, d);
      d2 = s.get_pixel(i);
      assert(d == d2);
    }
  cout << "test passed a" << endl;

  pointing p;
  for (i = 0; i <= s.max_pixel(); i++)
    {
      p = s.deproject(i);
      j = s.project(p);
      assert(i == j);
    }
  cout << "test passed b" << endl;
  
  for (i = 0; i < s.max_pixel(); i++)
    {
      j = s.get_next_pixel(i);
      assert(j == i + 1);
    }
  j = s.get_next_pixel(s.max_pixel());
  assert(j == s.max_pixel() + 1);
  cout << "test passed cc" << endl;
  
  return 0;
}

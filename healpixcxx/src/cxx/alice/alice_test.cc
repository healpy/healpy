#include "alice_utils.h"
#include "lsconstants.h"

#include "RectSkyMap.h"
#include "MollweideSkyMap.h"
#include "OrthogonalSkyMap.h"
#include "SoSSkyMap.h"

using namespace std;

void testMollweide()
  {
  MollweideSkyMap m;
  m.set_size(2048);

  int i, j, x, y, x2, y2;
  double xp, yp;

  for(i = 0; i <= m.max_pixel(); i++)
    {
    m.i2xy(i, x, y);
    m.xy2xpyp(x, y, xp, yp);
    m.xpyp2xy(xp, yp, x2, y2);
    planck_assert(x2==x, "test a failed");
    planck_assert(y2==y, "test a failed");
    m.xy2i(x2, y2, j);
    planck_assert(i==j, "test a failed");
    }
  cout << "test passed a" << endl;

  for(i=0; i<=m.max_pixel(); i++)
    if (m.is_valid_pixel(i))
      {
      pointing p = m.deproject(i);
      j = m.project(p);
      planck_assert(i==j,"test b failed");
      }
  cout << "test passed b" << endl;
  }

void testOrthogonal()
  {
  OrthogonalSkyMap orth(1024);

  int i, j, x, y, x2, y2;
  double xp, yp;

  for(i=0; i<=orth.max_pixel(); i++)
    {
    orth.i2xy(i, x, y);
    orth.xy2xpyp(x, y, xp, yp);
    orth.xpyp2xy(xp, yp, x2, y2);
    planck_assert(x2==x,"test a failed");
    planck_assert(y2==y,"test a failed");
    orth.xy2i(x2, y2, j);
    planck_assert(i==j,"test a failed");
    }
  cout << "test passed a" << endl;

  for(i = 0; i <= orth.max_pixel(); i++)
    if (orth.is_valid_pixel(i))
      {
      pointing p = orth.deproject(i);
      j = orth.project(p);
      planck_assert(i==j,"test b failed");
      }
  cout << "test passed b" << endl;
  }

void testSoS()
  {
  SoSSkyMap s;
  s.set_size(4096);
  for (int i = 0; i <= s.max_pixel(); i++)
    {
    int x, y, j;
    s.i2xy(i, x, y);
    s.xy2i(x, y, j);
    assert(i == j);
    }

  for (int i = 0; i <= s.max_pixel(); i++)
    {
    double d = double(i);
    s.set_pixel(i, d);
    double d2 = s.get_pixel(i);
    assert(d == d2);
    }
  cout << "test passed a" << endl;

  for (int i = 0; i <= s.max_pixel(); i++)
    {
    pointing p = s.deproject(i);
    int j = s.project(p);
    assert(i == j);
    }
  cout << "test passed b" << endl;

  for (int i = 0; i < s.max_pixel(); i++)
    {
    int j = s.get_next_pixel(i);
    assert(j == i + 1);
    }
  cout << "test passed c" << endl;
  }

/*! Contains test code to do sanity checks as I go along. */
void testMisc()
{
  vec3 location, north, east, old_axis, new_axis, dummy, direction, orthogonal;
  vec3 new_location;
  double epsilon = 1.0e-15;

  cout << "Hello World!" << endl;

  // alice_utils.h
  location = vec3(1.0, 0.0, 0.0);
  get_north_east(location, north, east);
  assert((north - vec3(0.0, 0.0, 1.0)).Length() < epsilon);
  assert((east - vec3(0.0, 1.0, 0.0)).Length() < epsilon);

  location = vec3(0.0, 1.0, 0.0);
  get_north_east(location, north, east);
  assert((north - vec3(0.0, 0.0, 1.0)).Length() < epsilon);
  assert((east - vec3(-1.0, 0.0, 0.0)).Length() < epsilon);

  location = vec3(1.0, 0.0, 0.0);
  get_qu_direction(location, 1.0, 0.0, direction, orthogonal);
  assert((direction - vec3(0.0, 0.0, -1.0)).Length() < epsilon);
  assert((orthogonal - vec3(0.0, 1.0, 0.0)).Length() < epsilon);

  get_qu_direction(location, 0.0, 1.0, direction, orthogonal);
  dummy = vec3(0.0, 1.0, -1.0);
  dummy.Normalize();
  assert((direction - dummy).Length() < epsilon);
  dummy = vec3(0.0, 1.0, 1.0);
  dummy.Normalize();
  assert((orthogonal - dummy).Length() < epsilon);

  location = vec3(1.0, 0.0, 0.0);
  old_axis = vec3(0.0, -1.0, 0.1);
  get_step(location, 1.0, 0.0, halfpi, old_axis, new_location, new_axis);
  assert((new_location - vec3(0.0, 0.0, 1.0)).Length() < epsilon);
  assert((new_axis - vec3(0.0, -1.0, 0.0)).Length() < epsilon);

  get_step(location, 1.0, 0.0, halfpi / 2.0, old_axis, new_location, new_axis);
  dummy = vec3(1.0, 0.0, 1.0);
  dummy.Normalize();
  assert((new_location - dummy).Length() < epsilon);
  assert((new_axis - vec3(0.0, -1.0, 0.0)).Length() < epsilon);

  old_axis = vec3(0.0, 0.1, 1.0);
  get_step(location, -1.0, 0.0, halfpi / 2.0, old_axis, new_location, new_axis);
  dummy = vec3(1.0, 1.0, 0.0);
  dummy.Normalize();
  assert((new_location - dummy).Length() < epsilon);
  assert((new_axis - vec3(0.0, 0.0, 1.0)).Length() < epsilon);

  cout << "Tests passed 2" << endl;
  }

int main()
  {
  testMollweide();
  testOrthogonal();
  testSoS();
  testMisc();
  }

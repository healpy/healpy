#include "alice_utils.h"
#include "lsconstants.h"

#include "RectSkyMap.h"

using namespace std;

/*! Contains test code to do sanity checks as I go along. */
int main(void) 
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

  return 0;
}

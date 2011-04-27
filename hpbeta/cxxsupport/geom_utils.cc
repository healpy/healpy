#include "geom_utils.h"
#include "arr.h"

using namespace std;

namespace {

void get_circle (const arr<vec3> &point, tsize q1, tsize q2, vec3 &center,
  double &cosrad)
  {
  center = (point[q1]+point[q2]).Norm();
  cosrad = dotprod(point[q1],center);
  for (tsize i=0; i<q1; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      {
      center=crossprod(point[q1]-point[i],point[q2]-point[i]).Norm();
      cosrad=dotprod(point[i],center);
      if (cosrad<0)
        { center.Flip(); cosrad=-cosrad; }
      }
  }
void get_circle (const arr<vec3> &point, tsize q, vec3 &center,
  double &cosrad)
  {
  center = (point[0]+point[q]).Norm();
  cosrad = dotprod(point[0],center);
  for (tsize i=1; i<q; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      get_circle(point,i,q,center,cosrad);
  }

} // unnamed namespace

void find_enclosing_circle (const arr<vec3> &point, vec3 &center,
  double &cosrad)
  {
  tsize np=point.size();
  planck_assert(np>=3,"too few points");
  center = (point[0]+point[1]).Norm();
  cosrad = dotprod(point[0],center);
  for (tsize i=2; i<np; ++i)
    if (dotprod(point[i],center)<cosrad) // point outside the current circle
      get_circle(point,i,center,cosrad);
  }

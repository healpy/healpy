#include "PolarizationHolder.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"

using namespace std;

void PolarizationHolder::load(const string& filename)
  {
  fitshandle fh;
  fh.open(filename);
  fh.goto_hdu(2);
  read_Healpix_map_from_fits(fh, Q, 2);
  read_Healpix_map_from_fits(fh, U, 3);
  }

void PolarizationHolder::getQU(const pointing& p, float& q, float& u) const
  {
  q = Q.interpolated_value(p);
  u = U.interpolated_value(p);
  }

float PolarizationHolder::getQUMagnitude(const pointing& p) const
  {
  float q = Q.interpolated_value(p);
  float u = U.interpolated_value(p);
  return sqrt(q*q + u*u);
  }

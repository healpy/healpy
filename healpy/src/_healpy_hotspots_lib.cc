#include "_healpy_hotspots_lib.h"

void hotspots(const Healpix_Map<double> & inmap,
              Healpix_Map<double> & outmap,
              std::vector<int> & min_pixels,
              std::vector<int> & max_pixels)
{
  outmap.Set(inmap.Order(),inmap.Scheme());

  const int npix = inmap.Npix();
  min_pixels.reserve(npix);
  max_pixels.reserve(npix);
// #pragma omp parallel for schedule(dynamic,10000)
  for (int m=0; m<npix; ++m) {
    const double value = inmap[m];
    if (approx<double>(value, Healpix_undef)) continue;
    fix_arr<int,8> nb;
    inmap.neighbors(m,nb);
    bool is_max = true, is_min = true;
    for (tsize n=0; n<nb.size(); ++n) {
      if (nb[n] == -1) continue;
      const float nbval = inmap[nb[n]];
      if (approx<double>(nbval, Healpix_undef)) continue;
      if (nbval>=value) is_max = false;
      if (nbval<=value) is_min = false;
    }
    outmap[m] = (is_max || is_min) ? value : Healpix_undef;
    if (is_min) min_pixels.push_back(m);
    if (is_max) max_pixels.push_back(m);
  }

}

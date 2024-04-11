#ifndef HEALPIX_HOTSPOTS_H
#define HEALPIX_HOTSPOTS_H

#include <vector>
#include "healpix_map.h"

void hotspots(const Healpix_Map<double> & inmap,
              Healpix_Map<double> & outmap,
              std::vector<int> & min_pixels,
              std::vector<int> & max_pixels);

#endif

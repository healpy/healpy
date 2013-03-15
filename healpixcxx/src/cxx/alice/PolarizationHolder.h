#ifndef POLARIZATION_HOLDER
#define POLARIZATION_HOLDER

#include "healpix_map.h"

class PolarizationHolder
  {
  public:
    // Load a polarized fits file, with Q and U as the second
    // and third columns (the standard form).
    void load(const std::string& filename);

    // Return the polarization at some pointing.
    void getQU(const pointing& p, float& q, float& u) const;
    void getQU(const pointing& p, double& q, double& u) const
      {
      float qf, uf;
      getQU(p, qf, uf);
      q = double(qf);
      u = double(uf);
      }

    // Return the magnitude of the polarization at some pointing.
    float getQUMagnitude(const pointing& p) const;

    // Fiddling with these by hand shouldn't be necessary.
    Healpix_Map<float> Q, U;
  };

#endif // POLARIZATION_HOLDER

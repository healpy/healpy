#ifndef TEXTURE_HOLDER
#define TEXTURE_HOLDER

#include "healpix_map.h"
#include "planck_rng.h"
#include <string>

class TextureHolder
  {
  public:
    TextureHolder();

    void setToWhiteNoise(int Nside);
    void setToEllNoise(int Nside, int ell);
    void load(const std::string& filename);

    // Return the texture at some pointing.
    float getTexture(const pointing& p) const
      { return d_texture.interpolated_value(p); }

    double getTextureDouble(const pointing &p) const
      { return double(d_texture.interpolated_value(p)); }

    // Fiddling with this by hand shouldn't be necessary.
    // It ends up in NEST ordering, but that shouldn't matter.
    Healpix_Map<float> d_texture;

  private:
    planck_rng d_rng;
  };

#endif // TEXTURE_HOLDER

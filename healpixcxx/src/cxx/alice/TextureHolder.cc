#include "TextureHolder.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "lsconstants.h"
#include <iostream>

using namespace std;

TextureHolder::TextureHolder()
  { d_rng.seed(); }

void TextureHolder::setToWhiteNoise(int Nside)
  {
  d_texture.SetNside(Nside, NEST);
  for(int i = 0; i < d_texture.Npix(); i++)
    d_texture[i] = d_rng.rand_uni() - 0.5;
  }

void TextureHolder::setToEllNoise(int Nside, int ellIn)
  {
  d_texture.SetNside(Nside, RING);
  // Use Ben's value of ell, if ell = -1.
  int ellBig = (ellIn == -1) ?
    int (floor(5 * sqrt(3*pi) * Nside / 8 - 1)) : ellIn;

  Alm< xcomplex< float > > a(ellBig + 5, ellBig + 5);
  a.SetToZero();

  cout << "Background texture using ell = " << ellBig << endl;

  int ell = ellBig;
  for(int m = 0; m <= ell; m++)
    {
      a(ell, m).re = d_rng.rand_gauss();
      a(ell, m).im = d_rng.rand_gauss();
    }

  alm2map(a, d_texture);
  float min, max;
  d_texture.minmax(min, max);
  cout << "min, max = " << min << ", " << max << endl;

  // Leave d_texture in NEST ordering.
  d_texture.swap_scheme();

  cout << "Leaving setToEllNoise" << endl;
  }

void TextureHolder::load(const std::string &filename)
  {
  fitshandle fh;
  fh.open(filename);
  fh.goto_hdu(2);
  read_Healpix_map_from_fits(fh, d_texture, 1);
  }

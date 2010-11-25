#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "planck_rng.h"
#include "lsconstants.h"
#include <iostream>
#include <cstdio>
#include <cstring>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc != 4)
    {
      cerr << "Usage:\n"
        "  generateTexture <nside> <ell> <seed>\n\n"
        "  example: generateTexture 512 500 1 \n\n"
        "  All arguments are positive integers.  nside is a power of 2.\n"
        "  ell is the multipole moment ell.  seed is a random seed.\n";
      exit(1);
    }

  cout << "Generating texture" << endl;
  system("date");

  int nside = atoi(argv[1]);
  int ellTexture = atoi(argv[2]);
  int seed = atoi(argv[3]);

  cout << "Ben's default ell = " << static_cast<int>(floor(5 * sqrt(3*pi) * nside / 8 - 1)) << endl;
  cout << "using ell = " << ellTexture << ", nside = " << nside << endl;

  Healpix_Map< float > texture;
  texture.SetNside(nside, RING);

  Alm< xcomplex< float > > a;
  a.Set(ellTexture + 5, ellTexture + 5);
  int ell, m;
  for(ell = 0; ell <= a.Lmax(); ell++)
    for(m = 0; m <= ell; m++)
      {
        a(ell, m).re = 0.0;
        a(ell, m).im = 0.0;
      }

  planck_rng rng(seed);
  ell = ellTexture;
  for(m = 0; m <= ell; m++)
    {
      a(ell, m).re = rng.rand_gauss();
      a(ell, m).im = rng.rand_gauss();
    }

  if (nside > 128)
    cout << "The spherical harmonic transform may take a few minutes..." << endl;
  alm2map(a, texture);

  // cout << "Swapping scheme" << endl;
  cout << "Switching to NEST ordering" << endl;
  system("date");
  texture.swap_scheme();

  cout << "Writing file" << endl;
  system("date");
  char filename[1000];
  char buffer[1000];
  sprintf(filename, "texture_nside%5.5d_ell%5.5d_seed%6.6d.fits", nside, ellTexture, seed);
  sprintf(buffer, "rm %s", filename);
  system(buffer);

  cout << "Writing out file: " << filename << endl;
  fitshandle fh;
  fh.create(filename);
  arr< string > colname(1);
  colname[0] = "texture";
  prepare_Healpix_fitsmap(fh, texture, PLANCK_FLOAT32, colname);
  fh.write_column(1, texture.Map());
  fh.close();

  system("date");

  return 0;
}

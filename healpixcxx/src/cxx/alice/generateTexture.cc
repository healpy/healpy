#include <iostream>
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "planck_rng.h"
#include "lsconstants.h"
#include "string_utils.h"

using namespace std;

int main(int argc, char* argv[])
  {
  planck_assert(argc==4,
        "Usage:\n"
        "  generateTexture <nside> <ell> <seed>\n\n"
        "  example: generateTexture 512 500 1 \n\n"
        "  All arguments are positive integers.  nside is a power of 2.\n"
        "  ell is the multipole moment ell.  seed is a random seed.\n");

  cout << "Generating texture" << endl;

  int nside = stringToData<int>(argv[1]);
  int ellTexture = stringToData<int>(argv[2]);
  int seed = stringToData<int>(argv[3]);

  cout << "Ben's default ell = " << int(floor(5 * sqrt(3*pi) * nside / 8 - 1)) << endl;
  cout << "using ell = " << ellTexture << ", nside = " << nside << endl;

  Healpix_Map< float > texture(nside,RING,SET_NSIDE);

  Alm< xcomplex< float > > a(ellTexture + 5, ellTexture + 5);
  a.SetToZero();
  int ell, m;

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
  texture.swap_scheme();

  cout << "Writing file" << endl;
  string filename = "!texture_nside"+intToString(nside,5)+"_ell"
    +intToString(ellTexture,5)+"_seed"+intToString(seed,6)+".fits";

  cout << "Writing out file: " << filename << endl;
  fitshandle fh;
  fh.create(filename);
  arr<string> colname(1);
  colname[0] = "texture";
  prepare_Healpix_fitsmap(fh, texture, PLANCK_FLOAT32, colname);
  fh.write_column(1, texture.Map());
  }

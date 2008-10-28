#include "TextureHolder.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include <iostream>

using namespace std;

//----------------------------------------------------------------------
TextureHolder::TextureHolder()
{
  d_rng.seed();
}
//----------------------------------------------------------------------
void TextureHolder::setToWhiteNoise(int Nside)
{
  d_texture.SetNside(Nside, NEST);
  for(int i = 0; i < d_texture.Npix(); i++)
    d_texture[i] = d_rng.rand_uni() - 0.5;
// d_texture[i] = d_rng.rand_gauss();
//   float min, max;
//   d_texture.minmax(min, max);
//   cout << "texture min, max = " << min << " " << max << endl;
}
//----------------------------------------------------------------------
void TextureHolder::setToEllNoise(int Nside, int ellIn)
{
  int ellBig;
  int ell;
  
  d_texture.SetNside(Nside, RING);
  // Use Ben's value of ell, if ell = -1.
  if (ellIn == -1)
    ellBig = static_cast<int>(floor(5 * sqrt(3*pi) * Nside / 8 - 1));
  else
    ellBig = ellIn;
  
  Alm< xcomplex< float > > a;		
  a.Set(ellBig + 5, ellBig + 5);
  for(ell = 0; ell <= a.Lmax(); ell++)
    for(int m = 0; m <= ell; m++)
      {
	a(ell, m).re = 0.0;
	a(ell, m).im = 0.0;
      }
  
  cout << "Background texture using ell = " << ellBig << endl;
  
  ell = ellBig;
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
  
  /*
    d_texture.minmax(min, max);
    cout << "min, max = " << min << ", " << max << endl;
    
    fitshandle fh;
    system("rm test2Texture.fits");
    fh.create("test2Texture.fits");
    arr< string > colname(1);
    colname[0] = "texture";
    prepare_Healpix_fitsmap(fh, d_texture, TFLOAT, colname);
    fh.write_column(1, d_texture.Map());
    fh.close();
  */
  
  cout << "Leaving setToEllNoise" << endl;
}
//----------------------------------------------------------------------
void TextureHolder::load(const std::string& filename)
{
  fitshandle fh;
  fh.open(filename, READONLY);
  fh.goto_hdu(2);
  read_Healpix_map_from_fits(fh, d_texture, 1);
  fh.close();
}
//----------------------------------------------------------------------

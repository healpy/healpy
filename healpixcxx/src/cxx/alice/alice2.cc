// ALICE = Amazing Line Integral Convolution Executable
// Programmer: David Larson
// Date: January 29, 2006
// (Code originally written last fall.)
//
// Using default emacs indentation.


#include <iostream>
#include <cassert>
#include "paramfile.h"
#include "ls_image.h"
#include "healpix_map_fitsio.h"
#include "lsconstants.h"
#include "arr.h"
#include "fitshandle.h"
#include "PolarizationHolder.h"
#include "TextureHolder.h"
#include "SoSSkyMap.h"
#include "MollweideSkyMap.h"
#include "OrthogonalSkyMap.h"
#include "color.h"
#include "alice_utils.h"
#include "vec3.h"
#include "string_utils.h"

using namespace std;

typedef struct setup_struct
{
  string temperature_file;
  string polarization_file;
  string output_base;
  string texture_file;
  int texture_nside;
  int texture_ell;
  int x_size;
  bool sin_kernel;
  int steps;
  int kernel_steps;
  double step_radian;
  bool threshold;
  bool orthogonal;
  float xmin, xmax, ymin, ymax;
  float alpha, beta, gamma;
  int color_scale;
  bool bar;
  string title;
  bool only_fits_output;
  string fits_output_file;
  bool make_movie;
  string texture_file2;
  int movie_frames;
  bool sos;
  float t_min;
  float t_max;
  float pol_min;
  float pol_max;
} SetupVars;

// ----------------------------------------------------------------------
void usage()
{
#include "alice_usage.h"
  planck_fail_quietly("");
}

// ----------------------------------------------------------------------
void parse_argv(int argc, char* argv[], SetupVars  &s)
{
  s.output_base = string("test");
  s.texture_nside = -1;
  s.texture_ell = -1;
  s.x_size = 500;
  s.sin_kernel = true;
  s.steps = 100;
  s.kernel_steps = 50;
  s.step_radian = -1.0;
  s.threshold = false;
  s.orthogonal = false;
  s.xmin = -1.0;
  s.xmax = 1.0;
  s.ymin = -1.0;
  s.ymax = 1.0;
  s.alpha = 0.0;
  s.beta = 0.0;
  s.gamma = 0.0;
  s.color_scale = 0;
  s.bar = false;
  s.only_fits_output = false;
  s.make_movie = false;
  s.movie_frames = 3;
  s.sos = false;
  s.t_min = 123456.0;
  s.t_max = 123456.0;
  s.pol_min = 123456.0;
  s.pol_max = 123456.0;

  int m = 1;
  while (m < argc)
    {
      int mstart = m;
      string arg = argv[m];
      if (arg == "-in")
        {
          s.temperature_file = argv[m+1];
          s.polarization_file = argv[m+1];
          m += 2;
        }
      if (arg == "-temperature") { s.temperature_file=argv[m+1]; m+=2; }
      if (arg == "-out") { s.output_base=argv[m+1]; m+=2; }
      if (arg == "-texture") { s.texture_file=argv[m+1]; m+=2; }
      if (arg == "-nside") { stringToData(argv[m+1],s.texture_nside); m+=2; }
      if (arg == "-ell") { stringToData(argv[m+1],s.texture_ell); m+=2; }
      if (arg == "-xsz") { stringToData(argv[m+1],s.x_size); m+=2; }
      if (arg == "-flat") { s.sin_kernel=false; m+=1; }
      if (arg == "-steps") { stringToData(argv[m+1],s.steps); m+=2; }
      if (arg == "-kernel_steps") { stringToData(argv[m+1],s.kernel_steps); m+=2; }
      if (arg == "-step_arcmin") { stringToData(argv[m+1],s.step_radian); s.step_radian*=pi/180/60.0; m+=2; }
      if (arg == "-step_radian") { stringToData(argv[m+1],s.step_radian); m+=2; }
      if (arg == "-threshold") { s.threshold=true; m+=1; }
      if (arg == "-orth")
        {
          stringToData(argv[m+1],s.xmax);
          s.xmin = -s.xmax;
          s.ymax = s.xmax;
          s.ymin = s.xmin;
          s.orthogonal = true;
          m += 2;
        }
      if (arg == "-alpha") { stringToData(argv[m+1],s.alpha); m+=2; }
      if (arg == "-beta") { stringToData(argv[m+1],s.beta); m+=2; }
      if (arg == "-gamma") { stringToData(argv[m+1],s.gamma); m+=2; }
      if (arg == "-col") { stringToData(argv[m+1],s.color_scale); m+=2; }
      if (arg == "-bar") { s.bar=true; m+=1; }
      if (arg == "-title") { stringToData(argv[m+1],s.title); m+=2; }
      if (arg == "-fitsout") { s.fits_output_file=argv[m+1]; s.only_fits_output=true; m+=2; }
      if (arg == "-texture2") { s.texture_file2=argv[m+1]; s.make_movie=true; m+=2; }
      if (arg == "-frames") { stringToData(argv[m+1],s.movie_frames); m+=2; }
      if (arg == "-sos") { s.sos=true; m+=1; }
      if (arg == "-min") { stringToData(argv[m+1],s.t_min); m+=2; }
      if (arg == "-max") { stringToData(argv[m+1],s.t_max); m+=2; }
      if (arg == "-polmin") { stringToData(argv[m+1],s.pol_min); m+=2; }
      if (arg == "-polmax") { stringToData(argv[m+1],s.pol_max); m+=2; }

      if (mstart == m)
        {
          cout << "unrecognized option: " + arg << endl;
          usage();
        }
    }

  if (s.kernel_steps > s.steps)
    {
      cout << "Error: kernel_steps is larger than steps" << endl;
      usage();
    }

  if (s.step_radian < 0.0)
    s.step_radian = pi / s.x_size;

  s.alpha *= pi / 180.0;
  s.beta *= pi / 180.0;
  s.gamma *= pi / 180.0;
}

// ----------------------------------------------------------------------
void print_setup(const SetupVars &s)
{
  cout << "Alice: temperature_file = " << s.temperature_file << endl;
  cout << "Alice: polarization_file = " << s.polarization_file << endl;
  cout << "Alice: output_base = " << s.output_base << endl;
  cout << "Alice: texture_file = " << s.texture_file << endl;
  cout << "Alice: texture_nside = " << s.texture_nside << endl;
  cout << "Alice: texture_ell = " << s.texture_ell << endl;
  cout << "Alice: x_size = " << s.x_size << endl;
  cout << "Alice: sin_kernel = " << s.sin_kernel << endl;
  cout << "Alice: steps = " << s.steps << endl;
  cout << "Alice: kernel_steps = " << s.kernel_steps << endl;
  cout << "Alice: step_radian = " << s.step_radian << " (" <<
    s.step_radian * 180 * 60 / pi << " arcmin)" << endl;
  cout << "Alice: threshold = " << s.threshold << endl;
  cout << "Alice: orthogonal = " << s.orthogonal << endl;
  cout << "Alice: xmin = " << s.xmin << endl;
  cout << "Alice: xmax = " << s.xmax << endl;
  cout << "Alice: ymin = " << s.ymin << endl;
  cout << "Alice: ymax = " << s.ymax << endl;
  cout << "Alice: alpha = " << s.alpha << endl;
  cout << "Alice: beta = " << s.beta << endl;
  cout << "Alice: gamma = " << s.gamma << endl;
  cout << "Alice: color_scale = " << s.color_scale << endl;
  cout << "Alice: bar = " << s.bar << endl;
  cout << "Alice: title = " << s.title << endl;
  cout << "Alice: only_fits_output (-fitsout) = " << s.only_fits_output << endl;
  cout << "Alice: fits_output_file = " << s.fits_output_file << endl;
  cout << "Alice: make_movie = " << s.make_movie << endl;
  cout << "Alice: texture_file2 = " << s.texture_file2 << endl;
  cout << "Alice: movie_frames = " << s.movie_frames << endl;
  cout << "Alice: sos = " << s.sos << endl;
  cout << "Alice: t_min = " << s.t_min << endl;
  cout << "Alice: t_max = " << s.t_max << endl;
  cout << "Alice: pol_min = " << s.pol_min << endl;
  cout << "Alice: pol_max = " << s.pol_max << endl;
}

// ----------------------------------------------------------------------
void colorscale1(double& red, double& green, double& blue, double x)
{
  double hue, saturation, value;
  assert(x >= 0);
  assert(x <= 1.0);

  hue = 240.0;
  if (x < 0.5)
    {
      saturation = 1.0;
      value = 0.2 + 1.6 * x;
    }
  else
    {
      saturation = 1.0 - 2 * (x - 0.5);
      value = 1.0;
    }
  hsvToRgb(red, green, blue, hue, saturation, value);
}

// ----------------------------------------------------------------------
int lic_function(SkyMap &hitcount, SkyMap &texture, const PolarizationHolder &ph,
         const TextureHolder &th, int steps, int kernel_steps, double step_radian)
{
  arr< double > kernel, convolution, rawtexture;
  kernel.alloc(kernel_steps);
  make_kernel(kernel);

  arr< pointing > curve;
  curve.alloc(steps);

  int num_curves = 0;
  int k;
  pointing p;

  for(int i = 0; i <= texture.max_pixel(); i++)
    {
      if (texture.is_valid_pixel(i))
        {
          p = texture.deproject(i);
          if(hitcount.get_pixel(i) < 1.0)
            {
              num_curves++;
              runge_kutta_2(p.to_vec3(), ph, step_radian, curve);
              pointings_to_textures(curve, th, rawtexture);
              convolve(kernel, rawtexture, convolution);
              for(tsize j = 0; j < convolution.size(); j++)
                {
                  p = curve[j + kernel.size()/2];
                  k = texture.project(p);
                  if (texture.is_valid_pixel(k))
                    {
                      texture.add_to_pixel(k, convolution[j]);
                      hitcount.add_to_pixel(k, 1.0);
                    }
                }
            }
        }
    }

  return num_curves;
}

// ----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  SetupVars s;

  announce("alice 0.2.1");
  if (argc < 2) usage();
  parse_argv(argc, argv, s);
  print_setup(s);

  Colour myColour;
  Palette palette4;
  palette4.setPredefined(4);

  // Temperature
  Healpix_Map< float > temperature_fits;
  read_Healpix_map_from_fits(s.temperature_file, temperature_fits, 1, 2);

  // Polariztion
  PolarizationHolder ph;
  ph.load(s.polarization_file);

  // Texture
  TextureHolder th;
  if (s.texture_file == string(""))
    if (s.texture_nside > 0 && s.texture_ell > 0)
      th.setToEllNoise(s.texture_nside, s.texture_ell);
    else if (s.texture_nside > 0)
      th.setToWhiteNoise(s.texture_nside);
    else
      {
        cout << "Unable to find or create a texture file." << endl;
        usage();
      }
  else
    th.load(s.texture_file);

  int i;
  pointing p;
  int num_curves = 0;
  SkyMap *hitcount=0, *texture=0, *magnitude=0, *temperature=0;
  RectSkyMap *modtexture=0;

  if (s.sos) // Science on a sphere projection
    {
      hitcount = new SoSSkyMap(s.x_size);
      texture = new SoSSkyMap(s.x_size);
      magnitude = new SoSSkyMap(s.x_size);
      temperature = new SoSSkyMap(s.x_size);
      modtexture = new SoSSkyMap(s.x_size);
    }
  else if (s.orthogonal) // Orthogonal projection
    {
      hitcount = new OrthogonalSkyMap(s.x_size, s.xmin, s.xmax, s.ymin, s.ymax);
      texture = new OrthogonalSkyMap(s.x_size, s.xmin, s.xmax, s.ymin, s.ymax);
      magnitude = new OrthogonalSkyMap(s.x_size, s.xmin, s.xmax, s.ymin, s.ymax);
      temperature = new OrthogonalSkyMap(s.x_size, s.xmin, s.xmax, s.ymin, s.ymax);
      modtexture = new OrthogonalSkyMap(s.x_size, s.xmin, s.xmax, s.ymin, s.ymax);
    }
  else // Mollweide by default
    {
      hitcount = new MollweideSkyMap(s.x_size);
      texture =  new MollweideSkyMap(s.x_size);
      magnitude = new MollweideSkyMap(s.x_size);
      temperature = new MollweideSkyMap(s.x_size);
      modtexture = new MollweideSkyMap(s.x_size);
    }

  // Set up the necessary rotations.
  hitcount->d_mat.Make_CPAC_Euler_Matrix(s.alpha, s.beta, s.gamma);
  texture->d_mat = hitcount->d_mat;
  magnitude->d_mat = hitcount->d_mat;
  temperature->d_mat = hitcount->d_mat;
  modtexture->d_mat = hitcount->d_mat;

  if (s.alpha == 0.0 && s.beta == 0.0 && s.gamma == 0.0)
    {
      hitcount->d_rotate = 0;
      texture->d_rotate = 0;
      magnitude->d_rotate = 0;
      temperature->d_rotate = 0;
      modtexture->d_rotate = 0;
    }
  else
    {
      hitcount->d_rotate = 1;
      texture->d_rotate = 1;
      magnitude->d_rotate = 1;
      temperature->d_rotate = 1;
      modtexture->d_rotate = 1;
    }

  for (i = 0; i < magnitude->max_pixel(); i++)
    if (magnitude->is_valid_pixel(i))
      {
        p = magnitude->deproject(i);
        p.normalize();

        float temp = temperature_fits[temperature_fits.ang2pix(p)];
        if (s.t_max != 123456.0) temp = (temp > s.t_max) ? s.t_max : temp;
        if (s.t_min != 123456.0) temp = (temp < s.t_min) ? s.t_min : temp;
        temperature->set_pixel(i, temp);

        temp = ph.getQUMagnitude(p);
        if (s.pol_max != 123456.0) temp = (temp > s.pol_max) ? s.pol_max : temp;
        if (s.pol_min != 123456.0) temp = (temp < s.pol_min) ? s.pol_min : temp;
        magnitude->set_pixel(i, temp);

        texture->set_pixel(i, th.getTexture(p));
      }

  // Print out the background texture before we start the convolution.
  if (!s.only_fits_output)
    {
      int x, y;
      RectSkyMap *rtexture = static_cast<RectSkyMap*>(texture);
      LS_Image image(rtexture->d_array.size1(), rtexture->d_array.size2());
      image.fill(Colour(1.0, 1.0, 1.0));
      float min, max, foo;
      rtexture->minmax(min, max);
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rtexture->i2xy(i, x, y);
            foo = (rtexture->d_array[x][y] - min) / (max - min);
            image.put_pixel(x, y, Colour(foo, foo, foo));
          }
      string filename = s.output_base + "_background.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);
    }

  // Reset the texture to zero.
  for (i = 0; i < texture->max_pixel(); i++)
    {
      texture->set_pixel(i, 0.0);
    }

  // This does all the work.
  num_curves = lic_function(*hitcount, *texture, ph, th,
                            s.steps, s.kernel_steps, s.step_radian);

  cout << "foo3" << endl;
  int num_pix = 0;
  float hitmin, hitmax;
  float foo;

  hitcount->minmax(hitmin, hitmax);
  for (i = 0; i <= hitcount->max_pixel(); i++)
    if (hitcount->is_valid_pixel(i)) num_pix++;
  cout << endl;
  cout << "number of curves calculated = " << num_curves << endl;
  foo = static_cast<float>(num_pix) / num_curves;
  cout << "number of pixels in image = " << num_pix << endl;
  cout << "average pixels per curve = " << foo << endl;
  cout << "average curve points per good pixel = " << s.steps / foo << endl;
  cout << "minimum number of hits = " << hitmin << endl;
  cout << "maximum number of hits = " << hitmax << endl;

  if (!s.only_fits_output)
    {
      RectSkyMap *rhitcount = static_cast<RectSkyMap*>(hitcount);
      RectSkyMap *rtexture = static_cast<RectSkyMap*>(texture);
      RectSkyMap *rmagnitude = static_cast<RectSkyMap*>(magnitude);
      RectSkyMap *rtemperature = static_cast<RectSkyMap*>(temperature);
      int x, y;
      float min, max;
      float mmin, mmax;
      float tmin, tmax;

      LS_Image image(rtexture->d_array.size1(), rtexture->d_array.size2());
      image.fill(Colour(1.0, 1.0, 1.0));

      // Hitpattern
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rhitcount->i2xy(i, x, y);
            foo = rhitcount->d_array[x][y] / static_cast<double>(hitmax);
            image.put_pixel(x, y, Colour(foo, foo, foo));
          }
      string filename = s.output_base + "_hitpattern.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);

      // Texture
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rhitcount->i2xy(i, x, y);
            rtexture->d_array[x][y] /= rhitcount->d_array[x][y];
            if (s.threshold) rtexture->d_array[x][y] = (rtexture->d_array[x][y] > 0) ? -1.0 : 1.0;
          }
      rtexture->minmax(min, max);
      cout << "min, max texture = " << min << ' ' << max << endl;
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rtexture->i2xy(i, x, y);
            foo = 1.0 - 0.95 * (rtexture->d_array[x][y] - min) / (max-min);
            image.put_pixel(x, y, Colour(foo, foo, foo));
          }
      filename = s.output_base + "_texture.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);

      // Texture modulated by polarization magnitude
      rtexture->minmax(min, max);
      rmagnitude->minmax(mmin, mmax);
      cout << "min, max polarization magnitude = " << mmin << ", " << mmax << endl;
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rhitcount->i2xy(i, x, y);
            modtexture->d_array[x][y] = (rtexture->d_array[x][y] - min) *
              (rmagnitude->d_array[x][y] - mmin);
          }
      modtexture->minmax(min, max);
      cout << "min, max = " << min << " " << max << endl;
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rhitcount->i2xy(i, x, y);
            foo = 1.0 - 0.95 * (modtexture->d_array[x][y] - min) / (max-min);
            image.put_pixel(x, y, Colour(foo, foo, foo));
          }
      filename = s.output_base + "_mod_texture.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);

      // Temperature
      rtemperature->minmax(min, max);
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rtexture->i2xy(i, x, y);
            foo = (rtemperature->d_array[x][y] - min) / (max-min);
            image.put_pixel(x, y, palette4.Get_Colour(foo));
          }
      filename = s.output_base + "_temperature.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);

      // Magnitude
      rmagnitude->minmax(min, max);
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rtexture->i2xy(i, x, y);
            foo = (rmagnitude->d_array[x][y] - min) / (max-min);
            image.put_pixel(x, y, palette4.Get_Colour(foo));
          }
      filename = s.output_base + "_magnitude.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);


      // Temperature and polarization
      double r, g, b;
      rtemperature->minmax(tmin, tmax);
      modtexture->minmax(min, max);
      for(i = 0; i <= rtexture->max_pixel(); i++)
        if (rtexture->is_valid_pixel(i))
          {
            rtexture->i2xy(i, x, y);
            myColour = palette4.Get_Colour((rtemperature->d_array[x][y] - tmin) / (tmax - tmin));
            r = myColour.r;
            g = myColour.g;
            b = myColour.b;
            rgbOverOperator(r, g, b, 0, 0, 0,
                            1.0 - 0.95 * (modtexture->d_array[x][y] - min) / (max - min));
            image.put_pixel(x, y, Colour(r, g, b));
          }
      filename = s.output_base + "_temperature_mod_texture.tga";
      cout  << "Writing image: " << filename << endl;
      image.write_TGA(filename);
    }
}

/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file tga_image.h
 *  Classes for creation and output of TGA image files
 *
 *  Copyright (C) 2003, 2006 Max-Planck-Society
 *  \author Martin Reinecke, David Larson
 */

#ifndef PLANCK_TGA_IMAGE_H
#define PLANCK_TGA_IMAGE_H

#include <string>
#include <vector>
#include <algorithm>
#include "arr.h"

/*! \defgroup imagegroup Image creation */
/*! \{ */

/*! A very simple class for storing RGB colours. */
class Colour
  {
  public:
    float r, /*!< the red component */
          g, /*!< the green component */
          b; /*!< the blue component */

    /*! Default constructor. Does not initialize \a r, \a g and \a b. */
    Colour() {}
    /*! Initializes the colour with \a R, \a G and \a B. */
    Colour (float R, float G, float B) : r(R), g(G), b(B) {}
    /*! Multiplies all components with \a fact. */
    const Colour &operator*= (float fact)
      { r*=fact; g*=fact; b*=fact; return *this; }
    /*! Returns the sum colour of \a *this and \a c2. */
    const Colour operator+ (const Colour &c2) const
      { return Colour(r+c2.r, g+c2.g, b+c2.b); }
    /*! Returns \a *this, scaled by \a f. */
    const Colour operator* (double f) const
      { return Colour(r*f, g*f, b*f); }
  };

/*! A class for mapping floting-point values into colours. */
class Palette
  {
  private:
    std::vector<Colour> cv;
    std::vector<float> fv;

  public:
    /*! Adds a new data point \a f, with the corresponding colour \a c.
        The additions must be done in the order of ascending \a f. */
    void add (float f, const Colour &c)
      {
      fv.push_back(f);
      cv.push_back(c);
      }
    /*! Sets the palette to the predefined palette \a num. */
    void setPredefined(int num);
    /*! Returns the colour corresponding to the value \a f. The colour is
        determined by linear interpolation between neighbouring data points. */
    Colour Get_Colour (float f) const
      {
      if (f<=fv[0]) return cv[0];
      if (f>=fv[fv.size()-1]) return cv[cv.size()-1];
      int i=0;
      while (f>fv[i]) ++i;
      return cv[i-1]*((fv[i]-f)/(fv[i]-fv[i-1]))
           + cv[i]*((f-fv[i-1])/(fv[i]-fv[i-1]));
      }
  };

class Colour8
  {
  private:
    void import (const Colour &col)
      {
      using namespace std;
      r = max(0,min(255,int(col.r*256)));
      g = max(0,min(255,int(col.g*256)));
      b = max(0,min(255,int(col.b*256)));
      }

  public:
    char r,g,b;

    Colour8() {}
    Colour8 (unsigned char R, unsigned char G, unsigned char B)
      : r(R), g(G), b(B) {}
    Colour8 (const Colour &col)
      { import (col); }
    const Colour8 &operator= (const Colour &col)
      { import (col); return *this; }
    bool operator== (const Colour8 &that)
      { return (r == that.r) && (g == that.g) && (b == that.b); }
    bool operator!= (const Colour8 &that)
      { return (r != that.r) || (g != that.g) || (b != that.b); }
  };

class Font
  {
  public:
    int offset, num_chars, xpix, ypix;
    const char *data;
  };

extern const Font medium_bold_font;
extern const Font giant_font;

/*! Class for creating and storing TGA image files. */
class TGA_Image
  {
  private:
    Font font;
    arr2<Colour8> pixel;

    void write_char (int xpos, int ypos, const Colour &col, char c,
                     int scale=1);

  public:
    /*! */
    TGA_Image ();
    /*! Creates an image object with a resolution of \a xres by \a yres. */
    TGA_Image (int xres, int yres);
    /*! */
    ~TGA_Image () {}

    /*! Fills the entire image with colour \a col. */
    void fill (const Colour &col) { pixel.fill(col); }
    /*! Sets the font used for annotations to \a fnt. */
    void set_font (const Font &fnt);
    /*! Outputs the string \a text in colour \a col.
        \a xpos, \a ypos is the lower left corner;
        the font is scaled by \a scale. */
    void annotate (int xpos, int ypos, const Colour &col,
      const std::string &text, int scale=1);
    /*! Outputs the string \a text centered at position \a xpos, \a ypos
        in colour \a col. The font is scaled by \a scale. */
    void annotate_centered (int xpos, int ypos, const Colour &col,
      const std::string &text, int scale=1);
    /*! Sets the pixel \a i, \a j, to the colour \a col. */
    void put_pixel (int i, int j, const Colour &col)
      {
      if ((i>=0) && (i<pixel.size1()) && (j>=0) && (j<pixel.size2()))
        pixel[i][j] = col;
      }
    /*! Returns the colour of the pixel \a i, \a j, or black if the pixel
        lies outside of the image. */
    Colour8 get_pixel (int i, int j)
      {
      if ((i>=0) && (i<pixel.size1()) && (j>=0) && (j<pixel.size2()))
        return pixel[i][j];
      else
        return Colour8(0, 0, 0);
      }

    /*! Writes the image to \a file. */
    void write (const std::string &file) const;
  };

/*! \} */

#endif

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

/*
 *  Classes for creation and output of TGA image files
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "tga_image.h"
#include <fstream>
#include "font_data.inc"

using namespace std;

const Font medium_bold_font = { 0, 128, 7, 13, medium_bold_font_data };
const Font giant_font = { 0, 128, 9, 15, giant_font_data };

void Palette::setPredefined (int num)
  {
  fv.clear(); cv.clear();
  switch(num)
    {
    case 0:
      add(0,Colour(0,0,0));
      add(1,Colour(1,1,1));
      break;
    case 1:
      add(0,Colour(0,0,0));
      add(0.4,Colour(0,0,0.5));
      add(0.75,Colour(0,0.6,1));
      add(1,Colour(1,1,1));
      break;
    case 4:
      add(0,Colour(0,0,.5));
      add(0.15,Colour(0,0,1));
      add(0.4,Colour(0,1,1));
      add(0.7,Colour(1,1,0));
      add(0.9,Colour(1,.33,0));
      add(1,Colour(.5,0,0));
      break;
    default:
      throw Message_error("Palette #"+dataToString(num)+" not yet supported.");
    }
  }

void TGA_Image::write_char (int xpos, int ypos, const Colour &col, char c,
  int scale)
  {
  for (int i=0; i<font.xpix; ++i)
    for (int j=0; j<font.ypix; ++j)
      {
      int ofs = (c-font.offset)*font.xpix*font.ypix + j*font.xpix + i;
      if (font.data[ofs]>0)
        for (int m=0; m<scale; ++m)
          for (int n=0; n<scale; ++n)
            put_pixel(xpos+scale*i+m,ypos+scale*j+n,col);
      }
  }

TGA_Image::TGA_Image ()
  : font(medium_bold_font) {}

TGA_Image::TGA_Image (int xres, int yres)
  : font(medium_bold_font), pixel(xres,yres)
  {
  pixel.fill(Colour(0,0,0));
  }

void TGA_Image::annotate (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  for (unsigned int m=0; m<text.length(); ++m)
    write_char(xpos+m*scale*font.xpix, ypos, col, text[m],scale);
  }

void TGA_Image::annotate_centered (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  xpos-=(scale*text.length()*font.xpix)/2;
  ypos-=scale*font.ypix/2;
  annotate (xpos,ypos,col,text,scale);
  }

void TGA_Image::set_font (const Font &fnt)
  {
  font = fnt;
  }

void TGA_Image::write (const string &file) const
  {
  int xres = pixel.size1();
  int yres = pixel.size2();

  const char header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    xres%256, xres/256, yres%256, yres/256, 24, 32 };

  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file " + file);

  out.write (header, 18);

  for (int j=0; j<yres; ++j)
    for (int i=0; i<xres; ++i)
      {
      out.write(&(pixel[i][j].b),1);
      out.write(&(pixel[i][j].g),1);
      out.write(&(pixel[i][j].r),1);
      }
  }

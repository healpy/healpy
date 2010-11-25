#include <math.h>
#include <assert.h>
#include <iostream>
#include "color.h"

using namespace std;

// Page 593, Computer Graphics: Principles and Practice,
// by James D. Foley et al.  Addison-Wesley, 1990
void hsvToRgb(double& red, double& green, double& blue, double hue, double saturation, double value)
  {
  int i;
  double f, p, q, t;

  // cout << "hsv = " << hue << " " << saturation << " " << value << endl;

  // Stop the largest problems
  assert(saturation >= -0.1);
  assert(saturation <= 1.1);
  assert(value >= -0.1);
  assert(value <= 1.1);
  assert(hue >= -1.0);
  assert(hue < 361.0);

  if (hue > 360.0)
    hue -= 360.0;
  if (hue < 0.0)
    {
    hue += 360.0;
    // Case where hue used to be really small, negative.
    if (hue >= 360.0) hue = 0.0;
    }

  if (hue >= 359.9) hue = 359.9;
  if (hue < 0.0001) hue = 0.0001;

  if (saturation > 1.0) saturation = 1.0;
  if (saturation < 0.0) saturation = 0.0;
  if (value > 1.0) value = 1.0;
  if (value < 0.0) value = 0.0;

  if (saturation == 0.0)
    {
    red = value;
    green = value;
    blue = value;
    }
  else
    {
    hue = hue / 60.0;
    i = static_cast<int>(floor(hue));
    f = hue - i;
    p = value * (1.0 - saturation);
    q = value * (1.0 - saturation * f);
    t = value * (1.0 - saturation * (1.0 - f));

    assert(i >= 0);
    assert(i < 6);

    switch(i)
      {
      case 0:
        red = value;
        green = t;
        blue = p;
        break;
      case 1:
        red = q;
        green = value;
        blue = p;
        break;
      case 2:
        red = p;
        green = value;
        blue = t;
        break;
      case 3:
        red = p;
        green = q;
        blue = value;
        break;
      case 4:
        red = t;
        green = p;
        blue = value;
        break;
      case 5:
        red = value;
        green = p;
        blue = q;
        break;
      }
    }
  }


// CS 319, powerpoint presentation by John C. Hart, refers to [Porter & Duff S’84]
// This is an implementation of transparancy.
void rgbOverOperator(double& red, double& green, double& blue, double redUnder, double greenUnder, double blueUnder, double alpha)
  {
  red = red * alpha + redUnder * (1 - alpha);
  green = green * alpha + greenUnder * (1 - alpha);
  blue = blue * alpha + blueUnder * (1 - alpha);
  }

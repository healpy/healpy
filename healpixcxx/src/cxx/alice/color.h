#ifndef ALICE_COLOR
#define ALICE_COLOR

// Page 593, Computer Graphics: Principles and Practice,
// by James D. Foley et al.  Addison-Wesley, 1990
void hsvToRgb(double& red, double& green, double& blue, double hue, double saturation, double value);

// CS 319, powerpoint presentation by John C. Hart, refers to [Porter & Duff S’84]
void rgbOverOperator(double& red, double& green, double& blue, double redUnder, double greenUnder, double blueUnder, double alpha);

#endif

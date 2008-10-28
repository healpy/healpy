#ifndef SOSSKYMAP 
#define SOSSKYMAP

#include "RectSkyMap.h"

class SoSSkyMap : public RectSkyMap
{
 public:
  double d_xscale, d_yscale;

  SoSSkyMap() {};
  SoSSkyMap(int x);
  void set_size(int x);

  int is_valid_pixel(int i) const
    {
      return 1;  /* All pixels are valid */
    }

  int project(pointing p) const;
  pointing deproject(int i) const;
};

#endif // SOSSKYMAP

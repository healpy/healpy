#ifndef MOLLWEIDESKYMAP
#define MOLLWEIDESKYMAP

#include "RectSkyMap.h"

class MollweideSkyMap : public RectSkyMap
  {
  private:
    double roottwo, d_halfx, d_halfy;

  public:
    MollweideSkyMap() {}
    MollweideSkyMap(int x) { set_size(x); }
    void set_size(int x);
    bool is_valid_pixel(int i) const;
    void xy2xpyp(int x, int y, double &xp, double &yp) const;
    void xpyp2xy(double xp, double yp, int &x, int &y) const;
    int project(pointing p) const;
    pointing deproject(int i) const;
  };

#endif // MOLLWEIDESKYMAP

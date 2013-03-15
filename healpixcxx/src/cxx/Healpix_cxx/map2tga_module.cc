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
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "safe_ptr.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "rotmatrix.h"
#include "pointing.h"
#include "ls_image.h"
#include "paramfile.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"
#include "string_utils.h"

using namespace std;

namespace {

void histo_eq (arr2<float> &img, float &minv, float &maxv, arr<double> &newpos)
  {
  const int nbins=100;
  arr<int> bincnt (nbins);
  bincnt.fill(0);
  int pixels=0;

  double fact = 1./(maxv-minv);
  for (tsize i=0; i<img.size1(); ++i)
    for (tsize j=0; j<img.size2(); ++j)
      if (img[i][j]>-1e30)
        {
        img[i][j] = float((img[i][j]-minv)*fact);
        int idx = int(img[i][j]*nbins);
        idx=max(0,min(idx,nbins-1));
        ++bincnt[idx];
        ++pixels;
        }

  newpos.alloc(nbins+1);
  int accu=0;
  for (int m=0; m<nbins; ++m)
    {
    newpos[m] = double(accu)/pixels;
    accu += bincnt[m];
    }
  newpos[nbins]=1.;

  for (tsize i=0; i<img.size1(); ++i)
    for (tsize j=0; j<img.size2(); ++j)
      if (img[i][j]>-1e30)
        {
        int idx = int(img[i][j]*nbins);
        idx=max(0,min(idx,nbins-1));
        double frac = nbins*img[i][j] - idx;
        img[i][j] = float((1-frac)*newpos[idx] + frac*newpos[idx+1]);
        img[i][j] = float(minv+(maxv-minv)*img[i][j]);
        }
  }

void pro_mollw (const Healpix_Map<float> &map, double lon0, double lat0,
  int xsize, arr2<float> &img, float &minv, float &maxv, bool smooth)
  {
  int ysize=xsize/2;
  img.alloc(xsize,ysize);
  img.fill(-1e35);
  double xc=(xsize-1)/2., yc=(ysize-1)/2.;

  rotmatrix rot;
  rot.Make_CPAC_Euler_Matrix(0,-lat0,-lon0);

  minv=1e30;
  maxv=-1e30;
  for (tsize i=0; i<img.size1(); ++i)
    for (tsize j=0; j<img.size2(); ++j)
      {
      double u = 2*(i-xc)/(xc/1.02);
      double v = (j-yc)/(yc/1.02);
      bool mask = ((u*u/4 + v*v) <= 1);
      if (mask)
        {
        pointing ptg (halfpi-(asin(2/pi*(asin(v) + v*sqrt((1-v)*(1+v))))),
                      -halfpi*u/max(sqrt((1-v)*(1+v)),1e-6));
        vec3 pnt = rot.Transform(ptg.to_vec3());
        if (smooth)
          img[i][j] = map.interpolated_value(pnt);
        else
          img[i][j] = map[map.ang2pix(pnt)];
        if (!approx<double>(img[i][j],Healpix_undef))
          {
          if (img[i][j]<minv) minv=img[i][j];
          if (img[i][j]>maxv) maxv=img[i][j];
          }
        }
      }
  }

void pro_gno (const Healpix_Map<float> &map, double lon0, double lat0,
  int xsize, double delta, arr2<float> &img, float &minv, float &maxv,
  bool smooth)
  {
  rotmatrix rot;
  rot.Make_CPAC_Euler_Matrix(lon0,-lat0,0);

  double start=-(xsize/2.)*delta;
  img.alloc(xsize,xsize);
  minv=1e30;
  maxv=-1e30;
  for (tsize i=0; i<img.size1(); ++i)
    for (tsize j=0; j<img.size2(); ++j)
      {
      vec3 pnt (1,-(start+i*delta), start+j*delta);
      pnt = rot.Transform(pnt);
        if (smooth)
          img[i][j] = map.interpolated_value(pnt);
        else
          img[i][j] = map[map.ang2pix(pnt)];
      if (!approx<double>(img[i][j],Healpix_undef))
        {
        if (img[i][j]<minv) minv=img[i][j];
        if (img[i][j]>maxv) maxv=img[i][j];
        }
      }
  }

void colorbar (LS_Image &img, const Palette &pal, int xmin, int xmax,
  int ymin, int ymax, bool flippal, const arr<double> &newpos)
  {
  int nbins = newpos.size()-1;
  for (int i=xmin; i<=xmax; ++i)
    {
    double val = (double(i)-xmin)/(xmax-xmin);
    if (nbins>0)
      {
      int idx = int(val*nbins);
      idx=max(0,min(idx,nbins-1));
      double frac = nbins*val - idx;
      val = (1-frac)*newpos[idx] + frac*newpos[idx+1];
      }
    if (flippal) val=1-val;
    Colour c = pal.Get_Colour(float(val));
    for (int j=ymin; j<=ymax; ++j)
      img.put_pixel(i,j,c);
    }
  }

string conv (double val)
  {
  ostringstream os;
  if (abs(val)>100 || abs(val)<0.01)
    {
    os << setw(10) << setprecision(3) << scientific << val;
    return os.str();
    }
  os << setw(10) << setprecision(6) << fixed << val;
  return trim(os.str());
  }

} // unnamed namespace

int map2tga_module (int argc, const char **argv)
  {
  module_startup ("map2tga",argc>=2,
    "\nUsage:\n"
    "  map2tga <init object>\n\n"
    "or:\n"
    "  map2tga <input file> <output file> [-sig <int>] [-pal <int>]\n"
    "    [-xsz <int>] [-bar] [-log] [-asinh] [-lon <float>] [-lat <float>]\n"
    "    [-mul <float>] [-add <float>] [-min <float>] [-max <float>]\n"
    "    [-res <float>] [-title <string>] [-flippal] [-gnomonic]\n"
    "    [-interpol] [-equalize] [-viewer <viewer>]\n\n");

  safe_ptr<paramfile> params;
  if (argc>2)
    {
    vector<string> leading;
    map<string,string> dict;
    leading.push_back("infile");
    leading.push_back("outfile");
    parse_cmdline_classic (argc,argv,leading,dict);
    if (dict.find("gnomonic")!=dict.end())
      {
      dict["pro"]="gno";
      dict.erase("gnomonic");
      }
    params = new paramfile(dict,false);
    }
  else
    params = new paramfile(argv[1]);

  string infile = params->find<string>("infile");
  string outfile = params->find<string>("outfile");
  int colnum = params->find<int>("sig",1);
  int palnr = params->find<int>("pal",4);
  bool flippal = params->find<bool>("flippal",false);
  int xres = params->find<int>("xsz",1024);
  bool bar = params->find<bool>("bar",false);
  bool logflag = params->find<bool>("log",false);
  bool eqflag = params->find<bool>("equalize",false);
  bool asinhflag = params->find<bool>("asinh",false);
  double lon0 = degr2rad*params->find<double>("lon",0);
  double lat0 = degr2rad*params->find<double>("lat",0);
  double factor = params->find<float>("mul",1);
  double offset = params->find<float>("add",0);
  float usermin = params->find<float>("min", -1e30);
  bool min_supplied = (usermin>-1e29);
  float usermax = params->find<float>("max", 1e30);
  bool max_supplied = (usermax<1e29);
  double res = arcmin2rad*params->find<double>("res",1);
  string title = params->find<string>("title","");
  string viewer = params->find<string>("viewer","");
  bool mollpro = (params->find<string>("pro","")!="gno");
  bool interpol = params->find<bool>("interpol",false);

  Healpix_Map<float> map(0,RING);
  read_Healpix_map_from_fits(infile,map,colnum,2);
  for (int m=0; m<map.Npix(); ++m)
    {
    if (!approx<double>(map[m],Healpix_undef))
      {
      map[m] = float((map[m]+offset)*factor);
      if (logflag)
        map[m] = (map[m]<=0) ? float(Healpix_undef)
                             : float(log(double(map[m]))/ln10);
      if (asinhflag)
        map[m] = (map[m]>=0) ?
          float( log(double( map[m]+sqrt(map[m]*map[m]+1)))) :
          float(-log(double(-map[m]+sqrt(map[m]*map[m]+1))));
      if (min_supplied) if (map[m] < usermin) map[m] = usermin;
      if (max_supplied) if (map[m] > usermax) map[m] = usermax;
      }
    }

  arr2<float> imgarr;
  float minv, maxv;
  mollpro ? pro_mollw (map,lon0,lat0,xres,imgarr,minv,maxv,interpol) :
            pro_gno (map,lon0,lat0,xres,res,imgarr,minv,maxv,interpol);

  arr<double> newpos;
  if (eqflag) histo_eq(imgarr,minv,maxv,newpos);

  if (min_supplied) minv = usermin;
  if (max_supplied) maxv = usermax;
  if (maxv==minv) maxv=minv+1e-10f;

  int xsz=imgarr.size1();
  int ysz=imgarr.size2();
  int yofs=ysz;
  int scale = max(1,xsz/800);
  if (bar) ysz+=60*scale;
  if (title != "") { ysz+=50*scale; yofs+=50*scale; }
  LS_Image img(xsz,ysz);
  img.fill(Colour(1,1,1));
  img.set_font (giant_font);
  Palette pal;
  pal.setPredefined(palnr);
  if (title != "")
    img.annotate_centered(xsz/2,25*scale,Colour(0,0,0),title,scale);
  for (tsize i=0; i<imgarr.size1(); ++i)
    for (tsize j=0; j<imgarr.size2(); ++j)
      {
      if (imgarr[i][j]>-1e32)
        {
        Colour c(0.5,0.5,0.5);
        if (!approx<double>(imgarr[i][j],Healpix_undef))
          {
          int col = int((imgarr[i][j]-minv)/(maxv-minv)*256);
          col = min(255,max(col,0));
          float colfrac = (imgarr[i][j]-minv)/(maxv-minv);
          if (flippal) colfrac = 1-colfrac;
          c = pal.Get_Colour(colfrac);
          }
        img.put_pixel(i,yofs-j-1,c);
        }
      }

  if (bar)
    {
    colorbar (img,pal,xsz/10,(xsz*9)/10,ysz-40*scale,ysz-20*scale,flippal,
      newpos);
    img.set_font (medium_bold_font);
    img.annotate_centered (xsz/20,ysz-30*scale,Colour(0,0,0),conv(minv),scale);
    img.annotate_centered ((xsz*19)/20,ysz-30*scale,Colour(0,0,0),conv(maxv),
      scale);
    }
  img.write_TGA(outfile);

  if (viewer!="")
    {
    int retcode = system((viewer+" "+outfile).c_str());
    if (retcode != 0)
      cout << "Warning: viewer return code was " << retcode << endl;
    }

  return 0;
  }

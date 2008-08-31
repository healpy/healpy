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
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "rotmatrix.h"
#include "pointing.h"
#include "tga_image.h"
#include "paramfile.h"

using namespace std;

void histo_eq (arr2<float> &img, float &minv, float &maxv, arr<double> &newpos)
  {
  const int nbins=100;
  arr<int> bincnt (nbins);
  bincnt.fill(0);
  int pixels=0;

  double fact = 1./(maxv-minv);
  for (int i=0; i<img.size1(); ++i)
    for (int j=0; j<img.size2(); ++j)
      if (img[i][j]>-1e30)
        {
        img[i][j] = (img[i][j]-minv)*fact;
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

  for (int i=0; i<img.size1(); ++i)
    for (int j=0; j<img.size2(); ++j)
      if (img[i][j]>-1e30)
        {
        int idx = int(img[i][j]*nbins);
        idx=max(0,min(idx,nbins-1));
        double frac = nbins*img[i][j] - idx;
        img[i][j] = (1-frac)*newpos[idx] + frac*newpos[idx+1];
        img[i][j] = minv+(maxv-minv)*img[i][j];
        }
  }

void pro_mollw (const Healpix_Map<float> &map, double lon0, double lat0,
  int xsize, arr2<float> &img, float &minv, float &maxv, bool smooth)
  {
  int ysize=xsize/2;
  img.alloc(xsize,ysize);
  img.fill(-1e35);
  double xc=(xsize-1)/2., yc=(ysize-1)/2.;
  double lon0rad = lon0*degr2rad;
  double lat0rad = lat0*degr2rad;

  rotmatrix rot;
  rot.Make_CPAC_Euler_Matrix(0,-lat0rad,-lon0rad);

  minv=1e30;
  maxv=-1e30;
  for (int i=0; i<img.size1(); ++i)
    for (int j=0; j<img.size2(); ++j)
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
          img[i][j]=map[map.ang2pix(pnt)];
        if (!approx<double>(img[i][j],Healpix_undef))
          {
          if (img[i][j]<minv) minv=img[i][j];
          if (img[i][j]>maxv) maxv=img[i][j];
          }
        }
      }
  }

void pro_gno (const Healpix_Map<float> &map, double lon0, double lat0,
  int xsize, double resgrid, arr2<float> &img, float &minv, float &maxv,
  bool smooth)
  {
  double lon0rad = lon0*degr2rad;
  double lat0rad = lat0*degr2rad;

  rotmatrix rot;
  rot.Make_CPAC_Euler_Matrix(lon0rad,-lat0rad,0);

  double delta=resgrid*degr2rad/60.;
  double start=-(xsize/2.)*delta;
  img.alloc(xsize,xsize);
  minv=1e30;
  maxv=-1e30;
  for (int i=0; i<img.size1(); ++i)
    for (int j=0; j<img.size2(); ++j)
      {
      vec3 pnt (1,-(start+i*delta), start+j*delta);
      pnt = rot.Transform(pnt);
        if (smooth)
          img[i][j] = map.interpolated_value(pnt);
        else
          img[i][j]=map[map.ang2pix(pnt)];
      if (!approx<double>(img[i][j],Healpix_undef))
        {
        if (img[i][j]<minv) minv=img[i][j];
        if (img[i][j]>maxv) maxv=img[i][j];
        }
      }
  }

void colorbar (TGA_Image &img, const Palette &pal, int xmin, int xmax,
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
    Colour c = pal.Get_Colour(val);
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

void usage()
  {
  cout <<
    "\nUsage:\n"
    "  map2tga <parameter file>\n\n"
    "or:\n"
    "  map2tga <input file> <output file> [-sig <int>] [-pal <int>]\n"
    "    [-xsz <int>] [-bar] [-log] [-asinh] [-lon <float>] [-lat <float>]\n"
    "    [-mul <float>] [-add <float>] [-min <float>] [-max <float>]\n"
    "    [-res <float>] [-title <string>] [-flippal] [-gnomonic]\n"
    "    [-interpol] [-equalize] [-viewer <viewer>]\n\n";
  throw Message_error();
  }

int main (int argc, const char **argv)
  {
PLANCK_DIAGNOSIS_BEGIN
  announce ("map2tga");
  if (argc<2) usage();

  string infile = "";
  string outfile = "";
  string title = "";
  int xres = 1024;
  bool bar = false, mollpro=true;
  double lat0=0, lon0=0, res=1;
  double usermin=-1e30, usermax=1e30, offset=0, factor=1;
  bool min_supplied=false, max_supplied=false;
  bool logflag=false, asinhflag=false, eqflag=false;
  int palnr = 4;
  int colnum=1;
  bool flippal = false;
  bool interpol = false;
  string viewer = "";

  if (argc>2)
    {
    infile = argv[1];
    outfile = argv[2];
    int m=3;
    while (m<argc)
      {
      int mstart = m;
      string arg = argv[m];
      if (arg == "-sig") { stringToData(argv[m+1],colnum); m+=2; }
      if (arg == "-pal") { stringToData(argv[m+1],palnr); m+=2; }
      if (arg == "-xsz") { stringToData(argv[m+1],xres); m+=2; }
      if (arg == "-bar") { bar=true; ++m; }
      if (arg == "-log") { logflag=true; ++m; }
      if (arg == "-equalize") { eqflag=true; ++m; }
      if (arg == "-asinh") { asinhflag=true; ++m; }
      if (arg == "-lon") { stringToData(argv[m+1],lon0); m+=2; }
      if (arg == "-lat") { stringToData(argv[m+1],lat0); m+=2; }
      if (arg == "-mul") { stringToData(argv[m+1],factor); m+=2; }
      if (arg == "-add") { stringToData(argv[m+1],offset); m+=2; }
      if (arg == "-min")
        {
        stringToData(argv[m+1],usermin);
        min_supplied=true;
        m+=2;
        }
      if (arg == "-max")
        {
        stringToData(argv[m+1],usermax);
        max_supplied=true;
        m+=2;
        }
      if (arg == "-res") { stringToData(argv[m+1],res); m+=2; }
      if (arg == "-title") { title = argv[m+1]; m+=2; }
      if (arg == "-viewer") { viewer = argv[m+1]; m+=2; }
      if (arg == "-flippal") { flippal=true; ++m; }
      if (arg == "-gnomonic") { mollpro=false; ++m; }
      if (arg == "-interpol") { interpol=true; ++m; }

      if (mstart==m)
        {
        cout << "unrecognized option: " + arg << endl;
        usage();
        }
      }
    }
  else
    {
    paramfile params (argv[1]);
    infile = params.find<string>("infile");
    outfile = params.find<string>("outfile");
    colnum = params.find<int>("sig",colnum);
    palnr = params.find<int>("pal",palnr);
    flippal = params.find<bool>("flippal",flippal);
    xres = params.find<int>("xsz",xres);
    bar = params.find<bool>("bar",bar);
    logflag = params.find<bool>("log",logflag);
    eqflag = params.find<bool>("equalize",eqflag);
    asinhflag = params.find<bool>("asinh",asinhflag);
    lon0 = params.find<double>("lon",lon0);
    lat0 = params.find<double>("lat",lat0);
    factor = params.find<double>("mul",factor);
    offset = params.find<double>("add",offset);
    usermin = params.find<double>("min", usermin);
    if (usermin>-1e29) min_supplied = true;
    usermax = params.find<double>("max", usermax);
    if (usermax<1e29) max_supplied = true;
    res = params.find<double>("res",res);
    title = params.find<string>("title","");
    viewer = params.find<string>("viewer","");
    string tmp = params.find<string>("pro","");
    if (tmp == "gno") mollpro=false;
    interpol = params.find<bool>("interpol");
    }

  Healpix_Map<float> map(0,RING);
  read_Healpix_map_from_fits(infile,map,colnum,2);
  for (int m=0; m<map.Npix(); ++m)
    {
    if (!approx<double>(map[m],Healpix_undef))
      {
      map[m] = (map[m]+offset)*factor;
      if (logflag)
        {
        if (map[m]<=0)
          map[m] = Healpix_undef;
        else
          map[m] = log(double(map[m]))/ln10;
        }
      if (asinhflag)
        {
        if (map[m]>=0)
          map[m] = log(double(map[m]+sqrt(map[m]*map[m]+1)));
        else
          map[m] = -log(double(-map[m]+sqrt(map[m]*map[m]+1)));
        }
      if (min_supplied) if (map[m] < usermin) map[m] = usermin;
      if (max_supplied) if (map[m] > usermax) map[m] = usermax;
      }
    }

  arr2<float> imgarr;
  float minv, maxv;
  if (mollpro)
    pro_mollw (map,lon0,lat0,xres,imgarr,minv,maxv,interpol);
  else
    pro_gno (map,lon0,lat0,xres,res,imgarr,minv,maxv,interpol);

  arr<double> newpos;
  if (eqflag) histo_eq(imgarr,minv,maxv,newpos);

  if (min_supplied) minv = usermin;
  if (max_supplied) maxv = usermax;
  if (maxv==minv) maxv=minv+1e-10;

  int xsz=imgarr.size1();
  int ysz=imgarr.size2();
  int yofs=ysz;
  int scale = max(1,xsz/800);
  if (bar) ysz+=60*scale;
  if (title != "") { ysz+=50*scale; yofs+=50*scale; }
  TGA_Image img(xsz,ysz);
  img.fill(Colour(1,1,1));
  img.set_font (giant_font);
  Palette pal;
  pal.setPredefined(palnr);
  if (title != "")
    img.annotate_centered(xsz/2,25*scale,Colour(0,0,0),title,scale);
  for (int i=0; i<imgarr.size1(); ++i)
    for (int j=0; j<imgarr.size2(); ++j)
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
  img.write(outfile);

  if (viewer!="")
    system ((viewer+" "+outfile).c_str());
PLANCK_DIAGNOSIS_END
  }

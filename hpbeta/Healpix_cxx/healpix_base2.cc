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

#include "healpix_base2.h"
#include "cxxutils.h"
#include "pointing.h"
#include "arr.h"
#include "geom_utils.h"

using namespace std;

short Healpix_Base2::ctab[];
short Healpix_Base2::utab[];

Healpix_Base2::Tablefiller::Tablefiller()
  {
  for (int m=0; m<0x100; ++m)
    {
    ctab[m] =
         (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
      | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
    utab[m] =
         (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
      | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
  }

Healpix_Base2::Tablefiller Healpix_Base2::Filler;

const int Healpix_Base2::jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
const int Healpix_Base2::jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

int64 Healpix_Base2::npix2nside (int64 npix)
  {
  int64 res=isqrt(npix/12);
  planck_assert (npix==res*res*12, "npix2nside: invalid argument");
  return res;
  }

void Healpix_Base2::nest2xyf (int64 pix, int &ix, int &iy, int &face_num)
  const
  {
  face_num = pix>>(2*order_);
  pix &= (npface_-1);
  int32 raw = ((pix&0x555500000000ull)>>16) 
             | ((pix&0x5555000000000000ull)>>31)
             | (pix&0x5555)
             | ((pix&0x55550000)>>15);
  ix =  ctab[raw&0xff]
     | (ctab[(raw>>8)&0xff]<<4)
     | (ctab[(raw>>16)&0xff]<<16)
     | (ctab[(raw>>24)&0xff]<<20);
  pix >>= 1;
  raw = ((pix&0x555500000000ull)>>16) 
             | ((pix&0x5555000000000000ull)>>31)
             | (pix&0x5555)
             | ((pix&0x55550000)>>15);
  iy =  ctab[raw&0xff]
     | (ctab[(raw>>8)&0xff]<<4)
     | (ctab[(raw>>16)&0xff]<<16)
     | (ctab[(raw>>24)&0xff]<<20);
  }

int64 Healpix_Base2::xyf2nest (int ix, int iy, int face_num) const
  {
  return (int64(face_num)<<(2*order_)) +
    (   (int64(utab[ ix     &0xff]))
      | (int64(utab[(ix>> 8)&0xff])<<16)
      | (int64(utab[(ix>>16)&0xff])<<32)
      | (int64(utab[(ix>>24)&0xff])<<48)
      | (int64(utab[ iy     &0xff])<<1)
      | (int64(utab[(iy>> 8)&0xff])<<17)
      | (int64(utab[(iy>>16)&0xff])<<33)
      | (int64(utab[(iy>>24)&0xff])<<49) ); 
  }

void Healpix_Base2::ring2xyf (int64 pix, int &ix, int &iy, int &face_num)
  const
  {
  int64 iring, iphi, kshift, nr;

  int64 nl2 = 2*nside_;

  if (pix<ncap_) // North Polar cap
    {
    iring = int64(0.5*(1+isqrt(1+2*pix))); //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    face_num=0;
    int64 tmp = iphi-1;
    if (tmp>=(2*iring))
      {
      face_num=2;
      tmp-=2*iring;
      }
    if (tmp>=iring) ++face_num;
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    int64 ip = pix - ncap_;
    if (order_>=0)
      {
      iring = (ip>>(order_+2)) + nside_; // counted from North pole
      iphi  = (ip&(4*nside_-1)) + 1;
      }
    else
      {
      iring = (ip/(4*nside_)) + nside_; // counted from North pole
      iphi  = (ip%(4*nside_)) + 1;
      }
    kshift = (iring+nside_)&1;
    nr = nside_;
    int64 ire = iring-nside_+1;
    int64 irm = nl2+2-ire;
    int64 ifm, ifp;
    if (order_>=0)
      {
      ifm = (iphi - ire/2 + nside_ -1) >> order_;
      ifp = (iphi - irm/2 + nside_ -1) >> order_;
      }
    else
      {
      ifm = (iphi - ire/2 + nside_ -1) / nside_;
      ifp = (iphi - irm/2 + nside_ -1) / nside_;
      }
    if (ifp == ifm) // faces 4 to 7
      face_num = (ifp==4) ? 4 : ifp+4;
    else if (ifp<ifm) // (half-)faces 0 to 3
      face_num = ifp;
    else // (half-)faces 8 to 11
      face_num = ifm + 8;
    }
  else // South Polar cap
    {
    int64 ip = npix_ - pix;
    iring = int64(0.5*(1+isqrt(2*ip-1))); //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    face_num=8;
    int64 tmp = iphi-1;
    if (tmp>=(2*nr))
      {
      face_num=10;
      tmp-=2*nr;
      }
    if (tmp>=nr) ++face_num;
    }

  int64 irt = iring - (jrll[face_num]*nside_) + 1;
  int64 ipt = 2*iphi- jpll[face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  ix =  (ipt-irt) >>1;
  iy =(-(ipt+irt))>>1;
  }

int64 Healpix_Base2::xyf2ring (int ix, int iy, int face_num) const
  {
  int64 nl4 = 4*nside_;
  int64 jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  int64 nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = npix_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  int64 jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }

double Healpix_Base2::ring2z (int64 ring) const
  {
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
  }

int64 Healpix_Base2::pix2ring (int64 pix) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      return int64(0.5*(1+isqrt(1+2*pix))); //counted from North pole
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      int64 ip  = pix - ncap_;
      return ip/(4*nside_) + nside_; // counted from North pole
      }
    else // South Polar cap
      {
      int64 ip = npix_ - pix;
      return 4*nside_ - int64(0.5*(1+isqrt(2*ip-1))); //counted from South pole
      }
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);
    return (int64(jrll[face_num])<<order_) - ix - iy - 1;
    }
  }

int64 Healpix_Base2::nest2ring (int64 pix) const
  {
  planck_assert(order_>=0, "nest2ring: need hierarchical map");
  int ix, iy, face_num;
  nest2xyf (pix, ix, iy, face_num);
  return xyf2ring (ix, iy, face_num);
  }

int64 Healpix_Base2::ring2nest (int64 pix) const
  {
  planck_assert(order_>=0, "ring2nest: need hierarchical map");
  int ix, iy, face_num;
  ring2xyf (pix, ix, iy, face_num);
  return xyf2nest (ix, iy, face_num);
  }

int64 Healpix_Base2::nest2peano (int64 pix) const
  {
  static const unsigned char subpix[8][4] = {
    { 0, 1, 3, 2 }, { 3, 0, 2, 1 }, { 2, 3, 1, 0 }, { 1, 2, 0, 3 },
    { 0, 3, 1, 2 }, { 1, 0, 2, 3 }, { 2, 1, 3, 0 }, { 3, 2, 0, 1 } };
  const unsigned char subpath[8][4] = {
    { 4, 0, 6, 0 }, { 7, 5, 1, 1 }, { 2, 4, 2, 6 }, { 3, 3, 7, 5 },
    { 0, 2, 4, 4 }, { 5, 1, 5, 3 }, { 6, 6, 0, 2 }, { 1, 7, 3, 7 } };
  static const unsigned char face2path[12] = {
    2, 5, 2, 5, 3, 6, 3, 6, 2, 3, 2, 3 };
  static const unsigned char face2peanoface[12] = {
    0, 5, 6, 11, 10, 1, 4, 7, 2, 3, 8, 9 };

  int face = pix>>(2*order_);
  unsigned char path = face2path[face];
  int64 result = 0;

  for (int shift=2*order_-2; shift>=0; shift-=2)
    {
    unsigned char spix = (pix>>shift) & 0x3;
    result <<= 2;
    result |= subpix[path][spix];
    path=subpath[path][spix];
    }

  return result + (int64(face2peanoface[face])<<(2*order_));
  }

int64 Healpix_Base2::peano2nest (int64 pix) const
  {
  static const unsigned char subpix[8][4] = {
    { 0, 1, 3, 2 }, { 1, 3, 2, 0 }, { 3, 2, 0, 1 }, { 2, 0, 1, 3 },
    { 0, 2, 3, 1 }, { 1, 0, 2, 3 }, { 3, 1, 0, 2 }, { 2, 3, 1, 0 } };
  static const unsigned char subpath[8][4] = {
    { 4, 0, 0, 6 }, { 5, 1, 1, 7 }, { 6, 2, 2, 4 }, { 7, 3, 3, 5 },
    { 0, 4, 4, 2 }, { 1, 5, 5, 3 }, { 2, 6, 6, 0 }, { 3, 7, 7, 1 } };
  static const unsigned char face2path[12] = {
    2, 6, 2, 3, 3, 5, 2, 6, 2, 3, 3, 5 };
  static const unsigned char peanoface2face[12] = {
    0, 5, 8, 9, 6, 1, 2, 7, 10, 11, 4, 3 };

  int face = pix>>(2*order_);
  unsigned char path = face2path[face];
  int64 result = 0;

  for (int shift=2*order_-2; shift>=0; shift-=2)
    {
    unsigned char spix = (pix>>shift) & 0x3;
    result <<= 2;
    result |= subpix[path][spix];
    path=subpath[path][spix];
    }

  return result + (int64(peanoface2face[face])<<(2*order_));
  }

int64 Healpix_Base2::ang2pix_z_phi (double z, double phi) const
  {
  double za = abs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; // in [0,4)

  if (scheme_==RING)
    {
    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*z*0.75;
      int64 jp = int64(temp1-temp2); // index of  ascending edge line
      int64 jm = int64(temp1+temp2); // index of descending edge line

      // ring number counted from z=2/3
      int64 ir = nside_ + 1 + jp - jm; // in {1,2n+1}
      int kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

      int64 ip = (jp+jm-nside_+kshift+1)/2; // in {0,4n-1}
      ip = imodulo<int64>(ip,4*nside_);

      return ncap_ + (ir-1)*4*nside_ + ip;
      }
    else  // North & South polar caps
      {
      double tp = tt-int(tt);
      double tmp = nside_*sqrt(3*(1-za));

      int64 jp = int64(tp*tmp); // increasing edge line index
      int64 jm = int64((1.0-tp)*tmp); // decreasing edge line index

      int64 ir = jp+jm+1; // ring number counted from the closest pole
      int64 ip = int64(tt*ir); // in {0,4*ir-1}
      ip = imodulo<int64>(ip,4*ir);

      if (z>0)
        return 2*ir*(ir-1) + ip;
      else
        return npix_ - 2*ir*(ir+1) + ip;
      }
    }
  else // scheme_ == NEST
    {
    int face_num, ix, iy;

    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*(z*0.75);
      int64 jp = int64(temp1-temp2); // index of  ascending edge line
      int64 jm = int64(temp1+temp2); // index of descending edge line
      int64 ifp = jp >> order_;  // in {0,4}
      int64 ifm = jm >> order_;
      if (ifp == ifm)           // faces 4 to 7
        face_num = (ifp==4) ? 4: ifp+4;
      else if (ifp < ifm)       // (half-)faces 0 to 3
        face_num = ifp;
      else                      // (half-)faces 8 to 11
        face_num = ifm + 8;

      ix = jm & (nside_-1);
      iy = nside_ - (jp & (nside_-1)) - 1;
      }
    else // polar region, za > 2/3
      {
      int ntt = int(tt);
      double tp = tt-ntt;
      double tmp = nside_*sqrt(3*(1-za));

      int64 jp = int64(tp*tmp); // increasing edge line index
      int64 jm = int64((1.0-tp)*tmp); // decreasing edge line index
      if (jp>=nside_) jp = nside_-1; // for points too close to the boundary
      if (jm>=nside_) jm = nside_-1;
      if (z >= 0)
        {
        face_num = ntt;  // in {0,3}
        ix = nside_ - jm - 1;
        iy = nside_ - jp - 1;
        }
      else
        {
        face_num = ntt + 8; // in {8,11}
        ix =  jp;
        iy =  jm;
        }
      }

    return xyf2nest(ix,iy,face_num);
    }
  }

void Healpix_Base2::pix2ang_z_phi (int64 pix, double &z, double &phi) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      {
      int64 iring = int64(0.5*(1+isqrt(1+2*pix))); //counted from North pole
      int64 iphi  = (pix+1) - 2*iring*(iring-1);

      z = 1.0 - (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      int64 ip  = pix - ncap_;
      int64 iring = ip/(4*nside_) + nside_; // counted from North pole
      int64 iphi  = ip%(4*nside_) + 1;
      // 1 if iring+nside is odd, 1/2 otherwise
      double fodd = ((iring+nside_)&1) ? 1 : 0.5;

      int64 nl2 = 2*nside_;
      z = (nl2-iring)*fact1_;
      phi = (iphi-fodd) * pi/nl2;
      }
    else // South Polar cap
      {
      int64 ip = npix_ - pix;
      int64 iring = int64(0.5*(1+isqrt(2*ip-1))); //counted from South pole
      int64 iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

      z = -1.0 + (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    }
  else
    {
    int64 nl4 = nside_*4;

    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);

    int64 jr = (int64(jrll[face_num])<<order_) - ix - iy - 1;

    int64 nr;
    int kshift;
    if (jr<nside_)
      {
      nr = jr;
      z = 1 - nr*nr*fact2_;
      kshift = 0;
      }
    else if (jr > 3*nside_)
      {
      nr = nl4-jr;
      z = nr*nr*fact2_ - 1;
      kshift = 0;
      }
    else
      {
      nr = nside_;
      z = (2*nside_-jr)*fact1_;
      kshift = (jr-nside_)&1;
      }

    int64 jp = (jpll[face_num]*nr + ix -iy + 1 + kshift) / 2;
    if (jp>nl4) jp-=nl4;
    if (jp<1) jp+=nl4;

    phi = (jp-(kshift+1)*0.5)*(halfpi/nr);
    }
  }

void Healpix_Base2::neighbors (int64 pix, fix_arr<int64,8> &result) const
  {
  static const int xoffset[] = { -1,-1, 0, 1, 1, 1, 0,-1 };
  static const int yoffset[] = {  0, 1, 1, 1, 0,-1,-1,-1 };
  static const int facearray[][12] =
        { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },   // S
          {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },   // SE
          { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },   // E
          {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },   // SW
          {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },   // center
          {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },   // NE
          { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },   // W
          {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },   // NW
          {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } }; // N
  static const int swaparray[][12] =
        { {  0,0,0,0,0,0,0,0,3,3,3,3 },   // S
          {  0,0,0,0,0,0,0,0,6,6,6,6 },   // SE
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // E
          {  0,0,0,0,0,0,0,0,5,5,5,5 },   // SW
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // center
          {  5,5,5,5,0,0,0,0,0,0,0,0 },   // NE
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // W
          {  6,6,6,6,0,0,0,0,0,0,0,0 },   // NW
          {  3,3,3,3,0,0,0,0,0,0,0,0 } }; // N

  int ix, iy, face_num;
  (scheme_==RING) ?
    ring2xyf(pix,ix,iy,face_num) : nest2xyf(pix,ix,iy,face_num);

  const int64 nsm1 = nside_-1;
  if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
    {
    if (scheme_==RING)
      for (int m=0; m<8; ++m)
        result[m] = xyf2ring(ix+xoffset[m],iy+yoffset[m],face_num);
    else
      for (int m=0; m<8; ++m)
        result[m] = xyf2nest(ix+xoffset[m],iy+yoffset[m],face_num);
    }
  else
    {
    for (int i=0; i<8; ++i)
      {
      int x=ix+xoffset[i];
      int y=iy+yoffset[i];
      int nbnum=4;
      if (x<0)
        { x+=nside_; nbnum-=1; }
      else if (x>=nside_)
        { x-=nside_; nbnum+=1; }
      if (y<0)
        { y+=nside_; nbnum-=3; }
      else if (y>=nside_)
        { y-=nside_; nbnum+=3; }

      int f = facearray[nbnum][face_num];
      if (f>=0)
        {
        if (swaparray[nbnum][face_num]&1) x=nside_-x-1;
        if (swaparray[nbnum][face_num]&2) y=nside_-y-1;
        if (swaparray[nbnum][face_num]&4) std::swap(x,y);
        result[i] = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
        }
      else
        result[i] = -1;
      }
    }
  }

double Healpix_Base2::max_pixrad() const
  {
  vec3 va,vb;
  va.set_z_phi (2./3., pi/(4*nside_));
  double t1 = 1.-1./nside_;
  t1*=t1;
  vb.set_z_phi (1-t1/3, 0);
  return v_angle(va,vb);
  }

/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Utilities for conversion between coordinates, Morton, and Peano indices
 *
 *  Copyright (C) 2015-2020 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ducc0/math/space_filling.h"

namespace ducc0 {

#ifndef DUCC0_USE_PDEP_PEXT

#if 1

namespace {

const uint16_t utab[] = {
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
X(0),X(1),X(4),X(5)
#undef X
#undef Y
#undef Z
};
static const uint16_t ctab[] = {
#define Z(a) a,a+1,a+256,a+257
#define Y(a) Z(a),Z(a+2),Z(a+512),Z(a+514)
#define X(a) Y(a),Y(a+4),Y(a+1024),Y(a+1028)
X(0),X(8),X(2048),X(2056)
#undef X
#undef Y
#undef Z
};

}

uint32_t spread_bits_2D_32 (uint32_t v)
  {
  using I=uint32_t;
  return  I(utab[ v     &0xff])     | (I(utab[(v>> 8)&0xff])<<16);
  }
uint64_t spread_bits_2D_64 (uint64_t v)
  {
  using I=uint64_t;
  return  I(utab[ v     &0xff])      | (I(utab[(v>> 8)&0xff])<<16)
       | (I(utab[(v>>16)&0xff])<<32) | (I(utab[(v>>24)&0xff])<<48);
  }
uint32_t block2morton2D_32 (uint32_t v)
  {
  using I=uint32_t;
  return  I(utab[ v     &0xff])     | (I(utab[(v>> 8)&0xff])<<16)
       | (I(utab[(v>>16)&0xff])<<1) | (I(utab[(v>>24)&0xff])<<17);
  }
uint32_t coord2morton2D_32 (std::array<uint32_t,2> xy)
  {
  using I=uint32_t;
  return  I(utab[xy[0]&0xff])     | (I(utab[(xy[0]>>8)&0xff])<<16)
       | (I(utab[xy[1]&0xff])<<1) | (I(utab[(xy[1]>>8)&0xff])<<17);
  }
uint32_t morton2block2D_32 (uint32_t v)
  {
  using I=uint32_t;
  I raw1 = v&0x55555555, raw2 = (v>>1)&0x55555555;
  raw1|=raw1>>15;
  raw2|=raw2>>15;
  return I(ctab[raw1&0xff]    ) | I(ctab[(raw1>>8)&0xff]<< 4)
       | I(ctab[raw2&0xff]<<16) | I(ctab[(raw2>>8)&0xff]<<20);
  }
std::array<uint32_t,2> morton2coord2D_32 (uint32_t v)
  {
  using I=uint32_t;
  I raw1 = v&0x55555555u, raw2 = (v>>1)&0x55555555u;
  raw1|=raw1>>15;
  raw2|=raw2>>15;
  return {I(ctab[raw1&0xff]) | I(ctab[(raw1>>8)&0xff]<<4),
          I(ctab[raw2&0xff]) | I(ctab[(raw2>>8)&0xff]<<4)};
  }
uint64_t block2morton2D_64 (uint64_t v)
  {
  using I=uint64_t;
  return  I(utab[ v     &0xff])      | (I(utab[(v>> 8)&0xff])<<16)
       | (I(utab[(v>>16)&0xff])<<32) | (I(utab[(v>>24)&0xff])<<48)
       | (I(utab[(v>>32)&0xff])<< 1) | (I(utab[(v>>40)&0xff])<<17)
       | (I(utab[(v>>48)&0xff])<<33) | (I(utab[(v>>56)&0xff])<<49);
  }
uint64_t coord2morton2D_64 (std::array<uint64_t,2> xy)
  {
  using I=uint64_t;
  return  I(utab[ xy[0]     &0xff])      | (I(utab[(xy[0]>> 8)&0xff])<<16)
       | (I(utab[(xy[0]>>16)&0xff])<<32) | (I(utab[(xy[0]>>24)&0xff])<<48)
       | (I(utab[ xy[1]     &0xff])<< 1) | (I(utab[(xy[1]>> 8)&0xff])<<17)
       | (I(utab[(xy[1]>>16)&0xff])<<33) | (I(utab[(xy[1]>>24)&0xff])<<49);
  }
uint64_t morton2block2D_64 (uint64_t v)
  {
  using I=uint64_t;
  I raw1 = v&0x5555555555555555, raw2 = (v>>1)&0x5555555555555555;
  raw1|=raw1>>15;
  raw2|=raw2>>15;
  return I(ctab[ raw1     &0xff])      | (I(ctab[(raw1>> 8)&0xff])<< 4)
      | (I(ctab[(raw1>>32)&0xff])<<16) | (I(ctab[(raw1>>40)&0xff])<<20)
      | (I(ctab[ raw2     &0xff])<<32) | (I(ctab[(raw2>> 8)&0xff])<<36)
      | (I(ctab[(raw2>>32)&0xff])<<48) | (I(ctab[(raw2>>40)&0xff])<<52);
  }
std::array<uint64_t,2> morton2coord2D_64 (uint64_t v)
  {
  using I=uint64_t;
  I raw1 = v&0x5555555555555555, raw2 = (v>>1)&0x5555555555555555;
  raw1|=raw1>>15;
  raw2|=raw2>>15;
  return { I(ctab[ raw1     &0xff])      | (I(ctab[(raw1>> 8)&0xff])<< 4)
        | (I(ctab[(raw1>>32)&0xff])<<16) | (I(ctab[(raw1>>40)&0xff])<<20),
           I(ctab[ raw2     &0xff])      | (I(ctab[(raw2>> 8)&0xff])<< 4)
        | (I(ctab[(raw2>>32)&0xff])<<16) | (I(ctab[(raw2>>40)&0xff])<<20)};
  }

#else
// alternative implementation, usually slower

static inline uint64_t spread2D_64 (uint64_t v)
  {
  v&=0xffffffffu;
  v = (v|(v<<16)) & 0x0000ffff0000ffffu;
  v = (v|(v<< 8)) & 0x00ff00ff00ff00ffu;
  v = (v|(v<< 4)) & 0x0f0f0f0f0f0f0f0fu;
  v = (v|(v<< 2)) & 0x3333333333333333u;
  v = (v|(v<< 1)) & 0x5555555555555555u;
  return v;
  }
static inline uint64_t compress2D_64 (uint64_t v)
  {
  v&=0x5555555555555555u;
  v = (v|(v>> 1)) & 0x3333333333333333u;
  v = (v|(v>> 2)) & 0x0f0f0f0f0f0f0f0fu;
  v = (v|(v>> 4)) & 0x00ff00ff00ff00ffu;
  v = (v|(v>> 8)) & 0x0000ffff0000ffffu;
  v = (v|(v>>16)) & 0x00000000ffffffffu;
  return v;
  }

uint32_t block2morton2D_32 (uint32_t v)
  { uint64_t t=spread2D_64(v); return (t | (t>>31)) & 0xffffffffu; }
uint32_t coord2morton2D_32 (uint32_t x, uint32_t y)
  {
  uint64_t t=spread2D_64((x&0xffff)|(y<<16));
  return (t | (t>>31)) & 0xffffffffu;
  }
uint32_t morton2block2D_32 (uint32_t v)
  { uint64_t t=v; t|=t<<31; t=compress2D_64(t); return t; }
void morton2coord2D_32 (uint32_t v, uint32_t *x, uint32_t *y)
  {
  uint64_t t=v; t|=t<<31; t=compress2D_64(t);
  *x = t&0xffff;
  *y = t>>16;
  }
uint64_t block2morton2D_64 (uint64_t v)
  { return spread2D_64(v) | (spread2D_64(v>>32)<<1); }
uint64_t coord2morton2D_64 (uint64_t x, uint64_t y)
  { return spread2D_64(x) | (spread2D_64(y)<<1); }
uint64_t morton2block2D_64 (uint64_t v)
  { return compress2D_64(v) | (compress2D_64(v>>1)<<32); }
void morton2coord2D_64 (uint64_t v, uint64_t *x, uint64_t *y)
  { *x = compress2D_64(v); *y = compress2D_64(v>>1); }

#endif

namespace {

inline uint32_t spread3D_32 (uint32_t v)
  {
  v&=0x3ff;
  v = (v|(v<< 8)|(v<<16)) & 0x0f00f00fu;
  v = (v|(v<< 4)) & 0xc30c30c3u;
  v = (v|(v<< 2)) & 0x49249249u;
  return v;
  }
inline uint32_t compress3D_32 (uint32_t v)
  {
  v&=0x9249249u;
  v = (v|(v>> 2)) & 0xc30c30c3u;
  v = (v|(v>> 4)) & 0x0f00f00fu;
  v = (v|(v>> 8)|(v>>16)) & 0x3ffu;
  return v;
  }
inline uint64_t spread3D_64 (uint64_t v)
  {
  v&=0x1fffff;
  v = (v|(v<<16)|(v<<32)) & 0x00ff0000ff0000ffu;
  v = (v|(v<< 8)) & 0xf00f00f00f00f00fu;
  v = (v|(v<< 4)) & 0x30c30c30c30c30c3u;
  v = (v|(v<< 2)) & 0x9249249249249249u;
  return v;
  }
inline uint64_t compress3D_64 (uint64_t v)
  {
  v&=0x1249249249249249u;
  v=(v|(v>> 2)) & 0x30c30c30c30c30c3u;
  v=(v|(v>> 4)) & 0xf00f00f00f00f00fu;
  v=(v|(v>> 8)) & 0x00ff0000ff0000ffu;
  v=(v|(v>>16)|(v>>32)) & 0x1fffffu;
  return v;
  }

}

uint32_t block2morton3D_32 (uint32_t v)
  {
  uint32_t v2=v&0xfffff;
  uint64_t v3=spread3D_64(v2);
  v3=(v3|(v3>>29))&0x1fffffff;
  return uint32_t(v3)|(spread3D_32(v>>20)<<2);
  }
uint32_t coord2morton3D_32 (std::array<uint32_t,3> xyz)
  {
  uint32_t v2=(xyz[0]&0x3ff)|((xyz[1]&0x3ff)<<10);
  uint64_t v3=spread3D_64(v2);
  v3=(v3|(v3>>29))&0x1fffffff;
  return uint32_t(v3)|(spread3D_32(xyz[2]&0x3ff)<<2);
  }
uint32_t morton2block3D_32 (uint32_t v)
  {
  return compress3D_32(v)
      | (compress3D_32(v>>1)<<10)
      | (compress3D_32(v>>2)<<20);
  }
std::array<uint32_t,3> morton2coord3D_32 (uint32_t v)
  {
  return {compress3D_32(v),
          compress3D_32(v>>1),
          compress3D_32(v>>2)};
  }

uint64_t block2morton3D_64 (uint64_t v)
  {
  return spread3D_64(v)
      | (spread3D_64(v>>21)<<1)
      | (spread3D_64(v>>42)<<2);
  }
uint64_t coord2morton3D_64 (std::array<uint64_t,3> xyz)
  {
  return spread3D_64(xyz[0])
      | (spread3D_64(xyz[1])<<1)
      | (spread3D_64(xyz[2])<<2);
  }
uint64_t morton2block3D_64 (uint64_t v)
  {
  return compress3D_64(v)
      | (compress3D_64(v>>1)<<21)
      | (compress3D_64(v>>2)<<42);
  }
std::array<uint64_t,3> morton2coord3D_64 (uint64_t v)
  {
  return { compress3D_64(v),
           compress3D_64(v>>1),
           compress3D_64(v>>2)};
  }

#endif

namespace {

const uint8_t m2p3D[24][8]={
{144,119, 97,110, 43, 44, 98,109},
{ 36, 35,101,106,127,152,102,105},
{ 96, 27,135, 28,161,162,182,181},
{170,169,189,190, 19,104, 20,143},
{134,137,151,112,133,138, 12, 11},
{130,141,  3,  4,129,142,120,159},
{174,173,185,186,103, 60,128, 59},
{ 52,111, 51,136,165,166,178,177},
{118,167,117, 92,153,168,154, 91},
{160,145, 83,146,175,126, 84,125},
{114, 75,113,176,157, 76,158,191},
{ 68,149,183,150, 67,122,184,121},
{ 16,139,  1,  2, 55,140, 14, 13},
{132, 63,  5,  6,131, 24, 10,  9},
{ 70,  7, 81, 32, 69,124, 82,123},
{116, 77,115, 90, 15, 78, 40, 89},
{ 38, 37, 23,108, 41, 42, 48,107},
{ 34, 33, 99, 56, 45, 46,100, 31},
{  0, 73, 39, 94,155, 74,156, 93},
{ 66,147, 85,148, 65,  8, 86, 47},
{ 72, 71,187,188, 17, 62, 18, 61},
{ 54, 25, 53, 26, 79, 64,180,179},
{172,171, 95, 80, 21, 58, 22, 57},
{ 50, 29, 49, 30,163,164, 88, 87}};

const uint8_t p2m3D[24][8]={
{144, 98,102, 44, 45,111,107,113},
{157,111,107, 33, 32, 98,102,124},
{ 96,164,165, 25, 27,183,182,130},
{109,169,168, 20, 22,186,187,143},
{115,137,141, 15, 14,132,128,146},
{126,132,128,  2,  3,137,141,159},
{134,186,187, 63, 61,169,168,100},
{139,183,182, 50, 48,164,165,105},
{173,156,158, 95, 91,114,112,161},
{160,145,147, 82, 86,127,125,172},
{179,114,112, 73, 77,156,158,191},
{190,127,125, 68, 64,145,147,178},
{ 16,  2,  3,137,141, 15, 14, 52},
{ 29, 15, 14,132,128,  2,  3, 57},
{ 35, 82, 86,127,125, 68, 64,  1},
{ 46, 95, 91,114,112, 73, 77, 12},
{ 54, 44, 45,111,107, 33, 32, 18},
{ 59, 33, 32, 98,102, 44, 45, 31},
{  0, 73, 77,156,158, 95, 91, 34},
{ 13, 68, 64,145,147, 82, 86, 47},
{ 72, 20, 22,186,187, 63, 61, 65},
{ 69, 25, 27,183,182, 50, 48, 76},
{ 83, 63, 61,169,168, 20, 22, 90},
{ 94, 50, 48,164,165, 25, 27, 87}};

const uint8_t m2p2D_1[4][4] = {
{ 4, 1, 11, 2},{0,15, 5, 6},{10,9,3,12},{14,7,13,8}};
uint8_t m2p2D_3[4][64];
const uint8_t p2m2D_1[4][4] = {
{ 4, 1, 3, 10},{0,6,7,13},{15,9,8,2},{11,14,12,5}};
uint8_t p2m2D_3[4][64];
bool peano2d_done=false;

void init_peano2d (void)
  {
  peano2d_done=true;

  for (unsigned d=0; d<4; ++d)
    for (uint32_t p=0; p<64; ++p)
      {
      unsigned rot = d;
      uint32_t v = p<<26;
      uint32_t res = 0;
      for (unsigned i=0; i<3; ++i)
        {
        unsigned tab=m2p2D_1[rot][v>>30];
        v<<=2;
        res = (res<<2) | (tab&0x3);
        rot = tab>>2;
        }
      m2p2D_3[d][p]=res|(rot<<6);
      }
  for (unsigned d=0; d<4; ++d)
    for (uint32_t p=0; p<64; ++p)
      {
      unsigned rot = d;
      uint32_t v = p<<26;
      uint32_t res = 0;
      for (unsigned i=0; i<3; ++i)
        {
        unsigned tab=p2m2D_1[rot][v>>30];
        v<<=2;
        res = (res<<2) | (tab&0x3);
        rot = tab>>2;
        }
      p2m2D_3[d][p]=res|(rot<<6);
      }
  }

}

uint32_t morton2peano3D_32(uint32_t v, unsigned bits)
  {
  unsigned rot = 0;
  uint32_t res = 0;
  v<<=3*(10-bits)+2;
  for (unsigned i=0; i<bits; ++i)
    {
    unsigned tab=m2p3D[rot][v>>29];
    v<<=3;
    res = (res<<3) | (tab&0x7);
    rot = tab>>3;
    }
  return res;
  }
uint32_t peano2morton3D_32(uint32_t v, unsigned bits)
  {
  unsigned rot = 0;
  uint32_t res = 0;
  v<<=3*(10-bits)+2;
  for (unsigned i=0; i<bits; ++i)
    {
    unsigned tab=p2m3D[rot][v>>29];
    v<<=3;
    res = (res<<3) | (tab&0x7);
    rot = tab>>3;
    }
  return res;
  }

uint64_t morton2peano3D_64(uint64_t v, unsigned bits)
  {
  unsigned rot = 0;
  uint64_t res = 0;
  v<<=3*(21-bits)+1;
  for (unsigned i=0; i<bits; ++i)
    {
    unsigned tab=m2p3D[rot][v>>61];
    v<<=3;
    res = (res<<3) | (tab&0x7);
    rot = tab>>3;
    }
  return res;
  }
uint64_t peano2morton3D_64(uint64_t v, unsigned bits)
  {
  unsigned rot = 0;
  uint64_t res = 0;
  v<<=3*(21-bits)+1;
  for (unsigned i=0; i<bits; ++i)
    {
    unsigned tab=p2m3D[rot][v>>61];
    v<<=3;
    res = (res<<3) | (tab&0x7);
    rot = tab>>3;
    }
  return res;
  }

uint32_t morton2peano2D_32(uint32_t v, unsigned bits)
  {
  if (!peano2d_done) init_peano2d();
  unsigned rot = 0;
  uint32_t res = 0;
  v<<=32-(bits<<1);
  unsigned i=0;
  for (; i+2<bits; i+=3)
    {
    unsigned tab=m2p2D_3[rot][v>>26];
    v<<=6;
    res = (res<<6) | (tab&0x3fu);
    rot = tab>>6;
    }
  for (; i<bits; ++i)
    {
    unsigned tab=m2p2D_1[rot][v>>30];
    v<<=2;
    res = (res<<2) | (tab&0x3);
    rot = tab>>2;
    }
  return res;
  }
uint32_t peano2morton2D_32(uint32_t v, unsigned bits)
  {
  if (!peano2d_done) init_peano2d();
  unsigned rot = 0;
  uint32_t res = 0;
  v<<=32-(bits<<1);
  unsigned i=0;
  for (; i+2<bits; i+=3)
    {
    unsigned tab=p2m2D_3[rot][v>>26];
    v<<=6;
    res = (res<<6) | (tab&0x3fu);
    rot = tab>>6;
    }
  for (; i<bits; ++i)
    {
    unsigned tab=p2m2D_1[rot][v>>30];
    v<<=2;
    res = (res<<2) | (tab&0x3);
    rot = tab>>2;
    }
  return res;
  }
uint64_t morton2peano2D_64(uint64_t v, unsigned bits)
  {
  if (!peano2d_done) init_peano2d();
  unsigned rot = 0;
  uint64_t res = 0;
  v<<=64-(bits<<1);
  unsigned i=0;
  for (; i+2<bits; i+=3)
    {
    unsigned tab=m2p2D_3[rot][v>>58];
    v<<=6;
    res = (res<<6) | (tab&0x3fu);
    rot = tab>>6;
    }
  for (; i<bits; ++i)
    {
    unsigned tab=m2p2D_1[rot][v>>62];
    v<<=2;
    res = (res<<2) | (tab&0x3);
    rot = tab>>2;
    }
  return res;
  }
uint64_t peano2morton2D_64(uint64_t v, unsigned bits)
  {
  if (!peano2d_done) init_peano2d();
  unsigned rot = 0;
  uint64_t res = 0;
  v<<=64-(bits<<1);
  unsigned i=0;
  for (; i+2<bits; i+=3)
    {
    unsigned tab=p2m2D_3[rot][v>>58];
    v<<=6;
    res = (res<<6) | (tab&0x3fu);
    rot = tab>>6;
    }
  for (; i<bits; ++i)
    {
    unsigned tab=p2m2D_1[rot][v>>62];
    v<<=2;
    res = (res<<2) | (tab&0x3);
    rot = tab>>2;
    }
  return res;
  }

}

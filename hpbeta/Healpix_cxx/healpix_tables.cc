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
 *  Copyright (C) 2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_tables.h"
#include "string_utils.h"

using namespace std;

const nside_dummy SET_NSIDE=nside_dummy();

Healpix_Ordering_Scheme string2HealpixScheme (const string &inp)
  {
  string tmp=trim(inp);
  if (equal_nocase(tmp,"RING")) return RING;
  if (equal_nocase(tmp,"NESTED")) return NEST;
  planck_fail ("bad Healpix ordering scheme '"+tmp+
               "': expected 'RING' or 'NESTED'");
  }

const uint16 Healpix_Tables::utab[] = {
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
X(0),X(1),X(4),X(5)
#undef X
#undef Y
#undef Z
};

const uint16 Healpix_Tables::ctab[] = {
#define Z(a) a,a+1,a+256,a+257
#define Y(a) Z(a),Z(a+2),Z(a+512),Z(a+514)
#define X(a) Y(a),Y(a+4),Y(a+1024),Y(a+1028)
X(0),X(8),X(2048),X(2056)
#undef X
#undef Y
#undef Z
};

const int Healpix_Tables::jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 },
          Healpix_Tables::jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

const uint8 Healpix_Tables::peano_subpix[2][8][4] =
  { { {0,1,3,2}, {3,0,2,1}, {2,3,1,0}, {1,2,0,3},
      {0,3,1,2}, {1,0,2,3}, {2,1,3,0}, {3,2,0,1} },
    { {0,1,3,2}, {1,3,2,0}, {3,2,0,1}, {2,0,1,3},
      {0,2,3,1}, {1,0,2,3}, {3,1,0,2}, {2,3,1,0} } };
const uint8 Healpix_Tables::peano_subpath[2][8][4] =
  { { {4,0,6,0}, {7,5,1,1}, {2,4,2,6}, {3,3,7,5},
      {0,2,4,4}, {5,1,5,3}, {6,6,0,2}, {1,7,3,7} },
    { {4,0,0,6}, {5,1,1,7}, {6,2,2,4}, {7,3,3,5},
      {0,4,4,2}, {1,5,5,3}, {2,6,6,0}, {3,7,7,1} } };
const uint8 Healpix_Tables::peano_face2path[2][12] =
  { { 2,5,2,5,3,6,3,6,2,3,2,3 }, { 2,6,2,3,3,5,2,6,2,3,3,5 } };
const uint8 Healpix_Tables::peano_face2face[2][12] =
  { { 0,5,6,11,10,1,4,7,2,3,8,9 }, { 0,5,8,9,6,1,2,7,10,11,4,3 } };

const int Healpix_Tables::nb_xoffset[] = { -1,-1, 0, 1, 1, 1, 0,-1 },
          Healpix_Tables::nb_yoffset[] = {  0, 1, 1, 1, 0,-1,-1,-1 };
const int Healpix_Tables::nb_facearray[][12] =
  { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },   // S
    {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },   // SE
    { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },   // E
    {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },   // SW
    {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },   // center
    {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },   // NE
    { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },   // W
    {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },   // NW
    {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } }; // N
const int Healpix_Tables::nb_swaparray[][3] =
  { { 0,0,3 },   // S
    { 0,0,6 },   // SE
    { 0,0,0 },   // E
    { 0,0,5 },   // SW
    { 0,0,0 },   // center
    { 5,0,0 },   // NE
    { 0,0,0 },   // W
    { 6,0,0 },   // NW
    { 3,0,0 } }; // N

const int Healpix_Tables::swap_clen[] =
  { 0,7,5,4,12,10,13,18,14,19,18,17,27,21 };
const int Healpix_Tables::swap_cycle[] =
  { 0,1,8,12,16,21,40,
    0,1,2,40,114,
    0,4,160,263,
    0,4,30,49,51,87,526,1027,1105,1387,1807,2637,
    0,8,10,18,39,74,146,307,452,4737,
    0,1,2,7,9,17,80,410,1526,1921,32859,33566,38931,
    0,5,6,10,12,24,27,95,372,494,924,1409,3492,4248,9137,66043,103369,156899,
    0,1,2,3,4,45,125,351,697,24337,102940,266194,341855,419857,
    0,1,2,3,9,16,1705,2082,2126,8177,12753,15410,52642,80493,83235,88387,99444,
      1675361,2495125,
    0,2,6,8,9,11,20,50,93,152,183,2137,13671,44794,486954,741908,4803258,
      5692573,
    0,1,5,6,44,53,470,2847,3433,4906,13654,14710,400447,1797382,2744492,
      18775974,23541521,
    0,4,9,10,16,33,83,117,318,451,5759,10015,128975,171834,211256,347608,
      1278690,2154097,2590798,3427694,5581717,21012301,27023976,72522811,
      95032729,139166747,171822389,
    0,5,10,267,344,363,2968,3159,9083,18437,76602,147614,1246902,1593138,
      2035574,6529391,9511830,11340287,29565945,281666026,677946848 };

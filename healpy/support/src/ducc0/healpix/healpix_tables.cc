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
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2011-2020 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ducc0/healpix/healpix_tables.h"
#include "ducc0/infra/string_utils.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_healpix {

using namespace std;

const nside_dummy SET_NSIDE=nside_dummy();

Ordering_Scheme string2HealpixScheme (const string &inp)
  {
  string tmp=trim(inp);
  if (equal_nocase(tmp,"RING")) return RING;
  if (equal_nocase(tmp,"NESTED")) return NEST;
  MR_fail ("bad Healpix ordering scheme '"+tmp+
           "': expected 'RING' or 'NESTED'");
  }

const int Healpix_Tables::jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 },
          Healpix_Tables::jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

const uint8_t Healpix_Tables::peano_arr2[] = {
        0, 35, 65, 66, 68,  5,103,  6,110,109, 15, 44, 72,  9,107, 10,
       31,126, 60,125, 81, 16, 82, 51,123, 88, 26, 25,119, 84, 22, 21,
       42, 75, 41,104, 12, 47, 77, 78, 38, 71, 37,100, 98, 97,  3, 32,
       53, 54,116, 87, 57, 58,120, 91, 19,114, 48,113, 93, 28, 94, 63,
       64,  1, 99,  2, 46, 79, 45,108,  4, 39, 69, 70,  8, 43, 73, 74,
       85, 20, 86, 55,115, 80, 18, 17, 89, 24, 90, 59, 61, 62,124, 95,
      106,105, 11, 40,102,101,  7, 36, 76, 13,111, 14, 34, 67, 33, 96,
      127, 92, 30, 29, 27,122, 56,121, 49, 50,112, 83, 23,118, 52,117,

      128,194,195,161,196,133,135,230,204,141,143,238,171,233,232,138,
      149,212,214,183,221,159,158,252,217,155,154,248,178,243,241,144,
      175,237,236,142,235,170,168,201,227,162,160,193,132,198,199,165,
      186,251,249,152,242,176,177,211,246,180,181,215,157,220,222,191,
      192,129,131,226,136,202,203,169,140,206,207,173,231,166,164,197,
      213,151,150,244,145,208,210,179,153,216,218,187,254,188,189,223,
      239,174,172,205,167,229,228,134,163,225,224,130,200,137,139,234,
      250,184,185,219,190,255,253,156,182,247,245,148,209,147,146,240 };
const uint8_t Healpix_Tables::peano_arr[] =
      { 16, 1,27, 2,31,20, 6, 5,10,19, 9,24,13,14,28,23,
         0,11,17,18,21, 4,22,15,26,25, 3, 8, 7,30,12,29,
        48,33,35,58,53,39,38,60,59,42,40,49,62,44,45,55,
        32,50,51,41,37,52,54,47,43,57,56,34,46,63,61,36 };
const uint8_t Healpix_Tables::peano_face2path[2][12] =
  { { 2,5,2,5,3,6,3,6,2,3,2,3 }, { 2,6,2,3,3,5,2,6,2,3,3,5 } };
const uint8_t Healpix_Tables::peano_face2face[2][12] =
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

const size_t Healpix_Tables::swap_clen[] =
  { 0,7,5,4,12,10,13,18,14,19,18,17,27,21 };
const size_t Healpix_Tables::swap_cycle[] =
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

}}

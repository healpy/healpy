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

/*! \file healpix_tables.h
 * 
 *  \copyright Copyright (C) 2011-2020 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_TABLES_H
#define HEALPIX_TABLES_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace ducc0 {

namespace detail_healpix {

/*! The two possible ordering schemes of a HEALPix map. */
enum Ordering_Scheme { RING, /*!< RING scheme */
                       NEST  /*!< NESTED scheme */
                     };

Ordering_Scheme string2HealpixScheme (const std::string &inp);

class nside_dummy {};
extern const nside_dummy SET_NSIDE;

class Healpix_Tables
  {
  protected:
    static const int jrll[], jpll[];

    static const uint8_t peano_face2path[2][12], peano_face2face[2][12],
                         peano_arr[],peano_arr2[];

    static const int nb_xoffset[], nb_yoffset[],
                     nb_facearray[][12], nb_swaparray[][3];

    static const size_t swap_clen[], swap_cycle[];
  };

}

using detail_healpix::RING;
using detail_healpix::NEST;
using detail_healpix::SET_NSIDE;

}

#endif

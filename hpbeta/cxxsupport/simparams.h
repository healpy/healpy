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
 *  Class for storing parameter information for later use
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_SIMPARAMS_H
#define PLANCK_SIMPARAMS_H

#include <string>
#include <vector>
#include <iostream>
#include "cxxutils.h"
class fitshandle;

class simparams
  {
  private:
    class Param
      {
      public:
        std::string key, shortkey, value, comment;

        Param (const std::string &Key, const std::string &Shortkey,
               const std::string &Value, const std::string &Comment)
          : key(Key), shortkey(Shortkey), value(Value), comment(Comment) {}
      };

    std::vector<Param> paramMap;
    std::vector<std::string> source_files;
    std::vector<int> hdus;

  public:
    void add_comment (const std::string &comment)
      { paramMap.push_back(Param("","","",comment)); }
    template<typename T> void add(const std::string &key,
      const std::string &shortkey, const T &value, const std::string &comment)
      {
      paramMap.push_back(Param(key, shortkey, dataToString(value), comment));
      }
    template<typename T> void add(const std::string &key,
      const std::string &shortkey, const T &value)
      {
      paramMap.push_back(Param(key, shortkey, dataToString(value), ""));
      }

    void add_source_file (const std::string &filename, int hdu=2)
      {
      source_files.push_back(filename);
      hdus.push_back(hdu);
      }

    void add_keys (std::ostream &os) const;
    void add_keys (fitshandle &out) const;
  };

#endif

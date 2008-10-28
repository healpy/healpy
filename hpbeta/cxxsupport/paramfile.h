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
 *  Class for parsing parameter files
 *
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <set>
#include <string>
#include <iostream>
#include "simparams.h"
#include "cxxutils.h"

class paramfile
  {
  private:
    typedef std::map<std::string,std::string> params_type;
    params_type params;
    mutable std::set<std::string> read_params;
    bool verbose;

    std::string get_valstr(const std::string &key) const
      {
      params_type::const_iterator loc=params.find(key);
      if (loc!=params.end()) return loc->second;
      throw Message_error ("Error: Cannot find the key \"" + key + "\".");
      }

  public:
    paramfile (const std::string &filename, bool verbose_=true)
      : verbose(verbose_)
      { parse_file (filename, params); }

    paramfile (const params_type &par)
      : params (par), verbose(true)
      {}

    ~paramfile ()
      {
      if (verbose)
        for (params_type::const_iterator loc=params.begin();
             loc!=params.end(); ++loc)
          if (read_params.find(loc->first)==read_params.end())
            std::cout << "Parser warning: unused parameter "
                      << loc->first << std::endl;
      }

    bool param_present(const std::string &key) const
      { return (params.find(key)!=params.end()); }

    template<typename T> T find (const std::string &key) const
      {
      T result;
      stringToData(get_valstr(key),result);
      if (verbose)
        std::cout << "Parser: " << key << " = " << dataToString(result)
                  << std::endl;
      read_params.insert(key);
      return result;
      }
    template<typename T> T find
      (const std::string &key, const T &deflt)
      {
      if (param_present(key)) return find<T>(key);
      if (verbose)
        std::cout << "Parser: " << key << " = " << dataToString(deflt)
                  << " <default>" << std::endl;
      params[key]=dataToString(deflt);
      read_params.insert(key);
      return deflt;
      }

    const params_type &getParams() const
      { return params; }

    template<typename T> void findParam
      (const std::string &key, T &value) const
      { value = find<T>(key); }

    template<typename T> void findHeaderParam(const std::string& key,
      T& value, simparams& headerParams, const std::string& headerKey,
      const std::string& headerComment) const
      {
      findParam(key, value);
      headerParams.add(key, headerKey, dataToString(value), headerComment);
      }
    void findSourceParam(const std::string& key, std::string& value,
      simparams& headerParams) const
      {
      findParam(key, value);
      headerParams.add_source_file(value);
      }
  };

#endif

/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Class for parsing parameter files
 *
 *  Copyright (C) 2003, 2004, 2005, 2008, 2009, 2010 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <set>
#include <string>
#include <iostream>
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
      planck_fail ("Cannot find the key '" + key + "'.");
      }

    bool param_unread (const std::string &key) const
      { return (read_params.find(key)==read_params.end()); }

  public:
    paramfile (const std::string &filename, bool verbose_=true)
      : verbose(verbose_)
      { parse_file (filename, params); }

    paramfile (const params_type &par, bool verbose_=true)
      : params (par), verbose(verbose_)
      {}

    ~paramfile ()
      {
      if (verbose)
        for (params_type::const_iterator loc=params.begin();
             loc!=params.end(); ++loc)
          if (param_unread(loc->first))
            std::cout << "Parser warning: unused parameter '"
                      << loc->first << "'" << std::endl;
      }

    void setVerbosity (bool verbose_)
      { verbose = verbose_; }

    bool getVerbosity () const
      { return verbose; }

    bool param_present(const std::string &key) const
      { return (params.find(key)!=params.end()); }

    template<typename T> T find (const std::string &key) const
      {
      T result;
      stringToData(get_valstr(key),result);
      std::string output = dataToString(result);
      if (planckType<T>()==PLANCK_STRING) output = "'"+output+"'";
      if (verbose && param_unread(key))
        std::cout << "Parser: " << key << " = " << output << std::endl;
      read_params.insert(key);
      return result;
      }
    template<typename T> T find
      (const std::string &key, const T &deflt)
      {
      if (param_present(key)) return find<T>(key);
      std::string output = dataToString(deflt);
      if (planckType<T>()==PLANCK_STRING) output = "'"+output+"'";
      if (verbose && param_unread(key))
        std::cout << "Parser: " << key << " = " << output
                  << " <default>" << std::endl;
      params[key]=dataToString(deflt);
      read_params.insert(key);
      return deflt;
      }

    const params_type &getParams() const
      { return params; }
  };

#endif

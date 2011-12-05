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
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <set>
#include <string>
#include "datatypes.h"
#include "string_utils.h"

class paramfile
  {
  private:
    typedef std::map<std::string,std::string> params_type;
    params_type params;
    mutable std::set<std::string> read_params;
    bool verbose;

    std::string get_valstr(const std::string &key) const;
    bool param_unread (const std::string &key) const;
    void findhelper (const std::string &key, const std::string &value, NDT type,
      bool deflt) const;

  public:
    paramfile (const std::string &filename, bool verbose_=true);
    paramfile (const params_type &par, bool verbose_=true);
    ~paramfile();

    void setVerbosity (bool verbose_)
      { verbose = verbose_; }

    bool getVerbosity () const
      { return verbose; }

    bool param_present(const std::string &key) const;

    template<typename T> T find (const std::string &key) const;
    template<typename T> T find
      (const std::string &key, const T &deflt);

    const params_type &getParams() const
      { return params; }

    void setParamString (const std::string &key, const std::string &value);
    template<typename T> void setParam (const std::string &key, const T &value)
      { setParamString(key,dataToString(value)); }
  };

#endif

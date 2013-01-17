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

#include <iostream>
#include "paramfile.h"
#include "string_utils.h"

using namespace std;

string paramfile::get_valstr(const string &key) const
  {
  params_type::const_iterator loc=params.find(key);
  if (loc!=params.end()) return loc->second;
  planck_fail ("Cannot find the key '" + key + "'.");
  }

bool paramfile::param_unread (const string &key) const
  { return (read_params.find(key)==read_params.end()); }

paramfile::paramfile (const string &filename, bool verbose_)
  : verbose(verbose_)
  { parse_file (filename, params); }

paramfile::paramfile (const params_type &par, bool verbose_)
  : params (par), verbose(verbose_)
  {}

paramfile::~paramfile()
  {
  if (verbose)
    for (params_type::const_iterator loc=params.begin();
      loc!=params.end(); ++loc)
      if (param_unread(loc->first))
      cout << "Parser warning: unused parameter '"
           << loc->first << "'" << endl;
  }

bool paramfile::param_present(const string &key) const
  { return (params.find(key)!=params.end()); }

void paramfile::findhelper (const string &key, const string &value, NDT type,
  bool deflt) const
  {
  string output = value;
  if (type==NAT_STRING) output = "'"+output+"'";
  if (verbose && param_unread(key))
    cout << "Parser: " << key << " = " << output
         << (deflt ? " <default>" : "") << endl;
  read_params.insert(key);
  }

template<typename T> T paramfile::find (const string &key) const
  {
  T result = stringToData<T>(get_valstr(key));
  findhelper (key, dataToString(result), nativeType<T>(), false);
  return result;
  }

template unsigned char paramfile::find (const string &key) const;
template signed char paramfile::find (const string &key) const;
template unsigned short paramfile::find (const string &key) const;
template short paramfile::find (const string &key) const;
template unsigned int paramfile::find (const string &key) const;
template int paramfile::find (const string &key) const;
template unsigned long paramfile::find (const string &key) const;
template long paramfile::find (const string &key) const;
template unsigned long long paramfile::find (const string &key) const;
template long long paramfile::find (const string &key) const;
template float paramfile::find (const string &key) const;
template double paramfile::find (const string &key) const;
template long double paramfile::find (const string &key) const;
template bool paramfile::find (const string &key) const;
template string paramfile::find (const string &key) const;

template<typename T> T paramfile::find (const string &key, const T &deflt)
  {
  if (param_present(key)) return find<T>(key);
  string sdeflt=dataToString(deflt);
  findhelper (key, sdeflt, nativeType<T>(), true);
  params[key]=sdeflt;
  return deflt;
  }

template unsigned char paramfile::find (const string &key,
  const unsigned char &deflt);
template signed char paramfile::find (const string &key,
  const signed char &deflt);
template unsigned short paramfile::find (const string &key,
  const unsigned short &deflt);
template short paramfile::find (const string &key, const short &deflt);
template unsigned int paramfile::find (const string &key,
  const unsigned int &deflt);
template int paramfile::find (const string &key, const int &deflt);
template unsigned long paramfile::find (const string &key,
  const unsigned long &deflt);
template long paramfile::find (const string &key, const long &deflt);
template unsigned long long paramfile::find (const string &key,
  const unsigned long long &deflt);
template long long paramfile::find (const string &key, const long long &deflt);
template float paramfile::find (const string &key, const float &deflt);
template double paramfile::find (const string &key, const double &deflt);
template long double paramfile::find (const string &key,
  const long double &deflt);
template bool paramfile::find (const string &key, const bool &deflt);
template string paramfile::find (const string &key, const string &deflt);

void paramfile::setParamString (const string &key, const string &value)
  {
  if (param_present(key))
    {
    if (params[key]!=value)
      {
      if (verbose)
        cout << "Parser: altering value of key'"<<key<<"' to '"<<value<<"'."
             << endl;
      params[key]=value;
      }
    }
  else
    {
    if (verbose)
      cout << "Parser: setting new key'"<<key<<"' to '"<<value<<"'."<<endl;
    params[key]=value;
    }
  }

paramfile getParamsFromCmdline (int argc, const char **argv, bool verbose)
  {
  planck_assert(argc>=2,"incorrect command line format");
  if ((argc==2) && (string(argv[1]).find("=")==string::npos))
    return paramfile(argv[1],verbose);
  map<string,string> pmap;
  parse_cmdline_equalsign(argc,argv,pmap);
  return paramfile(pmap,verbose);
  }

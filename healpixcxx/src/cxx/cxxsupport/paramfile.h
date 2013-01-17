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

/*! \file paramfile.h
 *  Class for parsing parameter files
 *
 *  Copyright (C) 2003-2012 Max-Planck-Society
 *  Authors: Martin Reinecke
 */

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <set>
#include <string>
#include "datatypes.h"
#include "string_utils.h"

/*! Class for storing and querying key/value pairs. The name is historical;
    the parameters can actually be obtained from othersources as well
    (e.g. the command line). */
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
    void setParamString (const std::string &key, const std::string &value);

  public:
    paramfile () {}
    /*! Constructs a paramfile object from the contents of \a filename.
        If  \a verbose_==true, diagnostic output is generated when calling
        methods on this object, otherwise not. */
    paramfile (const std::string &filename, bool verbose_=true);
    /*! Constructs a paramfile object from the contents of \a par.
        If  \a verbose_==true, diagnostic output is generated when calling
        methods on this object, otherwise not. */
    paramfile (const params_type &par, bool verbose_=true);
    ~paramfile();

    /*! Allows adjusting the verbosity. */
    void setVerbosity (bool verbose_)
      { verbose = verbose_; }

    /*! Returns the verbosity setting of the object. */
    bool getVerbosity () const
      { return verbose; }

    /*! Returns \c true, if a paremeter called \a key is stored in the object,
        else \c false. */
    bool param_present(const std::string &key) const;

    /*! Returns the value stored for the parameter name \a key, after converting
        it to the requested type. If \a key is not present, an exception is
        thrown. */
    template<typename T> T find (const std::string &key) const;
    /*! Returns the value stored for the parameter name \a key, after converting
        it to the requested type. If \a key is not present, \a deflt is returned
        instead, and is also entered into the parameter set. */
    template<typename T> T find
      (const std::string &key, const T &deflt);

    /*! Returns the entire set of currently stored parameters. */
    const params_type &getParams() const
      { return params; }

    /*! Sets the parameter with the name \a key to \a value. */
    template<typename T> void setParam (const std::string &key, const T &value)
      { setParamString(key,dataToString(value)); }
  };

/*! Tries to build a \a paramfile object from the contents of a command line.
    If \a argc==2 and \a argv[1] does not contain the character "=", the
    function tries to input parameters from the file \a argv[1]. Otherwise
    the function interprets each command line argument as a "key=value"
    statement. */
paramfile getParamsFromCmdline (int argc, const char **argv,
  bool verbose=true);

#endif

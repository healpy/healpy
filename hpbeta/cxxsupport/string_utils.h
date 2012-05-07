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

/*! \file string_utils.h
 *  Various functions for manipulating strings.
 *
 *  Copyright (C) 2002-2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_STRING_UTILS_H
#define PLANCK_STRING_UTILS_H

#include <vector>
#include <map>
#include "datatypes.h"

/*! \defgroup stringutilsgroup String handling helper functions */
/*! \{ */

/*! Returns the string \a orig without leading and trailing whitespace. */
std::string trim (const std::string &orig);

/*! Returns a string containing the text representation of \a x.
    Care is taken that no information is lost in the conversion. */
template<typename T> std::string dataToString(const T &x);
template<> std::string dataToString (const bool &x);
template<> std::string dataToString (const std::string &x);
template<> std::string dataToString (const float &x);
template<> std::string dataToString (const double &x);
template<> std::string dataToString (const long double &x);

/*! Returns a string containing the text representation of \a x, padded
    with leading zeroes to \a width characters. */
std::string intToString(int64 x, tsize width);

/*! Reads a value of a given datatype from a string */
template<typename T> void stringToData (const std::string &x, T &value);
template<> void stringToData (const std::string &x, std::string &value);
template<> void stringToData (const std::string &x, bool &value);

/*! Reads a value of a given datatype from a string */
template<typename T> inline T stringToData (const std::string &x)
  { T result; stringToData(x,result); return result; }

/*! Parses the file \a filename and returns the key/value pairs in \a dict. */
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

void parse_cmdline_classic (int argc, const char **argv,
  const std::vector<std::string> &leading_args,
  std::map<std::string,std::string> &dict);

void parse_cmdline_classic (int argc, const char **argv,
  std::map<std::string,std::string> &dict);

void parse_cmdline_equalsign (int argc, const char **argv,
  const std::vector<std::string> &leading_args,
  std::map<std::string,std::string> &dict);

void parse_cmdline_equalsign (int argc, const char **argv,
  std::map<std::string,std::string> &dict);

/*! Case-insensitive string comparison
    Returns \a true, if \a a and \a b differ only in capitalisation,
    else \a false. */
bool equal_nocase (const std::string &a, const std::string &b);

/*! Returns lowercase version of \a input. */
std::string tolower(const std::string &input);

/*! Tries to split \a inp into a white-space separated list of values of
    type \a T, and appends them to \a list. */
template<typename T> void split (const std::string &inp, std::vector<T> &list);

/*! Breaks the string \a inp into tokens separated by \a delim, and returns them
    in \a list. */
void tokenize (const std::string &inp, char delim,
  std::vector<std::string> &list);

/*! Reads all white-space separated strings from \a filename, and returns
    them in \a words. */
void parse_words_from_file (const std::string &filename,
  std::vector<std::string> &words);

/*! \} */

#endif

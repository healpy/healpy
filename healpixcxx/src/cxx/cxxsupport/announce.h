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

/*! \file announce.h
 *  Functions for printing information at startup.
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ANNOUNCE_H
#define PLANCK_ANNOUNCE_H

#include <string>

/*! Prints a banner containing \a name, as well as some information about the
    source code and the parallelisation techniques enabled. */
void announce (const std::string &name);

/*! Prints a banner containing \a name and checks if \a argc_valid is true.
    If not, the string \a usage is printed and the program is terminated. */
void module_startup (const std::string &name, bool argc_valid,
  const std::string &usage, bool verbose=true);

/*! Prints a banner containing \a name and checks if \a argc==argc_expected.
    If not, a usage description is given and the program is terminated. */
void module_startup (const std::string &name, int argc, const char **argv,
  int argc_expected, const std::string &argv_expected, bool verbose=true);

/*! Prints a banner containing \a name and checks if \a argc>=2.
    If not, a usage description is given and the program is terminated. */
void module_startup (const std::string &name, int argc, const char **argv,
  bool verbose=true);

#endif

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
 *  This file contains the implementation of various convenience functions
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

// if we are using g++, check for version 4.0 or higher
#ifdef __GNUC__
#if (__GNUC__<4)
#error your C++ compiler is too old. g++ version 4.0 or higher is required.
#endif
#endif

#include <iostream>
#include "announce.h"
#include "openmp_support.h"

using namespace std;

namespace {

void openmp_status()
  {
#ifndef _OPENMP
  cout << "OpenMP: not supported by this binary" << endl;
#else
  int threads = openmp_max_threads();
  if (threads>1)
    cout << "OpenMP active: max. " << threads << " threads." << endl;
  else
    cout << "OpenMP active, but running with 1 thread only." << endl;
#endif
  }

void vec_status()
  {
  cout << "Vector math: ";
#if(defined(__AVX__))
  cout << "AVX" << endl;
#elif(defined(__SSE2__))
  cout << "SSE2" << endl;
#elif(defined(__SSE__))
  cout << "SSE" << endl;
#else
  cout << "not supported by this binary" << endl;
#endif
  }

} //unnamed namespace

void announce (const string &name)
  {
  string version = "v3.1 (experimental)";
  string name2 = name+" "+version;
  cout << endl << "+-";
  for (tsize m=0; m<name2.length(); ++m) cout << "-";
  cout << "-+" << endl;
  cout << "| " << name2 << " |" << endl;
  cout << "+-";
  for (tsize m=0; m<name2.length(); ++m) cout << "-";
  cout << "-+" << endl << endl;
  vec_status();
  openmp_status();
  cout << endl;
  }

void module_startup (const string &name, bool argc_valid, const string &usage,
  bool verbose)
  {
  if (verbose) announce (name);
  if (argc_valid) return;
  if (verbose) cerr << usage << endl;
  planck_fail_quietly ("Incorrect usage");
  }

void module_startup (const string &name, int argc, const char **,
  int argc_expected, const string &argv_expected, bool verbose)
  {
  module_startup (name,argc==argc_expected,
    string("Usage: ")+name+" "+argv_expected, verbose);
  }

void module_startup (const std::string &name, int argc, const char ** /*argv*/,
  bool verbose)
  {
  module_startup (name,argc>=2,
    string("Usage:\n  ")+name+" <parameter file / init object>\nor:\n  "
                        +name+" par1=val1 par2=val2 ...", verbose);
  }

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
 *  Utilities for error reporting
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "error_handling.h"

using namespace std;

PlanckError::PlanckError(const string &message) : msg (message) {}
PlanckError::PlanckError(const char *message) : msg (message) {}

//virtual
PlanckError::~PlanckError() {}

void planck_failure__(const char *file, int line, const char *func,
  const string &msg)
  {
  cerr << "Error encountered at " << file << ", line " << line << endl;
  if (func) cerr << "(function " << func << ")" << endl;
  if (msg!="") cerr << endl << msg << endl;
  cerr << endl;
  }

void planck_failure__(const char *file, int line, const char *func,
  const char *msg)
  { planck_failure__ (file,line,func,string(msg)); }

void killjob__()
  { throw; }

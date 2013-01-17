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

/*! \file error_handling.h
 *  Utilities for error reporting
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_ERROR_HANDLING_H
#define PLANCK_ERROR_HANDLING_H

#include <string>
#include <iostream>

#if defined (__GNUC__)
#define PLANCK_FUNC_NAME__ __PRETTY_FUNCTION__
#else
#define PLANCK_FUNC_NAME__ 0
#endif

void planck_failure__(const char *file, int line, const char *func,
  const std::string &msg);
void planck_failure__(const char *file, int line, const char *func,
  const char *msg);
void killjob__();

class PlanckError
  {
  private:
    std::string msg;

  public:
    explicit PlanckError(const std::string &message);
    explicit PlanckError(const char *message);

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~PlanckError();
  };

/*! \defgroup errorgroup Error handling */
/*! \{ */

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg) \
do { planck_failure__(__FILE__,__LINE__,PLANCK_FUNC_NAME__,msg); \
throw PlanckError(msg); } while(0)

/*! Throws a PlanckError without diagnostic message. */
#define planck_fail_quietly(msg) \
do { throw PlanckError(msg); } while(0)

/*! Writes diagnostic output and exits with an error status if \a testval
    is \a false. */
#define planck_assert(testval,msg) \
do { if (testval); else planck_fail(msg); } while(0)

/*! Macro for improving error diagnostics. Should be placed immediately
    after the opening brace of \c main(). Must be used in conjunction with
    \c PLANCK_DIAGNOSIS_END. */
#define PLANCK_DIAGNOSIS_BEGIN try {
/*! Macro for improving error diagnostics. Should be placed immediately
    before the closing brace of \c main(). Must be used in conjunction with
    \c PLANCK_DIAGNOSIS_BEGIN. */
#define PLANCK_DIAGNOSIS_END \
} \
catch (PlanckError &) \
  { killjob__(); /* no need for further diagnostics; they were shown already */ } \
catch (std::exception &e) \
  { std::cerr << "std::exception: " << e.what() << std::endl; killjob__(); } \
catch (...) \
  { std::cerr << "Unknown exception" << std::endl; killjob__(); }

/*! \} */

#endif

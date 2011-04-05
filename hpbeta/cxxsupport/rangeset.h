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

/*! \file rangeset.h
 *  Class for storing sets of ranges of integer numbers
 *
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_RANGESET_H
#define PLANCK_RANGESET_H

#include <vector>
#include <utility>
#include "cxxutils.h"

/*! Class for storing sets of ranges of integer numbers */
template<typename T> class rangeset: public std::vector<std::pair<T,T> >
  {
  public:
    /*! Append the interval [\a vbegin; \a vend[.
        \note \a vend is not included in the interval. */
    void append(T vbegin, T vend)
      {
      if (vend<=vbegin) return;
      if (!(this->empty()) && (vbegin<=this->back().second))
        {
        planck_assert (vbegin==this->back().second,"bad append operation");
        this->back().second = vend;
        }
      else
        this->push_back(std::pair<T,T>(vbegin,vend));
      }
    /*! Append the one-number interval containing \a value. */
    void append(T value)
      { append (value, value+1); }

    T nval() const
      {
      T result=T(0);
      typename rangeset::const_iterator it;
      for (it=this->begin(); it!=this->end(); ++it)
        result+=it->second-it->first;
      return result;
      }

    void toVector (std::vector<T> &res)
      {
      res.clear();
      res.reserve(nval());
      typename rangeset::const_iterator it;
      for (it=this->begin(); it!=this->end(); ++it)
        for (T m=it->first; m<it->second; ++m)
          res.push_back(m);
      }
  };

#endif

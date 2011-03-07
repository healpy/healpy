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
 *  Class for storing sets of ranges of inteeger numbers
 *
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_RANGESET_H
#define PLANCK_RANGESET_H

#include <map>

/*! Class for storing sets of ranges of integer numbers */
template<typename T> class rangeset: public std::map<T,T>
  {
  public:
    /*! Add the interval [\a vbegin; \a vend[.
        \note \a vend is not included in the interval. */
    void add(T vbegin, T vend)
      {
      if (vend<=vbegin) return;
      rangeset<int>::iterator left=this->lower_bound(vbegin);
      if (left!=this->begin())
        {
        rangeset<int>::iterator preleft=left;
        --preleft;
        if (vbegin<=preleft->second) left=preleft;
        }
      rangeset<int>::iterator pastright=this->lower_bound(vend);
      if ((pastright!=this->end())&&(vend==pastright->first))++pastright;
      if (left==pastright) // between two intervals
        this->insert(std::pair<T,T>(vbegin,vend));
      else
        {
        rangeset<int>::iterator right = pastright;
        if (right!=this->begin()) --right;
        vbegin=std::min(vbegin,left->first);
        vend=std::max(vend,right->second);
        this->erase(left,pastright);
        this->insert(std::pair<T,T>(vbegin,vend));
        }
      }
    /*! Add the one-number interval containing \a value. */
    void add(T value)
      { add (value, value+1); }
  };

#endif

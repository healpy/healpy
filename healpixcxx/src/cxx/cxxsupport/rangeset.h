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
 *  Copyright (C) 2011, 2012 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_RANGESET_H
#define PLANCK_RANGESET_H

#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include "datatypes.h"
#include "error_handling.h"

/*! Class for storing sets of ranges of integer numbers */
template<typename T> class rangeset
  {
  private:
    typedef std::vector<T> rtype;
    typedef typename rtype::iterator iterator;
    typedef typename rtype::const_iterator c_iterator;
    rtype r;

    tdiff iiv (const T &val) const
      { return tdiff(std::upper_bound(r.begin(),r.end(),val)-r.begin())-1; }

    void addRemove (T a, T b, tdiff v)
      {
      tdiff pos1=iiv(a), pos2=iiv(b);
      if ((pos1>=0) && (r[pos1]==a)) --pos1;
      // first to delete is at pos1+1; last is at pos2
      bool insert_a = (pos1&1)==v;
      bool insert_b = (pos2&1)==v;
      int rmstart=pos1+1+(insert_a ? 1 : 0);
      int rmend  =pos2-(insert_b?1:0);

      planck_assert((rmend-rmstart)&1,"cannot happen");

      if (insert_a && insert_b && (pos1+1>pos2)) // insert
        {
        r.insert(r.begin()+pos1+1,2,a);
        r[pos1+2]=b;
        }
      else
        {
        if (insert_a) r[pos1+1]=a;
        if (insert_b) r[pos2]=b;
        r.erase(r.begin()+rmstart,r.begin()+rmend+1);
        }
      }

    static void generalUnion (const rtype &a, const rtype &b,
      bool flip_a, bool flip_b, rtype &c)
      {
      planck_assert((&c!=&a)&&(&c!=&b), "cannot overwrite the rangeset");
      c.clear();
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      tsize ia=0, ea=a.size(), ib=0, eb=b.size();
      bool runa = ia!=ea, runb = ib!=eb;
      while(runa||runb)
        {
        bool adv_a=false, adv_b=false;
        T val,va=T(),vb=T();
        if (runa) va = a[ia];
        if (runb) vb = b[ib];
        if (runa && (!runb || (va<=vb))) { adv_a=true; val=va; }
        if (runb && (!runa || (vb<=va))) { adv_b=true; val=vb; }
        if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
        if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
        if ((state_a||state_b)!=state_res)
          { c.push_back(val); state_res = !state_res; }
        }
      }

  public:
    /*! Removes all rangeset entries. */
    void clear() { r.clear(); }
    /*! Reserves space for \a n ranges. */
    void reserve(tsize n) { r.reserve(2*n); }
    /*! Returns the current number of ranges. */
    tsize size() const { return r.size()>>1; }
    /*! Returns the current vector of ranges. */
    const rtype &data() const { return r; }

    /*! Returns the first value of range \a i. */
    const T &ivbegin (tdiff i) const { return r[2*i]; }
    /*! Returns the one-past-last value of range \a i. */
    const T &ivend (tdiff i) const { return r[2*i+1]; }
    /*! Returns the length of range \a i. */
    T ivlen (tdiff i) const { return r[2*i+1]-r[2*i]; }

    /*! Appends \a [v1;v2[ to the rangeset. \a v1 must be larger
        than the minimum of the last range in the rangeset. */
    void append(const T &v1, const T &v2)
      {
      if (v2<=v1) return;
      if ((!r.empty()) && (v1<=r.back()))
        {
        planck_assert (v1>=r[r.size()-2],"bad append operation");
        if (v2>r.back()) r.back()=v2;
        }
      else
        { r.push_back(v1); r.push_back(v2); }
      }
    /*! Appends \a [v;v+1[ to the rangeset. \a v must be larger
        than the minimum of the last range in the rangeset. */
    void append(const T &v)
      { append(v,v+1); }

    /*! Appends \a other to the rangeset. All values in \a other must be larger
        than the minimum of the last range in the rangeset. */
    void append (const rangeset &other)
      {
      for (tsize j=0; j<other.size(); ++j)
        append(other.ivbegin(j),other.ivend(j));
      }

    /*! After this operation, the rangeset contains the union of itself
        with \a [v1;v2[. */
    void add(const T &v1, const T &v2) { addRemove(v1,v2,1); }
    /*! After this operation, the rangeset contains the union of itself
        with \a [v;v+1[. */
    void add(const T &v) { addRemove(v,v+1,1); }

    /*! Removes all values within \a [v1;v2[ from the rangeset. */
    void remove(const T &v1, const T &v2) { addRemove(v1,v2,0); }
    /*! Removes the value \a v from the rangeset. */
    void remove(const T &v) { addRemove(v,v+1,0); }

    /*! Removes all values not within \a [v1;v2[ from the rangeset. */
    void intersect (const T &a, const T &b)
      {
      tdiff pos1=iiv(a), pos2=iiv(b);
      if ((pos2>=0) && (r[pos2]==b)) --pos2;
      // delete all up to pos1 (inclusive); and starting from pos2+1
      bool insert_a = (pos1&1)==0;
      bool insert_b = (pos2&1)==0;

      // cut off end
      r.erase(r.begin()+pos2+1,r.end());
      if (insert_b) r.push_back(b);

      // erase start
      if (insert_a) r[pos1--]=a;
      if (pos1>=0)
        r.erase(r.begin(),r.begin()+pos1+1);
      }

    /*! Returns the total number of elements in the rangeset. */
    T nval() const
      {
      T result=T(0);
      for (tsize i=0; i<r.size(); i+=2)
        result+=r[i+1]-r[i];
      return result;
      }

    /*! After this opration, \a res contains all elements of the rangeset
        in ascending order. */
    void toVector (std::vector<T> &res) const
      {
      res.clear();
      res.reserve(nval());
      for (tsize i=0; i<r.size(); i+=2)
        for (T m(r[i]); m<r[i+1]; ++m)
          res.push_back(m);
      }

    /*! After this operation, the rangeset contains the union of itself
        and \a other. */
    void unite (const rangeset &other)
      {
      rtype tmp;
      generalUnion (r,other.r,false,false,tmp);
      std::swap(r,tmp);
      }
    /*! After this operation, the rangeset contains the intersection of itself
        and \a other. */
    void intersect (const rangeset &other)
      {
      rtype tmp;
      generalUnion (r,other.r,true,true,tmp);
      std::swap(r,tmp);
      }
    /*! After this operation, the rangeset contains the union of itself
        with the inverse of \a other. */
    void subtract (const rangeset &other)
      {
      rtype tmp;
      generalUnion (r,other.r,true,false,tmp);
      std::swap(r,tmp);
      }
    /*! After this operation, the rangeset contains the union of \a a
        and \a b. */
    void setToUnion (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,false,false,r); }
    /*! After this operation, the rangeset contains the intersection of \a a
        and \a b. */
    void setToIntersection (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,true,true,r); }
    /*! After this operation, the rangeset contains the union of \a a
        with the inverse of \a b. */
    void setToDifference (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,true,false,r); }

    /*! Returns the index of the interval containing \a v; if no such interval
        exists, -1 is returned. */
    tdiff findInterval (const T &v) const
      {
      tdiff res = iiv(v);
      return (res&1) ? -1 : res>>1;
      }

    /*! Returns \a true if the rangeset is identical to \a other, else \a false.
        */
    bool equals (const rangeset &other) const
      { return r==other.data(); }

    /*! Returns \a true if the rangeset contains all values in the range
        \a [a;b[, else \a false. */
    bool containsAll (T a,T b) const
      {
      tdiff res=iiv(a);
      if (res&1) return false;
      return (b<=r[res+1]);
      }
    /*! Returns \a true if the rangeset contains the value \a v,
        else \a false. */
    bool contains (T v) const
      { return !(iiv(v)&1); }
    /*! Returns \a true if the rangeset contains all values stored in \a other,
        else \a false. */
    bool contains (const rangeset &other) const
      {
      tsize im=0, em=r.size();
      for (tsize i=0; i<other.r.size(); i+=2)
        {
        T a=other.r[i], b=other.r[i+1];
        while ((im!=em) && (r[im+1] < a)) im+=2;
        if (im==em) return false;
        if ((r[im]>a) || (r[im+1]<b)) return false;
        }
      return true;
      }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const rangeset<T> &rs)
  {
  os << "{ ";
  for (tsize i=0; i<rs.size(); ++i)
    os << "["<<rs.ivbegin(i)<<";"<<rs.ivend(i)<<"[ ";
  return os << "}";
  }

#endif

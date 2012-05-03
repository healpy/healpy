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

/*! Class for holding a single interval of integer numbers. */
template<typename T> class interval
  {
  private:
    T A, B;
    void check() const
      { planck_assert(A<=B, "inconsistent interval"); }

  public:
    /*! Constructs the interval \a [v;v+1[. */
    explicit interval (const T &v) : A(v), B(v) { ++B; }
    /*! Constructs the interval \a [a_;b_[. */
    interval (const T &a_, const T &b_) : A(a_), B(b_) { check(); }
    /*! Copy constructor. */
    interval (const interval &other) : A(other.A), B(other.B) {}
    /*! Assignment operator. */
    const interval &operator= (const interval &other)
      { A=other.A; B=other.B; return *this; }

    /*! Returns the number of elements in the interval. */
    tsize size() const { return B-A; }
    /*! Returns \a true if the interval is empty, else \a false. */
    bool empty() const { return A==B; }

    /*! Returns \a true if \a v lies inside the interval, else \a false. */
    bool contains (const T &v) const { return (A<=v) && (v<B); }
    /*! Returns \a true if \a i lies completely inside the interval,
        else \a false. */
    bool contains (const interval &i) const { return (A<=i.A) && (i.B<=B); }
    /*! Returns \a true if \a i lies completely inside the interval and does
        not contain the limits of the interval, else \a false. */
    bool strictContains (const interval &i) const
      { return (A<i.A) && (i.B<B); }

    /*! Equality check. */
    bool operator== (const interval &i) const { return (A==i.A) && (B==i.B); }
//    bool operator< (const T &v) const { return B<v; }
    bool operator< (const interval &i) const { return B<i.A; }

    bool touching (const interval &i) const
      { return std::max(A,i.A) <= std::min(B, i.B); }

    void merge (const interval &i)
      { A=std::min(A,i.A); B=std::max(B,i.B); }

    void setA (const T &v)
      { A=v; check(); }
    void setB (const T &v)
      { B=v; check(); }

    operator bool() const { return A!=B; }

    T a() const { return A; }
    T b() const { return B; }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const interval<T> &iv)
  { return os << "[" << iv.a() << "," << iv.b() << ")"; }

/*! Class for storing sets of ranges of integer numbers */
template<typename T> class rangeset
  {
  private:
    typedef std::vector<interval<T> > rtype;
    typedef typename rtype::iterator iterator;
    typedef typename rtype::const_iterator c_iterator;
    rtype r;

    static void generalUnion (const rtype &a, const rtype &b,
      bool flip_a, bool flip_b, rtype &c)
      {
      planck_assert((&c!=&a)&&(&c!=&b), "cannot overwrite the rangeset");
      c.clear();
      tsize out=0;
      T vsave;
      bool state_a=flip_a, state_b=flip_b,
           state_res=state_a||state_b;
      tsize ia=0, ea=2*a.size(), ib=0, eb=2*b.size();
      bool runa = ia!=ea, runb = ib!=eb;
      while(runa||runb)
        {
        bool adv_a=false, adv_b=false;
        T val,va,vb;
        if (runa) va = (ia&1) ? a[ia>>1].b() : a[ia>>1].a();
        if (runb) vb = (ib&1) ? b[ib>>1].b() : b[ib>>1].a();
        if (runa && (!runb || (va<=vb))) { adv_a=true; val=va; }
        if (runb && (!runa || (vb<=va))) { adv_b=true; val=vb; }
        if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
        if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
        bool tmp=state_a||state_b;
        if (tmp!=state_res)
          {
          if ((++out)&1)
            vsave=val;
          else
            c.push_back(interval<T>(vsave,val));
          state_res = tmp;
          }
        }
      }

  public:
    /*! Removes all rangeset entries. */
    void clear() { r.clear(); }
    /*! Reserves space for \a n ranges. */
    void reserve(tsize n) { r.reserve(n); }
    /*! Returns the current number of ranges. */
    tsize size() const { return r.size(); }
    /*! Returns the current vector of ranges. */
    const rtype &data() const { return r; }

    /*! Returns the \a ith range. */
    const interval<T> &operator[] (tsize i) const { return r[i]; }

    /*! Appends \a iv to the rangeset. All values in \a iv must be larger
        than all values in the rangeset. */
    void append(const interval<T> &iv)
      {
      if (iv.empty()) return;
      if ((!r.empty()) && (iv.a()<=r.back().b()))
        {
        planck_assert (iv.a()==r.back().b(),"bad append operation");
        r.back().setB(iv.b());
        }
      else
        r.push_back(iv);
      }
    /*! Appends \a [v1;v2[ to the rangeset. \a v1 must be larger
        than all values in the rangeset. */
    void append(const T &v1, const T &v2) { append(interval<T>(v1,v2)); }
    /*! Appends \a [v;v+1[ to the rangeset. \a v must be larger
        than all values in the rangeset. */
    void append(const T &v) { append(interval<T>(v)); }

    /*! Appends \a iv to the rangeset. All values in \a iv must be larger
        than the minimum of the last range in the rangeset. */
    void appendRelaxed(const interval<T> &iv)
      {
      if (iv.empty()) return;
      if ((!r.empty()) && (iv.a()<=r.back().b()))
        {
        planck_assert (iv.a()>=r.back().a(),"bad append operation");
        if (iv.b()>r.back().b()) r.back().setB(iv.b());
        }
      else
        r.push_back(iv);
      }
    /*! Appends \a [v1;v2[ to the rangeset. \a v1 must be larger
        than the minimum of the last range in the rangeset. */
    void appendRelaxed(const T &v1, const T &v2)
      { appendRelaxed(interval<T>(v1,v2)); }
    /*! Appends \a [v;v+1[ to the rangeset. \a v must be larger
        than the minimum of the last range in the rangeset. */
    void appendRelaxed(const T &v)
      { appendRelaxed(interval<T>(v)); }

    /*! After this operation, the rangeset contains the union of itself
        with \a iv. */
    void add(const interval<T> &iv)
      {
      if (iv.empty()) return;
      iterator i=std::lower_bound (r.begin(),r.end(),iv);
      if (i==r.end() || !i->touching(iv))
        r.insert(i,iv);
      else
        {
        i->merge(iv);
        iterator j=i;
        while (++j!=r.end() && i->touching(*j))
          {}
        i->setB(std::max(i->b(),(--j)->b()));
        r.erase(++i,++j);
        }
      }
    /*! After this operation, the rangeset contains the union of itself
        with \a [v1;v2[. */
    void add(const T &v1, const T &v2) { add(interval<T>(v1,v2)); }
    /*! After this operation, the rangeset contains the union of itself
        with \a [v;v+1[. */
    void add(const T &v) { add(interval<T>(v)); }

    /*! Removes all values within \a iv from the rangeset. */
    void remove(const interval<T> &iv)
      {
      if (iv.empty()) return;
      iterator i=std::lower_bound (r.begin(),r.end(),iv);
      if (i==r.end() || i->a() >= iv.b())
        return;
      if (*i==iv)
        r.erase(i);
      else if (i->strictContains(iv))
        {
        interval<T> iv1(i->a(), iv.a()), iv2(iv.b(), i->b());
        *i=iv2;
        r.insert(i,iv1);
        }
      else
        {
        if (i->a() < iv.a())
          {
          i->setB(iv.a());
          ++i;
          }
        iterator j=i;
        while (j!=r.end() && iv.contains(*j))
          ++j;
        if (j!=r.end() && iv.b() > j->a())
          j->setA(iv.b());
        r.erase(i,j);
        }
      }
    /*! Removes all values within \a [v1;v2] from the rangeset. */
    void remove(const T &v1, const T &v2) { remove(interval<T>(v1,v2)); }
    /*! Removes the value \a v from the rangeset. */
    void remove(const T &v) { remove(interval<T>(v)); }

    /*! Returns the total number of elements in the rangeset. */
    T nval() const
      {
      T result=T(0);
      for (tsize i=0; i<r.size(); ++i)
        result+=r[i].size();
      return result;
      }

    /*! After this opration, \a res contains all elements of the rangeset
        in ascending order. */
    void toVector (std::vector<T> &res) const
      {
      res.clear();
      res.reserve(nval());
      for (tsize i=0; i<r.size(); ++i)
        for (T m(r[i].a()); m<r[i].b(); ++m)
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
    /*! After this operation, the rangeset contains the union of \a r1
        and \a r2. */
    void setToUnion (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,false,false,r); }
    /*! After this operation, the rangeset contains the intersection of \a r1
        and \a r2. */
    void setToIntersection (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,true,true,r); }
    /*! After this operation, the rangeset contains the union of \a r1
        with the inverse of \a r2. */
    void setToDifference (const rangeset &a, const rangeset &b)
      { generalUnion (a.r,b.r,true,false,r); }

    /*! Returns the index of the interval containing \a v; if no such interval
        exists, -1 is returned. */
    tdiff findInterval (const T &v) const
      {
      interval<T> iv(v);
      c_iterator i=std::lower_bound (r.begin(),r.end(),iv);
      if (i==r.end() || !i->contains(v))
        return -1;
      return i-r.begin();
      }

    /*! Returns \a true if the rangeset is identical to \a other, else \a false.
        */
    bool equals (const rangeset &other) const
      { return r==other.data(); }
    /*! Returns \a true if the rangeset contains all values stored in \a other,
        else \a false. */
    bool contains (const rangeset &other) const
      {
      c_iterator im=r.begin(), em=r.end();
      for (tsize i=0; i<other.size(); ++i)
        {
        T a=other[i].a(), b=other[i].b();
        while ((im!=em) && (im->b() < a)) ++im;
        if (im==em) return false;
        if ((im->a()>a) || (im->b()<b)) return false;
        }
      return true;
      }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const rangeset<T> &rs)
  {
  os << "{ ";
  for (tsize i=0; i<rs.size(); ++i)
    os << rs[i] << " ";
  return os << "}";
  }

#endif

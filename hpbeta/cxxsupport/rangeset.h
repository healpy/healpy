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

#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include "datatypes.h"
#include "error_handling.h"

template<typename T> class interval
  {
  private:
    T A, B;
    void check() const
      { planck_assert(A<=B, "inconsistent interval"); }

  public:
    explicit interval (const T &v) : A(v), B(v) { ++B; }
    interval (const T &a_, const T &b_) : A(a_), B(b_) { check(); }
    interval (const interval &other) : A(other.A), B(other.B) {}
    const interval &operator= (const interval &other)
      { A=other.A; B=other.B; return *this; }

    tsize size() const { return B-A; }
    bool empty() const { return A==B; }

    bool contains (const T &v) const { return (A<=v) && (v<B); }
    bool contains (const interval &i) const { return (A<=i.A) && (i.B<=B); }
    bool strictContains (const interval &i) const
      { return (A<i.A) && (i.B<B); }

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

  public:
    void clear() { r.clear(); }
    void reserve(tsize n) { r.reserve(n); }
    tsize size() const { return r.size(); }

    const interval<T> &operator[] (tsize i) const { return r[i]; }

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
    void append(const T &v1, const T &v2) { append(interval<T>(v1,v2)); }
    void append(const T &v) { append(interval<T>(v)); }

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
    void appendRelaxed(const T &v1, const T &v2)
      { appendRelaxed(interval<T>(v1,v2)); }
    void appendRelaxed(const T &v)
      { appendRelaxed(interval<T>(v)); }

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
    void add(const T &v1, const T &v2) { add(interval<T>(v1,v2)); }
    void add(const T &v) { add(interval<T>(v)); }

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
    void remove(const T &v1, const T &v2) { remove(interval<T>(v1,v2)); }
    void remove(const T &v) { remove(interval<T>(v)); }

    T nval() const
      {
      T result=T(0);
      for (tsize i=0; i<r.size(); ++i)
        result+=r[i].size();
      return result;
      }

    void toVector (std::vector<T> &res) const
      {
      res.clear();
      res.reserve(nval());
      for (tsize i=0; i<r.size(); ++i)
        for (T m(r[i].a()); m<r[i].b(); ++m)
          res.push_back(m);
      }

    void setToUnion (const rangeset &r1, const rangeset &r2)
      {
      c_iterator i1 = r1.r.begin(), e1 = r1.r.end(),
                 i2 = r2.r.begin(), e2 = r2.r.end();

      clear();

      bool run1 = i1!=e1, run2 = i2!=e2;

      while(run1 || run2)
        {
        if (run1 && (!run2 || i1->b()<i2->a()))
          {
          appendRelaxed(*i1);
          run1 = (++i1)!=e1;
          }
        else if (run2 && (!run1 || i2->b()<i1->a()))
          {
          appendRelaxed(*i2);
          run2 = (++i2)!=e2;
          }
        else if(run1 && run2)
          {
          appendRelaxed(min(i1->a(),i2->a()),max(i1->b(),i2->b()));
          run1 = (++i1)!=e1;
          run2 = (++i2)!=e2;
          }
        else
          {
          planck_fail("internal error");
          }
        }
      }

    tdiff findInterval (const T &v) const
      {
      interval<T> iv(v);
      c_iterator i=std::lower_bound (r.begin(),r.end(),iv);
      if (i==r.end() || !i->contains(v))
        return -1;
      return i-r.begin();
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

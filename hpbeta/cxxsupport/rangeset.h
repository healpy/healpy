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
      { planck_assert(a<=b, "inconsistent interval"); }

  public:
    const T &a, &b;

    explicit interval (const T &v) : A(v), B(v), a(A), b(B) { ++B; }
    interval (const T &a_, const T &b_) : A(a_), B(b_), a(A), b(B) { check(); }
    interval (const interval &other) : A(other.a), B(other.b), a(A), b(B) {}
    const interval &operator= (const interval &other)
      { A=other.a; B=other.b; return *this; }

    tsize size() const { return b-a; }
    bool empty() const { return a==b; }

    bool contains (const T &v) const { return (a<=v) && (v<b); }
    bool contains (const interval &i) const { return (a<=i.a) && (i.b<=b); }
    bool strictContains (const interval &i) const
      { return (a<i.a) && (i.b<b); }

    bool operator== (const interval &i) const { return (a==i.a) && (b==i.b); }
//    bool operator< (const T &v) const { return b<v; }
    bool operator< (const interval &i) const { return b<i.a; }

    bool touching (const interval &i) const
      { return std::max(a,i.a) <= std::min(b, i.b); }

    void merge (const interval &i)
      { A=std::min(a,i.a); B=std::max(b,i.b); }

    void setA (const T &v)
      { A=v; check(); }
    void setB (const T &v)
      { B=v; check(); }

    operator bool() const { return a!=b; }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const interval<T> &iv)
  { return os << "[" << iv.a << "," << iv.b << ") "; }

/*! Class for storing sets of ranges of integer numbers */
template<typename T> class rangeset
  {
  private:
    typedef std::vector<interval<T> > rtype;
    typedef typename rtype::iterator iterator;
    rtype r;

  public:
    void clear() { r.clear(); }
    void reserve(tsize n) { r.reserve(n); }
    tsize size() const { return r.size(); }

    const interval<T> &operator[] (tsize i) const { return r[i]; }

    void append(const interval<T> &iv)
      {
      if (iv.empty()) return;
      if ((!r.empty()) && (iv.a<=r.back().b))
        {
        planck_assert (iv.a==r.back().b,"bad append operation");
        r.back().setB(iv.b);
        }
      else
        r.push_back(iv);
      }
    void append(const T &v1, const T &v2) { append(interval<T>(v1,v2)); }
    void append(const T &v) { append (interval<T>(v)); }

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
        i->b=std::max(i->b,(--j)->b);
        r.erase(++i,++j);
        }
      }
    void add(const T &v1, const T &v2) { add(interval<T>(v1,v2)); }
    void add(const T &v) { add(interval<T>(v)); }

    void remove(const interval<T> &iv)
      {
      if (iv.empty()) return;
      iterator i=std::lower_bound (r.begin(),r.end(),iv);
      if (i==r.end() || i->a >= iv.b)
        return;
      if (*i==iv)
        r.erase(i);
      else if (i->strictContains(iv))
        {
        interval<T> iv1(i->a, iv.a), iv2(iv.b, i->b);
        *i=iv2;
        r.insert(i,iv1);
        }
      else
        {
        if (i->a < iv.a)
          {
          i->setB(iv.a);
          ++i;
          }
        iterator j=i;
        while (j!=r.end() && iv.contains(*j))
          ++j;
        if (j!=r.end() && iv.b > j->a)
          j->setA(iv.b);
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
        for (T m(r[i].a); m<r[i].b; ++m)
          res.push_back(m);
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

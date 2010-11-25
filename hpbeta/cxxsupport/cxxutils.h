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

/*! \file cxxutils.h
 *  Various convenience functions used by the Planck LevelS package.
 *
 *  Copyright (C) 2002 - 2010 Max-Planck-Society
 *  \author Martin Reinecke \author Reinhard Hell
 */

#ifndef PLANCK_CXXUTILS_H
#define PLANCK_CXXUTILS_H

#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include "error_handling.h"
#include "datatypes.h"

/*! \defgroup mathutilsgroup Mathematical helper functions */
/*! \{ */

/*! Returns \e true if | \a a-b | < \a epsilon * | \a b |, else \e false. */
template<typename F> inline bool approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < (epsilon*abs(b));
  }

/*! Returns \e true if | \a a-b | < \a epsilon, else \e false. */
template<typename F> inline bool abs_approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < epsilon;
  }

/*! Returns the largest integer which is smaller than (or equal to) \a arg. */
template<typename I, typename F> inline I ifloor (F arg)
  {
  using namespace std;
  return I(floor(arg));
  }

/*! Returns the integer which is nearest to \a arg. */
template<typename I, typename F> inline I nearest (F arg)
  { return ifloor<I>(arg+0.5); }

/*! Returns \a v1+v2 if \a v1<0, \a v1-v2 if \a v1>=v2, else \a v1.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename T> inline T weak_modulo (T v1, T v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1-v2)) : (v1+v2); }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline double fmodulo (double v1, double v2)
  {
  using namespace std;
  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
  }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename I> inline I imodulo (I v1, I v2)
  { I v=v1%v2; return (v>=0) ? v : v+v2; }

/*! Returns -1 if \a signvalue is negative, else +1. */
template<typename T> inline T sign (const T& signvalue)
  { return (signvalue>=0) ? 1 : -1; }

/*! Returns \a val*pow(-1,m) */
template<typename T, typename I> inline T xpow (I m, T val)
  { return (m&1) ? -val : val; }

template <typename I, bool g4> struct isqrt_helper__
  {};
template <typename I> struct isqrt_helper__ <I, false>
  {
  static uint32 isqrt (I arg)
    {
    using namespace std;
    return uint32 (sqrt(arg+0.5));
    }
  };
template <typename I> struct isqrt_helper__ <I, true>
  {
  static uint32 isqrt (I arg)
    {
    using namespace std;
    long double arg2 = static_cast<long double>(arg)+0.5;
    return uint32 (sqrt(arg2));
    }
  };

/*! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1). */
template<typename I> inline uint32 isqrt (I arg)
  { return isqrt_helper__<I,(sizeof(I)>4)>::isqrt(arg); }

/*! Returns the largest integer \a n that fulfills \a 2^n<=arg. */
template<typename I> inline unsigned int ilog2 (I arg)
  {
  unsigned int res=0;
  while (arg > 0x0000FFFF) { res+=16; arg>>=16; }
  if (arg > 0x000000FF) { res|=8; arg>>=8; }
  if (arg > 0x0000000F) { res|=4; arg>>=4; }
  if (arg > 0x00000003) { res|=2; arg>>=2; }
  if (arg > 0x00000001) { res|=1; }
  return res;
  }

/*! Returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0. */
inline double safe_atan2 (double y, double x)
  {
  using namespace std;
  return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
  }

/*! Helper function for linear interpolation (or extrapolation).
    The array must be ordered in ascending order; no two values may be equal. */
template<typename T, typename Iter, typename Comp> inline void interpol_helper
  (const Iter &begin, const Iter &end, const T &val, Comp comp, tsize &idx,
  T &frac)
  {
  using namespace std;
  planck_assert((end-begin)>1,"sequence too small for interpolation");
  idx = lower_bound(begin,end,val,comp)-begin;
  if (idx>0) --idx;
  idx = min(tsize(end-begin-2),idx);
  frac = (val-begin[idx])/(begin[idx+1]-begin[idx]);
  }

/*! Helper function for linear interpolation (or extrapolation).
    The array must be ordered in ascending order; no two values may be equal. */
template<typename T, typename Iter> inline void interpol_helper
  (const Iter &begin, const Iter &end, const T &val, tsize &idx, T &frac)
  { interpol_helper (begin,end,val,std::less<T>(),idx,frac); }

/*! \} */

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c)
  { return (a==b) && (a==c); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d)
  { return (a==b) && (a==c) && (a==d); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d, const T &e)
  { return (a==b) && (a==c) && (a==d) && (a==e); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d, const T &e, const T &f)
  { return (a==b) && (a==c) && (a==d) && (a==e) && (a==f); }

template<typename It, typename Comp> class IdxComp__
  {
  private:
    It begin;
    Comp comp;
  public:
    IdxComp__ (It begin_, Comp comp_): begin(begin_), comp(comp_) {}
    bool operator() (std::size_t a, std::size_t b) const
      { return comp(*(begin+a),*(begin+b)); }
  };

/*! Performs an indirect sort on the supplied iterator range and returns in
    \a idx a \a vector containing the indices of the smallest, second smallest,
    third smallest, etc. element, according to \a comp. */
template<typename It, typename T2, typename Comp>
  inline void buildIndex (It begin, It end, std::vector<T2> &idx, Comp comp)
  {
  using namespace std;
  T2 num=end-begin;
  idx.resize(num);
  for (T2 i=0; i<num; ++i) idx[i] = i;
  sort (idx.begin(),idx.end(),IdxComp__<It,Comp>(begin,comp));
  }

/*! Performs an indirect sort on the supplied iterator range and returns in
    \a idx a \a vector containing the indices of the smallest, second smallest,
    third smallest, etc. element. */
template<typename It, typename T2> inline void buildIndex (It begin, It end,
  std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  buildIndex(begin,end,idx,less<T>());
  }

/*! Sorts the supplied iterator range according to the order given by \a idx.
    The operation is done out of place and requires temporary extra storage. */
template<typename It, typename T2> inline void sortByIndex (It begin, It end,
  const std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  T2 num=end-begin;
  T *tmp= new T[num];
  for (T2 i=0; i<num; ++i) tmp[i]=*(begin+i);
  for (T2 i=0; i<num; ++i) *(begin+i) = tmp[idx[i]];
  delete[] tmp;
  }

/*! Sorts the supplied iterator range according to the order given by \a idx.
    The operation is done in place. */
template<typename It, typename T2> inline void sortByIndex_inplace
  (It begin, It end, const std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  T2 num=end-begin;
  vector<bool> done(num,false);
  T2 cnt=0;
  while (cnt<num)
    {
    if (!done[cnt]) // new cycle
      {
      T tmp(*(begin+cnt));
      T2 cnt2 = cnt;
      T2 cnt3 = idx[cnt];
      while (cnt3!=cnt)
        {
        done[cnt2]=true;
        *(begin+cnt2)=*(begin+cnt3);
        cnt2=cnt3;
        cnt3=idx[cnt3];
        }
      *(begin+cnt2) = tmp;
      }
    ++cnt;
    }
  }

template<typename It, typename Comp> inline void indirectSort (It begin, It end,
  Comp comp)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  vector<std::size_t> idx;
  buildIndex (begin,end,idx,comp);
  sortByIndex (begin,end,idx);
  }

template<typename It> inline void indirectSort (It begin, It end)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  indirectSort(begin,end,less<T>());
  }

/*! \defgroup stringutilsgroup String handling helper functions */
/*! \{ */

/*! Returns the string \a orig without leading and trailing whitespace. */
std::string trim (const std::string &orig);

/*! Returns a string containing the text representation of \a x.
    Care is taken that no information is lost in the conversion. */
template<typename T> std::string dataToString(const T &x);
template<> std::string dataToString (const bool &x);
template<> std::string dataToString (const std::string &x);
template<> std::string dataToString (const float &x);
template<> std::string dataToString (const double &x);

/*! Returns a string containing the text representation of \a x, padded
    with leading zeroes to \a width characters. */
std::string intToString(int64 x, tsize width);

/*! Reads a value of a given datatype from a string */
template<typename T> void stringToData (const std::string &x, T &value);
template<> void stringToData (const std::string &x, std::string &value);
template<> void stringToData (const std::string &x, bool &value);

/*! Reads a value of a given datatype from a string */
template<typename T> inline T stringToData (const std::string &x)
  { T result; stringToData(x,result); return result; }

/*! Parses the file \a filename and returns the key/value pairs in \a dict. */
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

/*! Case-insensitive string comparison
    Returns \a true, if \a a and \a b differ only in capitalisation,
    else \a false. */
bool equal_nocase (const std::string &a, const std::string &b);

/*! Returns lowercase version of \a input. */
std::string tolower(const std::string &input);

/*! \} */

/*! Prints a banner containing \a name, as well as some information about the
    source code and the parallelisation techniques enabled. */
void announce (const std::string &name);

/*! Prints a banner containing \a name and checks if \a argc==argc_expected.
    If not, a usage description is given and the program is terminated. */
void module_startup (const std::string &name, int argc, const char **argv,
  int argc_expected, const std::string &argv_expected, bool verbose=true);

/*! Divides the index range [\a glo; \a ghi) into \a nshares approximately
    equal parts, and returns the sub-range [\a lo; \a hi) of the
    part with the number \a myshare (first part has index 0). */
void calcShareGeneral (int64 glo, int64 ghi, int64 nshares, int64 myshare,
  int64 &lo, int64 &hi);

/*! Helper class for dividing a range of work items into chunks of specified
    size. */
class chunkMaker
  {
  private:
    uint64 s_full, s_chunk, offset;

  public:
    /*! Creates an object that produces chunk information for \a s_full_
        work items and a desired chunk size of \a s_chunk_. */
    chunkMaker (uint64 s_full_, uint64 s_chunk_)
      : s_full(s_full_), s_chunk(s_chunk_), offset(0) {}

    /*! Returns the total number of chunks. */
    uint64 nchunks() const
      { return (s_full+s_chunk-1)/s_chunk; }

    /*! Returns the start index of the next chunk in \a start, and its size
        in \a size. If all chunks have been processed already, the return
        value is \a false, else \a true. */
    bool getNext (uint64 &start, uint64 &size)
      {
      using namespace std;
      if (offset>=s_full) return false;
      start=offset;
      size=min(s_chunk,s_full-offset);
      offset+=s_chunk;
      return true;
      }
  };

/*! Tries to split \a inp into a white-space separated list of values of
    type \a T, and appends them to \a list. */
template<typename T> void split (const std::string &inp, std::vector<T> &list);

/*! Resizes \a container to zero and releases its memory. Typically used for
    std::vector.
    Taken from http://www.gotw.ca/gotw/054.htm */
template<typename T> inline void releaseMemory (T &container)
  { T().swap(container); }

/*! Releases all unused memory that \a container might have. Typically used for
    std::vector.
    Taken from http://www.gotw.ca/gotw/054.htm */
template<typename T> inline void shrinkToFit (T &container)
  { T(container).swap(container); }

/*! Breaks the string \a inp into tokens separated by \a delim, and returns them
    in \a list. */
void tokenize (const std::string &inp, char delim,
  std::vector<std::string> &list);

#endif

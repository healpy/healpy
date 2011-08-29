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

/*! \file arr.h
 *  Various high-performance array classes used by the Planck LevelS package.
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ARR_H
#define PLANCK_ARR_H

#include <algorithm>
#include <cstdlib>
#include "alloc_utils.h"
#include "datatypes.h"
#include "math_utils.h"

/*! \defgroup arraygroup Array classes */
/*! \{ */

/*! View of a 1D array */
template <typename T> class arr_ref
  {
  protected:
    tsize s;
    T *d;

  public:
    /*! Constructs an \a arr_ref of size \a s_, starting at \a d_. */
    arr_ref(T *d_, tsize s_) : s(s_),d(d_) {}

    /*! Returns the current array size. */
    tsize size() const { return s; }

    /*! Writes \a val into every element of the array. */
    void fill (const T &val)
      { for (tsize m=0; m<s; ++m) d[m]=val; }

    /*! Returns a reference to element \a n */
    template<typename T2> T &operator[] (T2 n) {return d[n];}
    /*! Returns a constant reference to element \a n */
    template<typename T2> const T &operator[] (T2 n) const {return d[n];}

    /*! Returns a pointer to the first array element, or NULL if the array
        is zero-sized. */
    T *begin() { return d; }
    /*! Returns a pointer to the one-past-last array element, or NULL if the
        array is zero-sized. */
    T *end() { return d+s; }
    /*! Returns a constant pointer to the first array element, or NULL if the
        array is zero-sized. */
    const T *begin() const { return d; }
    /*! Returns a constant pointer to the one-past-last array element, or NULL
        if the array is zero-sized. */
    const T *end() const { return d+s; }

    /*! Copies all array elements to \a ptr. */
    template<typename T2> void copyToPtr (T *ptr) const
      { for (tsize m=0; m<s; ++m) ptr[m]=d[m]; }

    /*! Sorts the elements in the array, in ascending order. */
    void sort()
      { std::sort (d,d+s); }

    /*! Sorts the elements in the array, such that \a comp(d[i],d[j])==true
        for \a i<j. */
    template<typename Comp> void sort(Comp comp)
      { std::sort (d,d+s,comp); }

    /*! Helper function for linear interpolation (or extrapolation).
        \a idx and \a val are computed such that
        \a val=d[idx]+frac*(d[idx+1]-d[idx]). If \a val<d[0], \a frac will be
        negative, if \a val>d[s-1], frac will be larger than 1. In all other
        cases \a 0<=frac<=1.

        The array must be ordered in ascending order; no two values may be
        equal. */
    void interpol_helper (const T &val, tsize &idx, double &frac) const
      { ::interpol_helper (d, d+s, val, idx, frac); }

    /*! Helper function for linear interpolation (or extrapolation).
        \a idx and \a val are computed such that
        \a val=d[idx]+frac*(d[idx+1]-d[idx]). If \a comp(val,d[0])==true,
        \a frac will be negative, if \a comp(val,d[s-1])==false, frac will be
        larger than 1. In all other cases \a 0<=frac<=1.

        The array must be ordered such that \a comp(d[i],d[j])==true
        for \a i<j; no two values may be equal. */
    template<typename Comp> void interpol_helper (const T &val, Comp comp,
      tsize &idx, double &frac) const
      { ::interpol_helper (d, d+s, val, comp, idx, frac); }

    /*! Returns the minimum and maximum entry in \a minv and \a maxv,
        respectively. Throws an exception if the array is zero-sized. */
    void minmax (T &minv, T &maxv) const
      {
      planck_assert(s>0,"trying to find min and max of a zero-sized array");
      minv=maxv=d[0];
      for (tsize m=1; m<s; ++m)
        {
        if (d[m]<minv) minv=d[m];
        else if (d[m]>maxv) maxv=d[m];
        }
      }

    /*! Returns \a true, if \a val is found in the array, else \a false. */
    bool contains (const T &val) const
      {
      for (tsize m=0; m<s; ++m)
        if (d[m]==val) return true;
      return false;
      }

    /*! Returns the index of the first occurrence of \a val in the array.
        If it is not found, an exception is thrown. */
    tsize find (const T &val) const
      {
      for (tsize m=0; m<s; ++m)
        if (d[m]==val) return m;
      planck_fail ("entry not found in array");
      }

    /*! Returns \a true if the array has the same size as \a other and all
        elements of both arrays are equal, else \a false. */
    bool contentsEqual(const arr_ref &other) const
      {
      if (s!=other.s) return false;
      for (tsize i=0; i<s; ++i)
        if (d[i]!=other.d[i]) return false;
      return true;
      }
  };

/*! An array whose size is known at compile time. Very useful for storing
    small arrays on the stack, without need for \a new and \a delete(). */
template <typename T, tsize sz> class fix_arr
  {
  private:
    T d[sz];

  public:
    /*! Returns the size of the array. */
    tsize size() const { return sz; }

    /*! Returns a reference to element \a n */
    template<typename T2> T &operator[] (T2 n) {return d[n];}
    /*! Returns a constant reference to element \a n */
    template<typename T2> const T &operator[] (T2 n) const {return d[n];}
  };


/*! One-dimensional array type, with selectable storage management. */
template <typename T, typename storageManager> class arrT: public arr_ref<T>
  {
  private:
    storageManager stm;
    bool own;

    void reset()
      { this->d=0; this->s=0; own=true; }

  public:
    /*! Creates a zero-sized array. */
    arrT() : arr_ref<T>(0,0), own(true) {}
    /*! Creates an array with \a sz entries. */
    arrT(tsize sz) : arr_ref<T>(stm.alloc(sz),sz), own(true) {}
    /*! Creates an array with \a sz entries, and initializes them with
        \a inival. */
    arrT(tsize sz, const T &inival) : arr_ref<T>(stm.alloc(sz),sz), own(true)
      { this->fill(inival); }
    /*! Creates an array with \a sz entries, which uses the memory pointed
        to by \a ptr.
        \note \a ptr will <i>not</i> be deallocated by the destructor.
        \warning Only use this if you REALLY know what you are doing.
        In particular, this is only safely usable if
          <ul>
          <li>\a T is a POD type</li>
          <li>\a ptr survives during the lifetime of the array object</li>
          <li>\a ptr is not subject to garbage collection</li>
          </ul>
        Other restrictions may apply. You have been warned. */
    arrT (T *ptr, tsize sz): arr_ref<T>(ptr,sz), own(false) {}
    /*! Creates an array which is a copy of \a orig. The data in \a orig
        is duplicated. */
    arrT (const arrT &orig): arr_ref<T>(stm.alloc(orig.s),orig.s), own(true)
      { for (tsize m=0; m<this->s; ++m) this->d[m] = orig.d[m]; }
    /*! Frees the memory allocated by the object. */
    ~arrT() { if (own) stm.dealloc(this->d); }

    /*! Allocates space for \a sz elements. The content of the array is
        undefined on exit. \a sz can be 0. If \a sz is the
        same as the current size, no reallocation is performed. */
    void alloc (tsize sz)
      {
      if (sz==this->s) return;
      if (own) stm.dealloc(this->d);
      this->s = sz;
      this->d = stm.alloc(sz);
      own = true;
      }
    /*! Allocates space for \a sz elements. If \a sz is the
        same as the current size, no reallocation is performed.
        All elements are set to \a inival. */
    void allocAndFill (tsize sz, const T &inival)
      { alloc(sz); this->fill(inival); }
    /*! Deallocates the memory held by the array, and sets the array size
        to 0. */
    void dealloc() {if (own) stm.dealloc(this->d); reset();}
    /*! Resizes the array to hold \a sz elements. The existing content of the
        array is copied over to the new array to the extent possible.
        \a sz can be 0. If \a sz is the same as the current size, no
        reallocation is performed. */
    void resize (tsize sz)
      {
      using namespace std;
      if (sz==this->s) return;
      T *tmp = stm.alloc(sz);
      for (tsize m=0; m<min(sz,this->s); ++m)
        tmp[m]=this->d[m];
      if (own) stm.dealloc(this->d);
      this->s = sz;
      this->d = tmp;
      own = true;
      }

    /*! Changes the array to be a copy of \a orig. */
    arrT &operator= (const arrT &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s);
      for (tsize m=0; m<this->s; ++m) this->d[m] = orig.d[m];
      return *this;
      }

    /*! Reserves space for \a sz elements, then copies \a sz elements
        from \a ptr into the array. */
    template<typename T2> void copyFromPtr (const T2 *ptr, tsize sz)
      {
      alloc(sz);
      for (tsize m=0; m<this->s; ++m) this->d[m]=ptr[m];
      }

    /*! Assigns the contents and size of \a other to the array.
        \note On exit, \a other is zero-sized! */
    void transfer (arrT &other)
      {
      if (own) stm.dealloc(this->d);
      this->d=other.d;
      this->s=other.s;
      own=other.own;
      other.reset();
      }
    /*! Swaps contents and size with \a other. */
    void swap (arrT &other)
      {
      std::swap(this->d,other.d);
      std::swap(this->s,other.s);
      std::swap(own,other.own);
      }
  };

/*! One-dimensional array type. */
template <typename T>
  class arr: public arrT<T,normalAlloc__<T> >
  {
  public:
    /*! Creates a zero-sized array. */
    arr() : arrT<T,normalAlloc__<T> >() {}
    /*! Creates an array with \a sz entries. */
    arr(tsize sz) : arrT<T,normalAlloc__<T> >(sz) {}
    /*! Creates an array with \a sz entries, and initializes them with
        \a inival. */
    arr(tsize sz, const T &inival) : arrT<T,normalAlloc__<T> >(sz,inival) {}
    /*! Creates an array with \a sz entries, which uses the memory pointed
        to by \a ptr.
        \note \a ptr will <i>not</i> be deallocated by the destructor.
        \warning Only use this if you REALLY know what you are doing.
        In particular, this is only safely usable if
          <ul>
          <li>\a T is a POD type</li>
          <li>\a ptr survives during the lifetime of the array object</li>
          <li>\a ptr is not subject to garbage collection</li>
          </ul>
        Other restrictions may apply. You have been warned. */
    arr (T *ptr, tsize sz): arrT<T,normalAlloc__<T> >(ptr,sz) {}
  };

/*! One-dimensional array type, with selectable storage alignment. */
template <typename T, int align>
  class arr_align: public arrT<T,alignAlloc__<T,align> >
  {
  public:
    /*! Creates a zero-sized array. */
    arr_align() : arrT<T,alignAlloc__<T,align> >() {}
    /*! Creates an array with \a sz entries. */
    arr_align(tsize sz) : arrT<T,alignAlloc__<T,align> >(sz) {}
    /*! Creates an array with \a sz entries, and initializes them with
        \a inival. */
    arr_align(tsize sz, const T &inival)
      : arrT<T,alignAlloc__<T,align> >(sz,inival) {}
  };


/*! Two-dimensional array type, with selectable storage management.
    The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by double dereferencing.
    The indices start at zero. */
template <typename T, typename storageManager> class arr2T
  {
  private:
    tsize s1, s2;
    arrT<T, storageManager> d;

  public:
    /*! Creates a zero-sized array. */
    arr2T() : s1(0), s2(0) {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2T(tsize sz1, tsize sz2)
      : s1(sz1), s2(sz2), d(s1*s2) {}
    /*! Creates an array with the dimensions  \a sz1 and \a sz2
        and initializes them with \a inival. */
    arr2T(tsize sz1, tsize sz2, const T &inival)
      : s1(sz1), s2(sz2), d (s1*s2)
      { fill(inival); }
    /*! Creates the array as a copy of \a orig. */
    arr2T(const arr2T &orig)
      : s1(orig.s1), s2(orig.s2), d(orig.d) {}
    /*! Frees the memory associated with the array. */
    ~arr2T() {}

    /*! Returns the first array dimension. */
    tsize size1() const { return s1; }
    /*! Returns the second array dimension. */
    tsize size2() const { return s2; }
    /*! Returns the total array size, i.e. the product of both dimensions. */
    tsize size () const { return s1*s2; }

    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit.
        \a sz1 or \a sz2 can be 0. If \a sz1*sz2 is the same as the
        currently allocated space, no reallocation is performed. */
    void alloc (tsize sz1, tsize sz2)
      {
      if (sz1*sz2 != d.size())
        d.alloc(sz1*sz2);
      s1=sz1; s2=sz2;
      }
    /*! Allocates space for an array with \a sz1*sz2 elements.
        All elements are set to \a inival.
        \a sz1 or \a sz2 can be 0. If \a sz1*sz2 is the same as the
        currently allocated space, no reallocation is performed. */
    void allocAndFill (tsize sz1, tsize sz2, const T &inival)
      { alloc(sz1,sz2); fill(inival); }
    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit.
        \a sz1 or \a sz2 can be 0. If \a sz1*sz2 is smaller than the
        currently allocated space, no reallocation is performed. */
    void fast_alloc (tsize sz1, tsize sz2)
      {
      if (sz1*sz2<=d.size())
        { s1=sz1; s2=sz2; }
      else
        alloc(sz1,sz2);
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); s1=0; s2=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { for (tsize m=0; m<s1*s2; ++m) d[m]=val; }

    /*! Changes the array to be a copy of \a orig. */
    arr2T &operator= (const arr2T &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2);
      d = orig.d;
      return *this;
      }

    /*! Returns a pointer to the beginning of slice \a n. */
    template<typename T2> T *operator[] (T2 n) {return &d[n*s2];}
    /*! Returns a constant pointer to the beginning of slice \a n. */
    template<typename T2> const T *operator[] (T2 n) const {return &d[n*s2];}

    /*! Returns a reference to the element with the indices \a n1 and \a n2. */
    template<typename T2, typename T3> T &operator() (T2 n1, T3 n2)
      {return d[n1*s2 + n2];}
    /*! Returns a constant reference to the element with the indices
        \a n1 and \a n2. */
    template<typename T2, typename T3> const T &operator() (T2 n1, T3 n2) const
      {return d[n1*s2 + n2];}

    /*! Returns the minimum and maximum entry in \a minv and \a maxv,
        respectively. Throws an exception if the array is zero-sized. */
    void minmax (T &minv, T &maxv) const
      {
      planck_assert(s1*s2>0,
        "trying to find min and max of a zero-sized array");
      minv=maxv=d[0];
      for (tsize m=1; m<s1*s2; ++m)
        {
        if (d[m]<minv) minv=d[m];
        if (d[m]>maxv) maxv=d[m];
        }
      }

    /*! Swaps contents and sizes with \a other. */
    void swap (arr2T &other)
      {
      d.swap(other.d);
      std::swap(s1,other.s1);
      std::swap(s2,other.s2);
      }

    /*! Returns \c true if the array and \a other have the same dimensions,
        else \c false. */
    template<typename T2, typename T3> bool conformable
      (const arr2T<T2,T3> &other) const
      { return (other.size1()==s1) && (other.size2()==s2); }
  };

/*! Two-dimensional array type. The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by double dereferencing.
    The indices start at zero. */
template <typename T>
  class arr2: public arr2T<T,normalAlloc__<T> >
  {
  public:
    /*! Creates a zero-sized array. */
    arr2() : arr2T<T,normalAlloc__<T> > () {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2(tsize sz1, tsize sz2) : arr2T<T,normalAlloc__<T> > (sz1,sz2) {}
    /*! Creates an array with the dimensions  \a sz1 and \a sz2
        and initializes them with \a inival. */
    arr2(tsize sz1, tsize sz2, const T &inival)
      : arr2T<T,normalAlloc__<T> > (sz1,sz2,inival) {}
  };

/*! Two-dimensional array type, with selectable storage alignment.
    The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by double dereferencing.
    The indices start at zero. */
template <typename T, int align>
  class arr2_align: public arr2T<T,alignAlloc__<T,align> >
  {
  public:
    /*! Creates a zero-sized array. */
    arr2_align() : arr2T<T,alignAlloc__<T,align> > () {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2_align(tsize sz1, tsize sz2)
      : arr2T<T,alignAlloc__<T,align> > (sz1,sz2) {}
    /*! Creates an array with the dimensions  \a sz1 and \a sz2
        and initializes them with \a inival. */
    arr2_align(tsize sz1, tsize sz2, const T &inival)
      : arr2T<T,alignAlloc__<T,align> > (sz1,sz2,inival) {}
  };

/*! Two-dimensional array type. An entry is located by double dereferencing,
    i.e. via an array of pointers. The indices start at zero. */
template <typename T> class arr2b
  {
  private:
    tsize s1, s2;
    arr<T> d;
    arr<T *> d1;

    void fill_d1()
      { for (tsize m=0; m<s1; ++m) d1[m] = &d[m*s2]; }

  public:
    /*! Creates a zero-sized array. */
    arr2b() : s1(0), s2(0), d(0), d1(0) {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2b(tsize sz1, tsize sz2)
      : s1(sz1), s2(sz2), d(s1*s2), d1(s1)
      { fill_d1(); }
    /*! Creates the array as a copy of \a orig. */
    arr2b(const arr2b &orig)
      : s1(orig.s1), s2(orig.s2), d(orig.d), d1(s1)
      { fill_d1(); }
    /*! Frees the memory associated with the array. */
    ~arr2b() {}

    /*! Returns the first array dimension. */
    tsize size1() const { return s1; }
    /*! Returns the second array dimension. */
    tsize size2() const { return s2; }
    /*! Returns the total array size, i.e. the product of both dimensions. */
    tsize size () const { return s1*s2; }

    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit. */
    void alloc (tsize sz1, tsize sz2)
      {
      if ((s1==sz1) && (s2==sz2)) return;
      s1=sz1; s2=sz2;
      d.alloc(s1*s2);
      d1.alloc(s1);
      fill_d1();
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); d1.dealloc(); s1=0; s2=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { d.fill(val); }

    /*! Changes the array to be a copy of \a orig. */
    arr2b &operator= (const arr2b &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2);
      for (tsize m=0; m<s1*s2; ++m) d[m] = orig.d[m];
      return *this;
      }

    /*! Returns a pointer to the beginning of slice \a n. */
    template<typename T2> T *operator[] (T2 n) {return d1[n];}
    /*! Returns a constant pointer to the beginning of slice \a n. */
    template<typename T2> const T *operator[] (T2 n) const {return d1[n];}

    /*! Returns a pointer to the beginning of the pointer array. */
    T **p0() {return &d1[0];}
  };


/*! Three-dimensional array type. The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by multiple dereferencing.
    The indices start at zero. */
template <typename T> class arr3
  {
  private:
    tsize s1, s2, s3, s2s3;
    arr<T> d;

  public:
    /*! Creates a zero-sized array. */
    arr3() : s1(0), s2(0), s3(0), s2s3(0), d(0) {}
    /*! Creates an array with the dimensions \a sz1, \a sz2 and \a sz3. */
    arr3(tsize sz1, tsize sz2, tsize sz3)
      : s1(sz1), s2(sz2), s3(sz3), s2s3(s2*s3), d(s1*s2*s3) {}
    /*! Creates the array as a copy of \a orig. */
    arr3(const arr3 &orig)
      : s1(orig.s1), s2(orig.s2), s3(orig.s3), s2s3(orig.s2s3), d(orig.d) {}
    /*! Frees the memory associated with the array. */
    ~arr3() {}

    /*! Returns the first array dimension. */
    tsize size1() const { return s1; }
    /*! Returns the second array dimension. */
    tsize size2() const { return s2; }
    /*! Returns the third array dimension. */
    tsize size3() const { return s3; }
    /*! Returns the total array size, i.e. the product of all dimensions. */
    tsize size () const { return s1*s2*s3; }

    /*! Allocates space for an array with \a sz1*sz2*sz3 elements.
        The content of the array is undefined on exit. */
    void alloc (tsize sz1, tsize sz2, tsize sz3)
      {
      d.alloc(sz1*sz2*sz3);
      s1=sz1; s2=sz2; s3=sz3; s2s3=s2*s3;
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); s1=0; s2=0; s3=0; s2s3=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { d.fill(val); }

    /*! Changes the array to be a copy of \a orig. */
    arr3 &operator= (const arr3 &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2, orig.s3);
      d = orig.d;
      return *this;
      }

    /*! Returns a reference to the element with the indices
        \a n1, \a n2 and \a n3. */
    template<typename T2, typename T3, typename T4> T &operator()
      (T2 n1, T3 n2, T4 n3)
      {return d[n1*s2s3 + n2*s3 + n3];}
    /*! Returns a constant reference to the element with the indices
        \a n1, \a n2 and \a n3. */
    template<typename T2, typename T3, typename T4> const T &operator()
      (T2 n1, T3 n2, T4 n3) const
      {return d[n1*s2s3 + n2*s3 + n3];}

    /*! Swaps contents and sizes with \a other. */
    void swap (arr3 &other)
      {
      d.swap(other.d);
      std::swap(s1,other.s1);
      std::swap(s2,other.s2);
      std::swap(s3,other.s3);
      std::swap(s2s3,other.s2s3);
      }

    /*! Returns \c true if the array and \a other have the same dimensions,
        else \c false. */
    template<typename T2> bool conformable (const arr3<T2> &other) const
      { return (other.size1()==s1)&&(other.size2()==s2)&&(other.size3()==s3); }
  };

/*! \} */

#endif

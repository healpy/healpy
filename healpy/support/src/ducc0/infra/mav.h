/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*! \file ducc0/infra/mav.h
 *  Classes for dealing with multidimensional arrays
 *
 *  \copyright Copyright (C) 2019-2021 Max-Planck-Society
 *  \author Martin Reinecke
 *  */

#ifndef DUCC0_MAV_H
#define DUCC0_MAV_H

#include <array>
#include <vector>
#include <memory>
#include <numeric>
#include <cstddef>
#include <functional>
#include <tuple>
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/aligned_array.h"
#include "ducc0/infra/misc_utils.h"

namespace ducc0 {

namespace detail_mav {

using namespace std;

struct uninitialized_dummy {};
constexpr uninitialized_dummy UNINITIALIZED;

template<typename T> class membuf
  {
  protected:
    shared_ptr<vector<T>> ptr;
    shared_ptr<aligned_array<T>> rawptr;
    const T *d;
    bool rw;

    membuf(const T *d_, membuf &other)
      : ptr(other.ptr), rawptr(other.rawptr), d(d_), rw(other.rw) {}
    membuf(const T *d_, const membuf &other)
      : ptr(other.ptr), rawptr(other.rawptr), d(d_), rw(false) {}

    // externally owned data pointer
    membuf(T *d_, bool rw_=false)
      : d(d_), rw(rw_) {}
    // externally owned data pointer, nonmodifiable
    membuf(const T *d_)
      : d(d_), rw(false) {}
    // share another memory buffer, but read-only
    membuf(const membuf &other)
      : ptr(other.ptr), d(other.d), rw(false) {}
#if defined(_MSC_VER)
    // MSVC is broken
    membuf(membuf &other)
      : ptr(other.ptr), d(other.d), rw(other.rw) {}
    membuf(membuf &&other)
      : ptr(move(other.ptr)), d(move(other.d)), rw(move(other.rw)) {}
#else
    // share another memory buffer, using the same read/write permissions
    membuf(membuf &other) = default;
    // take over another memory buffer
    membuf(membuf &&other) = default;
#endif

  public:
    // allocate own memory
    membuf() : d(nullptr), rw(false) {}
    membuf(size_t sz)
      : ptr(make_shared<vector<T>>(sz)), d(ptr->data()), rw(true) {}
    membuf(size_t sz, uninitialized_dummy)
      : rawptr(make_shared<aligned_array<T>>(sz)), d(rawptr->data()), rw(true) {}
    void assign(membuf &other)
      {
      ptr = other.ptr;
      rawptr = other.rawptr;
      d = other.d;
      rw = other.rw;
      }
    void assign(const membuf &other)
      {
      ptr = other.ptr;
      rawptr = other.rawptr;
      d = other.d;
      rw = false;
      }
    // read/write access to element #i
    template<typename I> T &vraw(I i)
      {
      MR_assert(rw, "array is not writable");
      return const_cast<T *>(d)[i];
      }
    // read access to element #i
    template<typename I> const T &craw(I i) const
      { return d[i]; }
    // read/write access to data area
    const T *cdata() const
      { return d; }
    // read access to data area
    T *vdata()
      {
      MR_assert(rw, "array is not writable");
      return const_cast<T *>(d);
      }
    bool writable() const { return rw; }
  };

constexpr size_t MAXIDX=~(size_t(0));

/// Helper class containing shape and stride information of an `fmav` object
class fmav_info
  {
  public:
    /// vector of nonnegative integers for storing the array shape
    using shape_t = vector<size_t>;
    /// vector of integers for storing the array strides
    using stride_t = vector<ptrdiff_t>;

  protected:
    shape_t shp;
    stride_t str;
    size_t sz;

    static stride_t shape2stride(const shape_t &shp)
      {
      auto ndim = shp.size();
      stride_t res(ndim);
      if (ndim==0) return res;
      res[ndim-1]=1;
      for (size_t i=2; i<=ndim; ++i)
        res[ndim-i] = res[ndim-i+1]*ptrdiff_t(shp[ndim-i+1]);
      return res;
      }
    template<typename... Ns> ptrdiff_t getIdx(size_t dim, size_t n, Ns... ns) const
      { return str[dim]*ptrdiff_t(n) + getIdx(dim+1, ns...); }
    ptrdiff_t getIdx(size_t dim, size_t n) const
      { return str[dim]*ptrdiff_t(n); }
    ptrdiff_t getIdx(size_t /*dim*/) const
      { return 0; }

  public:
    /// Constructs a 1D object with all extents and strides set to zero.
    fmav_info() : shp(1,0), str(1,0), sz(0) {}
    /// Constructs an object with the given shape and stride.
    fmav_info(const shape_t &shape_, const stride_t &stride_)
      : shp(shape_), str(stride_), sz(accumulate(shp.begin(),shp.end(),size_t(1),multiplies<>()))
      {
      MR_assert(shp.size()==str.size(), "dimensions mismatch");
      }
    /// Constructs an object with the given shape and computes the strides
    /// automatically, assuming a C-contiguous memory layout.
    fmav_info(const shape_t &shape_)
      : fmav_info(shape_, shape2stride(shape_)) {}
    void assign(const fmav_info &other)
      {
      shp = other.shp;
      str = other.str;
      sz = other.sz;
      }
    /// Returns the dimensionality of the object.
    size_t ndim() const { return shp.size(); }
    /// Returns the total number of entries in the object.
    size_t size() const { return sz; }
    /// Returns the shape of the object.
    const shape_t &shape() const { return shp; }
    /// Returns the common broadcast shape of *this and \a shp2
    shape_t bcast_shape(const shape_t &shp2) const
      {
      shape_t res(max(shp.size(), shp2.size()), 1);
      for (size_t i=0; i<shp.size(); ++i)
        res[i+res.size()-shp.size()] = shp[i];
      for (size_t i=0; i<shp2.size(); ++i)
        {
        size_t i2 = i+res.size()-shp2.size();
        if (res[i2]==1)
          res[i2] = shp2[i];
        else
          MR_assert((res[i2]==shp2[i])||(shp2[i]==1),
            "arrays cannot be broadcast together");
        }
      return res;
      }
    void bcast_to_shape(const shape_t &shp2)
      {
      MR_assert(shp2.size()>=shp.size(), "cannot reduce dimensionallity");
      stride_t newstr(shp2.size(), 0);
      for (size_t i=0; i<shp.size(); ++i)
        {
        size_t i2 = i+shp2.size()-shp.size();
        if (shp[i]!=1)
          {
          MR_assert(shp[i]==shp2[i2], "arrays cannot be broadcast together");
          newstr[i2] = str[i];
          }
        }
      shp = shp2;
      str = newstr;
      }
    void prepend_dim()
      {
      shape_t shp2(shp.size()+1);
      stride_t str2(str.size()+1);
      shp2[0] = 1;
      str2[0] = 0;
      for (size_t i=0; i<shp.size(); ++i)
        {
        shp2[i+1] = shp[i];
        str2[i+1] = str[i];
        }
      shp = shp2;
      str = str2;
      }
    /// Returns the length along dimension \a i.
    size_t shape(size_t i) const { return shp[i]; }
    /// Returns the strides of the object.
    const stride_t &stride() const { return str; }
    /// Returns the stride along dimension \a i.
    const ptrdiff_t &stride(size_t i) const { return str[i]; }
    /// Returns true iff the last dimension has stride 1.
    /**  Typically used for optimization purposes. */
    bool last_contiguous() const
      { return ((ndim()==0) || (str.back()==1)); }
    /** Returns true iff the object is C-contiguous, i.e. if the stride of the
     *  last dimension is 1, the stride for the next-to-last dimension is the
     *  shape of the last dimension etc. */
    bool contiguous() const
      {
      auto ndim = shp.size();
      ptrdiff_t stride=1;
      for (size_t i=0; i<ndim; ++i)
        {
        if (str[ndim-1-i]!=stride) return false;
        stride *= ptrdiff_t(shp[ndim-1-i]);
        }
      return true;
      }
    /// Returns true iff this->shape and \a other.shape match.
    bool conformable(const fmav_info &other) const
      { return shp==other.shp; }
    /// Returns the one-dimensional index of an entry from the given
    /// multi-dimensional index tuple, taking strides into account.
    template<typename... Ns> ptrdiff_t idx(Ns... ns) const
      {
      MR_assert(ndim()==sizeof...(ns), "incorrect number of indices");
      return getIdx(0, ns...);
      }
  };

/// Helper class containing shape and stride information of a `mav` object
template<size_t ndim> class mav_info
  {
  public:
    /// Fixed-size array of nonnegative integers for storing the array shape
    using shape_t = array<size_t, ndim>;
    /// Fixed-size array of integers for storing the array strides
    using stride_t = array<ptrdiff_t, ndim>;

  protected:
    shape_t shp;
    stride_t str;
    size_t sz;

    static stride_t shape2stride(const shape_t &shp)
      {
      stride_t res;
      if (ndim==0) return res;
      res[ndim-1]=1;
      for (size_t i=2; i<=ndim; ++i)
        res[ndim-i] = res[ndim-i+1]*ptrdiff_t(shp[ndim-i+1]);
      return res;
      }
    template<typename... Ns> ptrdiff_t getIdx(size_t dim, size_t n, Ns... ns) const
      { return str[dim]*n + getIdx(dim+1, ns...); }
    ptrdiff_t getIdx(size_t dim, size_t n) const
      { return str[dim]*n; }
    ptrdiff_t getIdx(size_t /*dim*/) const
      { return 0; }

  public:
    /// Constructs an object with all extents and strides set to zero.
    mav_info() : sz(0)
      {
      for (size_t i=0; i<ndim; ++i)
        { shp[i]=0; str[i]=0; }
      }
    /// Constructs an object with the given shape and stride.
    mav_info(const shape_t &shape_, const stride_t &stride_)
      : shp(shape_), str(stride_), sz(accumulate(shp.begin(),shp.end(),size_t(1),multiplies<>())) {}
    /// Constructs an object with the given shape and computes the strides
    /// automatically, assuming a C-contiguous memory layout.
    mav_info(const shape_t &shape_)
      : mav_info(shape_, shape2stride(shape_)) {}
    void assign(const mav_info &other)
      {
      shp = other.shp;
      str = other.str;
      sz = other.sz;
      }
    /// Returns the total number of entries in the object.
    size_t size() const { return sz; }
    /// Returns the shape of the object.
    const shape_t &shape() const { return shp; }
    /// Returns the length along dimension \a i.
    size_t shape(size_t i) const { return shp[i]; }
    /// Returns the strides of the object.
    const stride_t &stride() const { return str; }
    /// Returns the stride along dimension \a i.
    const ptrdiff_t &stride(size_t i) const { return str[i]; }
    /// Returns true iff the last dimension has stride 1.
    /**  Typically used for optimization purposes. */
    bool last_contiguous() const
      { return ((ndim==0) || (str.back()==1)); }
    /** Returns true iff the object is C-contiguous, i.e. if the stride of the
     *  last dimension is 1, the stride for the next-to-last dimension is the
     *  shape of the last dimension etc. */
    bool contiguous() const
      {
      ptrdiff_t stride=1;
      for (size_t i=0; i<ndim; ++i)
        {
        if (str[ndim-1-i]!=stride) return false;
        stride *= ptrdiff_t(shp[ndim-1-i]);
        }
      return true;
      }
    /// Returns true iff this->shape and \a other.shape match.
    bool conformable(const mav_info &other) const
      { return shp==other.shp; }
    /// Returns true iff this->shape and \a other match.
    bool conformable(const shape_t &other) const
      { return shp==other; }
    /// Returns the one-dimensional index of an entry from the given
    /// multi-dimensional index tuple, taking strides into account.
    template<typename... Ns> ptrdiff_t idx(Ns... ns) const
      {
      static_assert(ndim==sizeof...(ns), "incorrect number of indices");
      return getIdx(0, ns...);
      }
  };


class FmavIter
  {
  private:
    fmav_info::shape_t pos;
    fmav_info arr;
    ptrdiff_t p;
    size_t rem;

  public:
    FmavIter(const fmav_info &arr_)
      : pos(arr_.ndim(), 0), arr(arr_), p(0), rem(arr_.size()) {}
    void advance()
      {
      --rem;
      for (int i_=int(pos.size())-1; i_>=0; --i_)
        {
        auto i = size_t(i_);
        p += arr.stride(i);
        if (++pos[i] < arr.shape(i))
          return;
        pos[i] = 0;
        p -= ptrdiff_t(arr.shape(i))*arr.stride(i);
        }
      }
    ptrdiff_t ofs() const { return p; }
    size_t remaining() const { return rem; }
  };


/// Class for storing (or referring to) multi-dimensional arrays with a
/// dimensionality that is not known at compile time.
/** "fmav" stands for "flexible multidimensional array view".
 *  The shape must consist of non-negative integers (zeros are allowed).
 *  Strides may be positive or negative; stride values of zero are accepted and
 *  may be useful in specific circumstances (e.g. read-only arrays with the same
 *  value everywhere.
 *
 *  An fmav may "own" or "not own" the memory holding its array data. If it does
 *  not own the memory, it will not be deallocated when the mav is destroyed.
 *  If it owns the memory, this "ownership" may be shared with other fmav objects.
 *  Memory is only deallocated if the last fmav object owning it is destroyed. */
template<typename T> class fmav: public fmav_info, public membuf<T>
  {
  protected:
    using tbuf = membuf<T>;
    using tinfo = fmav_info;

  public:
    using typename tinfo::shape_t;
    using typename tinfo::stride_t;

  protected:
    template<typename Func, typename T2> void applyHelper(size_t idim,
      ptrdiff_t idx, ptrdiff_t idx2, const fmav<T2> &other, Func func)
      {
      auto ndim = tinfo::ndim();
      if (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<Func>(idim+1, idx+i*str[idim], idx2+i*other.stride(idim), other, func);
      else
        {
        T *d1 = vdata();
        const T2 *d2 = other.cdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d1[idx+i*str[idim]], d2[idx2+i*other.stride(idim)]);
        }
      }
    template<typename Func, typename T2> void applyHelper(size_t idim,
      ptrdiff_t idx, ptrdiff_t idx2, const fmav<T2> &other, Func func) const
      {
      auto ndim = tinfo::ndim();
      if (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<Func>(idim+1, idx+i*str[idim], idx2+i*other.stride(idim), other, func);
      else
        {
        const T *d1 = cdata();
        const T2 *d2 = other.cdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d1[idx+i*str[idim]], d2[idx2+i*other.stride(idim)]);
        }
      }
    template<typename Func> void applyHelper(size_t idim, ptrdiff_t idx, Func func)
      {
      auto ndim = tinfo::ndim();
      if (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<Func>(idim+1, idx+i*str[idim], func);
      else
        {
        T *d2 = vdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d2[idx+i*str[idim]]);
        }
      }
    template<typename Func> void applyHelper(size_t idim, ptrdiff_t idx, Func func) const
      {
      auto ndim = tinfo::ndim();
      if (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<Func>(idim+1, idx+i*str[idim], func);
      else
        {
        const T *d2 = cdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d2[idx+i*str[idim]]);
        }
      }

    auto subdata(const shape_t &i0, const shape_t &extent) const
      {
      auto ndim = tinfo::ndim();
      shape_t nshp(ndim);
      stride_t nstr(ndim);
      ptrdiff_t nofs;
      MR_assert(i0.size()==ndim, "bad dimensionality");
      MR_assert(extent.size()==ndim, "bad dimensionality");
      size_t n0=0;
      for (auto x:extent) if (x==0) ++n0;
      nofs=0;
      nshp.resize(ndim-n0);
      nstr.resize(ndim-n0);
      for (size_t i=0, i2=0; i<ndim; ++i)
        {
        MR_assert(i0[i]<shp[i], "bad subset");
        nofs+=i0[i]*str[i];
        if (extent[i]!=0)
          {
          auto ext = extent[i];
          if (ext==MAXIDX)
            ext = shp[i]-i0[i];
          MR_assert(i0[i]+ext<=shp[i], "bad subset");
          nshp[i2]=ext; nstr[i2]=str[i];
          ++i2;
          }
        }
      return make_tuple(nshp, nstr, nofs);
      }

  public:
    using tbuf::vraw, tbuf::craw, tbuf::vdata, tbuf::cdata;
    /// Constructs a 1D fmav with size and stride zero and no data content.
    fmav() {}
    /** Constructs a read-only fmav with its first data entry at \a d
     *  and the given shape and strides. The fmav does not own the memory. */
    fmav(const T *d_, const shape_t &shp_, const stride_t &str_)
      : tinfo(shp_, str_), tbuf(d_) {}
    /** Constructs a read-only fmav with its first data entry at \a d
     *  and the given shape. The array is assumed to be C-contiguous.
     *  The fmav does not own the memory. */
    fmav(const T *d_, const shape_t &shp_)
      : tinfo(shp_), tbuf(d_) {}
    /** Constructs an fmav with its first data entry at \a d
     *  and the given shape and strides. The fmav does not own the memory.
     *  Iff \a rw_ is true, write accesses to the array are allowed. */
    fmav(T *d_, const shape_t &shp_, const stride_t &str_, bool rw_)
      : tinfo(shp_, str_), tbuf(d_,rw_) {}
    /** Constructs an fmav with its first data entry at \a d and the given shape.
     *  The array is assumed to be C-contiguous.
     *  The fmav does not own the memory.
     *  Iff \a rw_ is true, write accesses to the array are allowed. */
    fmav(T *d_, const shape_t &shp_, bool rw_)
      : tinfo(shp_), tbuf(d_,rw_) {}
    /** Constructs a C-contiguous read/write fmav with the given shape.
     *  The array contents are default-initialized.
     *  The fmav owns the array memory. */
    fmav(const shape_t &shp_)
      : tinfo(shp_), tbuf(size()) {}
    /** Constructs a C-contiguous read/write fmav with the given shape.
     *  The array contents are not initialized.
     *  The fmav owns the array memory. */
    fmav(const shape_t &shp_, uninitialized_dummy)
      : tinfo(shp_), tbuf(size(), UNINITIALIZED) {}
    fmav(const shape_t &shp_, const stride_t &str_, uninitialized_dummy)
      : tinfo(shp_, str_), tbuf(size(), UNINITIALIZED)
      {
      ptrdiff_t ofs=0;
      for (size_t i=0; i<ndim(); ++i)
        ofs += (ptrdiff_t(shp[i])-1)*str[i];
      MR_assert(ofs+1==ptrdiff_t(size()), "array is not compact");
      }
    fmav(const shape_t &shp_, const stride_t &str_)
      : tinfo(shp_, str_), tbuf(size())
      {
      ptrdiff_t ofs=0;
      for (size_t i=0; i<ndim(); ++i)
        ofs += (ptrdiff_t(shp[i])-1)*str[i];
      MR_assert(ofs+1==ptrdiff_t(size()), "array is not compact");
      }
    fmav(const T* d_, const tinfo &info)
      : tinfo(info), tbuf(d_) {}
    fmav(T* d_, const tinfo &info, bool rw_=false)
      : tinfo(info), tbuf(d_, rw_) {}
#if defined(_MSC_VER)
    // MSVC is broken
    fmav(const fmav &other) : tinfo(other), tbuf(other) {}
    fmav(fmav &other) : tinfo(other), tbuf(other) {}
    fmav(fmav &&other) : tinfo(other), tbuf(other) {}
#else
    /** Constructs a read-only fmav with the same shape and strides as \a other,
     *  pointing to the same memory. Ownership is shared. */
    fmav(const fmav &other) = default;
    /** Constructs an fmav with the same read-write status, shape and strides
     *  as \a other, pointing to the same memory. Ownership is shared. */
    fmav(fmav &other) = default;
    fmav(fmav &&other) = default;
#endif
    fmav(tbuf &buf, const shape_t &shp_, const stride_t &str_)
      : tinfo(shp_, str_), tbuf(buf) {}
    fmav(const tbuf &buf, const shape_t &shp_, const stride_t &str_)
      : tinfo(shp_, str_), tbuf(buf) {}
    fmav(const shape_t &shp_, const stride_t &str_, const T *d_, tbuf &buf)
      : tinfo(shp_, str_), tbuf(d_, buf) {}
    fmav(const shape_t &shp_, const stride_t &str_, const T *d_, const tbuf &buf)
      : tinfo(shp_, str_), tbuf(d_, buf) {}

    void assign(fmav &other)
      {
      fmav_info::assign(other);
      membuf<T>::assign(other);
      }
    void assign(const fmav &other)
      {
      fmav_info::assign(other);
      membuf<T>::assign(other);
      }

    /// Returns the data entry at the given set of indices.
    template<typename... Ns> const T &operator()(Ns... ns) const
      { return craw(idx(ns...)); }
    /// Returns the data entry at the given set of indices.
    template<typename... Ns> const T &c(Ns... ns) const
      { return craw(idx(ns...)); }
    /** Returns a writable reference to the data entry at the given set of
     *  indices. This call will throw an exception if the fmav is read-only. */
    template<typename... Ns> T &v(Ns... ns)
      { return vraw(idx(ns...)); }

    /** Returns an fmav (of the same or smaller dimensionality) representing a
     *  sub-array of *this. \a i0 indicates the starting indices, and \a extent
     *  the number of entries along this dimension. If any extent is 0, this
     *  dimension will be omitted in the output array.
     *  Specifying an extent of MAXIDX will make the extent as large as possible.
     *  if *this is writable, the returned fmav will also be writable. */
    fmav subarray(const shape_t &i0, const shape_t &extent)
      {
      auto [nshp, nstr, nofs] = subdata(i0, extent);
      return fmav(nshp, nstr, tbuf::d+nofs, *this);
      }
    /** Returns an fmav (of the same or smaller dimensionality) representing a
     *  sub-array of *this. \a i0 indicates the starting indices, and \a extent
     *  the number of entries along this dimension. If any extent is 0, this
     *  dimension will be omitted in the output array.
     *  Specifying an extent of MAXIDX will make the extent as large as possible.
     *  The returned fmav is read-only. */
    fmav subarray(const shape_t &i0, const shape_t &extent) const
      {
      auto [nshp, nstr, nofs] = subdata(i0, extent);
      return fmav(nshp, nstr, tbuf::d+nofs, *this);
      }
    /** Calls \a func for every entry in the array, passing a reference to it. */
    template<typename Func> void apply(Func func)
      {
      if (contiguous()) // covers 0-d case
        {
        T *d2 = vdata();
        for (auto v=d2; v!=d2+size(); ++v)
          func(*v);
        return;
        }
      applyHelper<Func>(0, 0, func);
      }
    /** Calls \a func for every entry in the array, passing a constant
     *  reference to it. */
    template<typename Func> void apply(Func func) const
      {
      if (contiguous()) // covers 0-d case
        {
        const T *d2 = cdata();
        for (auto v=d2; v!=d2+size(); ++v)
          func(*v);
        return;
        }
      applyHelper<Func>(0, 0, func);
      }
    /** Calls \a func for every entry in the array and the corresponding entry
     *  in \a other, passing constant references. */
    template<typename Func, typename T2> void apply(const fmav<T2> &other, Func func)
      {
      MR_assert(conformable(other), "fmavs are not conformable");
      if (ndim() == 0)
        func(cdata(), other.cdata());
      else
        applyHelper<Func>(0, 0, 0, other, func);
      }
    /** Calls \a func for every entry in the array and the corresponding entry
     *  in \a other, passing a nonconstant reference to the entry in this array
     *  and a constant one for the entry in \a other. */
    template<typename Func, typename T2> void apply(const fmav<T2> &other, Func func) const
      {
      MR_assert(conformable(other), "fmavs are not conformable");
      if (ndim() == 0)
        func(vdata(), other.cdata());
      else
        applyHelper<Func>(0, 0, 0, other, func);
      }
    vector<T> dump() const
      {
      vector<T> res(sz);
      size_t ii=0;
      apply([&](const T&v){res[ii++]=v;});
      return res;
      }
    void load (const vector<T> &v)
      {
      MR_assert(v.size()==sz, "bad input data size");
      size_t ii=0;
      apply([&](T &val){val=v[ii++];});
      }
  };

template<typename T> fmav<T> subarray
  (fmav<T> &arr, const typename fmav<T>::shape_t &i0, const typename fmav<T>::shape_t &extent)  
  { return arr.subarray(i0, extent); }

template<typename T> fmav<T> subarray
  (const fmav<T> &arr, const typename fmav<T>::shape_t &i0, const typename fmav<T>::shape_t &extent)  
  { return arr.subarray(i0, extent); }


// template<typename Func, typename T0, typename Ts...> void fmav_pointwise_op(Func func, T0 & arg0, Ts&... args)
//   {
//   MR_assert(multiequal(arg0.shape()==args.shape()...), "fmav shape mismatch");
//   if (multiequal(true, arg0.stride()==args.stride()...)) // equal strides, we can make simplifications
//     {
//     if (arg0.compact()) // even better, we can go through everything in a single loop
//       {
//       for (size_t i=0; i<arg0.size(); ++i)
//         func(arg0.ptr[i], args.ptr[i]...);
//       }
//     else
//   }

/// Class for storing (or referring to) multi-dimensional arrays with a
/// dimensionality known at compile time.
/** "mav" stands for "multidimensional array view".
 *  The shape must consist of non-negative integers (zeros are allowed).
 *  Strides may be positive or negative; stride values of zero are accepted and
 *  may be useful in specific circumstances (e.g. read-only arrays with the same
 *  value everywhere.
 *
 *  A mav may "own" or "not own" the memory holding its array data. If it does
 *  not own the memory, it will not be deallocated when the mav is destroyed.
 *  If it owns the memory, this "ownership" may be shared with other mav objects.
 *  Memory is only deallocated if the last mav object owning it is destroyed. */
template<typename T, size_t ndim> class mav: public mav_info<ndim>, public membuf<T>
  {
  protected:
    using tinfo = mav_info<ndim>;
    using tbuf = membuf<T>;
    using tinfo::shp, tinfo::str;

  public:
    using typename tinfo::shape_t;
    using typename tinfo::stride_t;

  protected:
    template<size_t idim, typename Func> void applyHelper(ptrdiff_t idx, Func func)
      {
      if constexpr (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<idim+1, Func>(idx+i*str[idim], func);
      else
        {
        T *d2 = vdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d2[idx+i*str[idim]]);
        }
      }
    template<size_t idim, typename Func> void applyHelper(ptrdiff_t idx, Func func) const
      {
      if constexpr (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<idim+1, Func>(idx+i*str[idim], func);
      else
        {
        const T *d2 = cdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d2[idx+i*str[idim]]);
        }
      }
    template<size_t idim, typename T2, typename Func>
      void applyHelper(ptrdiff_t idx, ptrdiff_t idx2,
                       const mav<T2,ndim> &other, Func func)
      {
      if constexpr (idim==0)
        MR_assert(conformable(other), "dimension mismatch");
      if constexpr (idim+1<ndim)
        for (size_t i=0; i<shp[idim]; ++i)
          applyHelper<idim+1, T2, Func>(idx+i*str[idim],
                                        idx2+i*other.stride(idim), other, func);
      else
        {
        T *d2 = vdata();
        const T2 *d3 = other.cdata();
        for (size_t i=0; i<shp[idim]; ++i)
          func(d2[idx+i*str[idim]],d3[idx2+i*other.stride(idim)]);
        }
      }

    template<size_t nd2> auto subdata(const shape_t &i0, const shape_t &extent) const
      {
      array<size_t, nd2> nshp;
      array<ptrdiff_t, nd2> nstr;
      size_t n0=0;
      for (auto x:extent) if (x==0) ++n0;
      MR_assert(n0+nd2==ndim, "bad extent");
      ptrdiff_t nofs=0;
      for (size_t i=0, i2=0; i<ndim; ++i)
        {
        MR_assert(i0[i]<shp[i], "bad subset");
        nofs+=i0[i]*str[i];
        if (extent[i]!=0)
          {
          auto ext = extent[i];
          if (ext==MAXIDX)
            ext = shp[i]-i0[i];
          MR_assert(i0[i]+ext<=shp[i], "bad subset");
          nshp[i2]=ext; nstr[i2]=str[i];
          ++i2;
          }
        }
      return make_tuple(nshp, nstr, nofs);
      }

  public:
    using tbuf::vraw, tbuf::craw, tbuf::vdata, tbuf::cdata;
    using tinfo::contiguous, tinfo::size, tinfo::idx, tinfo::conformable;

    /// Constructs a mav with size and stride zero in all dimensions and no
    /// data content.
    mav() {}
    /** Constructs a read-only mav with its first data entry at \a d
     *  and the given shape and strides. The mav does not own the memory. */
    mav(const T *d_, const shape_t &shp_, const stride_t &str_)
      : tinfo(shp_, str_), tbuf(d_) {}
    /** Constructs a mav with its first data entry at \a d
     *  and the given shape and strides. The mav does not own the memory.
     *  Iff \a rw_ is true, write accesses to the array are allowed. */
    mav(T *d_, const shape_t &shp_, const stride_t &str_, bool rw_=false)
      : tinfo(shp_, str_), tbuf(d_, rw_) {}
    /** Constructs a read-only mav with its first data entry at \a d
     *  and the given shape. The array is assumed to be C-contiguous.
     *  The mav does not own the memory. */
    mav(const T *d_, const shape_t &shp_)
      : tinfo(shp_), tbuf(d_) {}
    /** Constructs a mav with its first data entry at \a d and the given shape.
     *  The array is assumed to be C-contiguous.
     *  The mav does not own the memory.
     *  Iff \a rw_ is true, write accesses to the array are allowed. */
    mav(T *d_, const shape_t &shp_, bool rw_=false)
      : tinfo(shp_), tbuf(d_, rw_) {}
    /** Constructs a C-contiguous read/write mav with the given shape.
     *  The array contents are default-initialized.
     *  The mav owns the array memory. */
    mav(const shape_t &shp_)
      : tinfo(shp_), tbuf(size()) {}
    /** Constructs a C-contiguous read/write mav with the given shape.
     *  The array contents are not initialized.
     *  The mav owns the array memory. */
    mav(const shape_t &shp_, uninitialized_dummy)
      : tinfo(shp_), tbuf(size(), UNINITIALIZED) {}
#if defined(_MSC_VER)
    // MSVC is broken
    mav(const mav &other) : tinfo(other), tbuf(other) {}
    mav(mav &other): tinfo(other), tbuf(other) {}
    mav(mav &&other): tinfo(other), tbuf(other) {}
#else
    /** Constructs a read-only mav with the same shape and strides as \a other,
     *  pointing to the same memory. Ownership is shared. */
    mav(const mav &other) = default;
    /** Constructs a mav with the same read-write status, shape and strides
     *  as \a other, pointing to the same memory. Ownership is shared. */
    mav(mav &other) = default;
    mav(mav &&other) = default;
#endif
    void assign(mav &other)
      {
      mav_info<ndim>::assign(other);
      membuf<T>::assign(other);
      }
    void assign(const mav &other)
      {
      mav_info<ndim>::assign(other);
      membuf<T>::assign(other);
      }
    mav(const shape_t &shp_, const stride_t &str_, const T *d_, membuf<T> &mb)
      : mav_info<ndim>(shp_, str_), membuf<T>(d_, mb) {}
    mav(const shape_t &shp_, const stride_t &str_, const T *d_, const membuf<T> &mb)
      : mav_info<ndim>(shp_, str_), membuf<T>(d_, mb) {}
    operator fmav<T>() const
      {
      return fmav<T>(*this, {shp.begin(), shp.end()}, {str.begin(), str.end()});
      }
    operator fmav<T>()
      {
      return fmav<T>(*this, {shp.begin(), shp.end()}, {str.begin(), str.end()});
      }
    /// Returns the data entry at the given set of indices.
    template<typename... Ns> const T &operator()(Ns... ns) const
      { return craw(idx(ns...)); }
    /// Returns the data entry at the given set of indices.
    template<typename... Ns> const T &c(Ns... ns) const
      { return craw(idx(ns...)); }
    /** Returns a writable reference to the data entry at the given set of
     *  indices. This call will throw an exception if the mav is read-only. */
    template<typename... Ns> T &v(Ns... ns)
      { return vraw(idx(ns...)); }
    /** Calls \a func for every entry in the array, passing a reference to it. */
    template<typename Func> void apply(Func func)
      {
      if (contiguous()) // covers 0-d case
        {
        T *d2 = vdata();
        for (auto v=d2; v!=d2+size(); ++v)
          func(*v);
        return;
        }
      applyHelper<0,Func>(0,func);
      }
    /** Calls \a func for every entry in the array, passing a constant
     *  reference to it. */
    template<typename Func> void apply(Func func) const
      {
      if (contiguous()) // covers 0-d case
        {
        const T *d2 = cdata();
        for (auto v=d2; v!=d2+size(); ++v)
          func(*v);
        return;
        }
      applyHelper<0,Func>(0,func);
      }
    /** Calls \a func for every entry in the array and the corresponding entry
     *  in \a other, passing a nonconstant reference to the entry in this array
     *  and a constant one for the entry in \a other. */
    template<typename T2, typename Func> void apply
      (const mav<T2, ndim> &other,Func func)
      {
      if constexpr (ndim==0)
        func(vdata(), other.cdata());
      else
        applyHelper<0,T2,Func>(0,0,other,func);
      }
    /// Sets every entry of the array to \a val.
    void fill(const T &val)
      { apply([val](T &v){v=val;}); }
    /** Returns a mav (of the same or smaller dimensionality) representing a
     *  sub-array of *this. \a i0 indicates the starting indices, and \a extent
     *  the number of entries along this dimension. If any extent is 0, this
     *  dimension will be omitted in the output array.
     *  Specifying an extent of MAXIDX will make the extent as large as possible.
     *  if *this is writable, the returned mav will also be writable. */
    template<size_t nd2> mav<T,nd2> subarray(const shape_t &i0, const shape_t &extent)
      {
      auto [nshp, nstr, nofs] = subdata<nd2> (i0, extent);
      return mav<T,nd2> (nshp, nstr, tbuf::d+nofs, *this);
      }
    /** Returns a mav (of the same or smaller dimensionality) representing a
     *  sub-array of *this. \a i0 indicates the starting indices, and \a extent
     *  the number of entries along this dimension. If any extent is 0, this
     *  dimension will be omitted in the output array.
     *  Specifying an extent of MAXIDX will make the extent as large as possible.
     *  The returned mav is read-only. */
    template<size_t nd2> mav<T,nd2> subarray(const shape_t &i0, const shape_t &extent) const
      {
      auto [nshp, nstr, nofs] = subdata<nd2> (i0, extent);
      return mav<T,nd2> (nshp, nstr, tbuf::d+nofs, *this);
      }

    /// Returns a zero-extent mav with no associatd data.
    static mav build_empty()
      {
      shape_t nshp;
      nshp.fill(0);
      return mav(static_cast<const T *>(nullptr), nshp);
      }

    /** Returns a read-only mav with the specified shape, filled with \a value.
     *  This is stored as a single value (by using strides of 0) and is
     *  therefore very memory efficient. */
    static mav build_uniform(const shape_t &shape, const T &value)
      {
      membuf<T> buf(1);
      buf.vraw(0) = value;
      stride_t nstr;
      nstr.fill(0);
      return mav(shape, nstr, buf.cdata(), buf);
      }

    /** Returns a writable mav with the specified shape.
     *  The strides are chosen in such a way that critical strides (multiples
     *  of 4096 bytes) along any dimension are avoided, by enlarging the
     *  allocated memory slightly if necessary.
     *  The array data is default-initialized. */
    static mav build_noncritical(const shape_t &shape)
      {
      auto shape2 = noncritical_shape(shape, sizeof(T));
      mav tmp(shape2);
      return tmp.subarray<ndim>(shape_t(), shape);
      }
    /** Returns a writable mav with the specified shape.
     *  The strides are chosen in such a way that critical strides (multiples
     *  of 4096 bytes) along any dimension are avoided, by enlarging the
     *  allocated memory slightly if necessary.
     *  The array data is not initialized. */
    static mav build_noncritical(const shape_t &shape, uninitialized_dummy)
      {
      if (ndim<=1) return mav(shape, UNINITIALIZED);
      auto shape2 = noncritical_shape(shape, sizeof(T));
      mav tmp(shape2, UNINITIALIZED);
      return tmp.subarray<ndim>(shape_t(), shape);
      }
  };

/** Returns a mav (of the same or smaller dimensionality) representing a
 *  sub-array of \a arr. \a i0 indicates the starting indices, and \a extent
 *  the number of entries along this dimension. If any extent is 0, this
 *  dimension will be omitted in the output array.
 *  Specifying an extent of MAXIDX will make the extent as large as possible.
 *  if *thi is writable, the returned mav will also be writable. */
template<size_t nd2, typename T, size_t ndim> mav<T,nd2> subarray
  (mav<T, ndim> &arr, const typename mav<T, ndim>::shape_t &i0, const typename mav<T, ndim>::shape_t &extent)  
  { return arr.template subarray<nd2>(i0, extent); }

/** Returns a mav (of the same or smaller dimensionality) representing a
 *  sub-array of \a arr. \a i0 indicates the starting indices, and \a extent
 *  the number of entries along this dimension. If any extent is 0, this
 *  dimension will be omitted in the output array.
 *  Specifying an extent of MAXIDX will make the extent as large as possible.
 *  The returned mav is read-only. */
template<size_t nd2, typename T, size_t ndim> mav<T,nd2> subarray
  (const mav<T, ndim> &arr, const typename mav<T, ndim>::shape_t &i0, const typename mav<T, ndim>::shape_t &extent)  
  { return arr.template subarray<nd2>(i0, extent); }

template<typename T, size_t ndim> class MavIter
  {
  protected:
    fmav<T> mav;
    array<size_t, ndim> shp;
    array<ptrdiff_t, ndim> str;
    fmav_info::shape_t pos;
    ptrdiff_t idx_;
    bool done_;

    template<typename... Ns> ptrdiff_t getIdx(size_t dim, size_t n, Ns... ns) const
      { return str[dim]*n + getIdx(dim+1, ns...); }
    ptrdiff_t getIdx(size_t dim, size_t n) const
      { return str[dim]*n; }
    ptrdiff_t getIdx(size_t /*dim*/) const
      { return 0; }

  public:
    MavIter(const fmav<T> &mav_)
      : mav(mav_), pos(mav.ndim()-ndim,0), idx_(0), done_(false)
      {
      for (size_t i=0; i<ndim; ++i)
        {
        shp[i] = mav.shape(mav.ndim()-ndim+i);
        str[i] = mav.stride(mav.ndim()-ndim+i);
        }
      }
    MavIter(fmav<T> &mav_)
      : mav(mav_), pos(mav.ndim()-ndim,0), idx_(0), done_(false)
      {
      for (size_t i=0; i<ndim; ++i)
        {
        shp[i] = mav.shape(mav.ndim()-ndim+i);
        str[i] = mav.stride(mav.ndim()-ndim+i);
        }
      }
    bool done() const
      { return done_; }
    void inc()
      {
      for (ptrdiff_t i=mav.ndim()-ndim-1; i>=0; --i)
        {
        idx_+=mav.stride(i);
        if (++pos[i]<mav.shape(i)) return;
        pos[i]=0;
        idx_-=mav.shape(i)*mav.stride(i);
        }
      done_=true;
      }
    size_t shape(size_t i) const { return shp[i]; }
    template<typename... Ns> ptrdiff_t idx(Ns... ns) const
      {
      static_assert(ndim==sizeof...(ns), "incorrect number of indices");
      return idx_ + getIdx(0, ns...);
      }
    template<typename... Ns> const T &operator()(Ns... ns) const
      { return mav.craw(idx(ns...)); }
    template<typename... Ns> const T &c(Ns... ns) const
      { return mav.craw(idx(ns...)); }
    template<typename... Ns> T &v(Ns... ns)
      { return mav.vraw(idx(ns...)); }
  };

}

using detail_mav::UNINITIALIZED;
using detail_mav::fmav_info;
using detail_mav::fmav;
using detail_mav::mav_info;
using detail_mav::mav;
using detail_mav::FmavIter;
using detail_mav::MavIter;
using detail_mav::MAXIDX;
using detail_mav::subarray;

}

#endif

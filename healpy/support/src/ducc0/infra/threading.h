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

/** \file ducc0/infra/threading.h
 *  Mulithreading support, similar to functionality provided by OpenMP
 *
 * \copyright Copyright (C) 2019-2021 Peter Bell, Max-Planck-Society
 * \authors Peter Bell, Martin Reinecke
 */

#ifndef DUCC0_THREADING_H
#define DUCC0_THREADING_H

#include <cstddef>
#include <functional>

namespace ducc0 {

namespace detail_threading {

using std::size_t;

/// Index range describing a chunk of work inside a parallellized loop
struct Range
  {
  size_t lo, //< first index of the chunk
         hi; //< one-past-last index of the chunk
  Range() : lo(0), hi(0) {}
  Range(size_t lo_, size_t hi_) : lo(lo_), hi(hi_) {}
  /// Returns true iff the chunk is not empty
  operator bool() const { return hi>lo; }
  };

/// Class supplied to parallel regions, which allows them to determine their
/// work chunks.
class Scheduler
  {
  public:
    virtual ~Scheduler() {}
    /// Returns the number of threads working in this parallel region
    virtual size_t num_threads() const = 0;
    /// Returns the number of this thread, from the range 0 to num_threads()-1.
    virtual size_t thread_num() const = 0;
    /// Returns information about the next chunk of work.
    /// If this chunk is empty, the work on this thread is done.
    virtual Range getNext() = 0;
  };

/// Returns the maximum number of threads that are supported by the hardware.
/** More threads can be used, but this will probably hurt performance. */
size_t max_threads();
void set_default_nthreads(size_t new_default_nthreads);
size_t get_default_nthreads();

/// Execute \a func over \a nwork work items, on a single thread.
void execSingle(size_t nwork,
  std::function<void(Scheduler &)> func);
/// Execute \a func over \a nwork work items, on \a nthreads threads.
/** Chunks will have the size \a chunksize, except for the last one which
 *  may be smaller.
 * 
 *  Chunks are statically assigned to threads at startup. */
void execStatic(size_t nwork, size_t nthreads, size_t chunksize,
  std::function<void(Scheduler &)> func);
/// Execute \a func over \a nwork work items, on \a nthreads threads.
/** Chunks will have the size \a chunksize, except for the last one which
 *  may be smaller.
 * 
 *  Chunks are assigned dynamically to threads;whenever a thread is finished
 *  with its current chunk, it will obtain the next one from the list of
 *  remaining chunks. */
void execDynamic(size_t nwork, size_t nthreads, size_t chunksize,
  std::function<void(Scheduler &)> func);
void execGuided(size_t nwork, size_t nthreads, size_t chunksize_min,
  double fact_max, std::function<void(Scheduler &)> func);
/// Execute \a func on \a nthreads threads.
/** Work subdivision must be organized within \a func. */
void execParallel(size_t nthreads, std::function<void(Scheduler &)> func);
/// Execute \a func on work items [\a lo; \a hi[ over \a nthreads threads.
/** Work items are subdivided fairly among threads. */
void execParallel(size_t work_lo, size_t work_hi, size_t nthreads,
  std::function<void(size_t, size_t)> func);
/// Execute \a func on work items [0; \a nwork[ over \a nthreads threads.
/** Work items are subdivided fairly among threads. */
inline void execParallel(size_t nwork, size_t nthreads,
  std::function<void(size_t, size_t)> func)
  { execParallel(0, nwork, nthreads, func); }
/// Execute \a func on work items [\a lo; \a hi[ over \a nthreads threads.
/** The first argument to \a func is the thread number.
 *
 *  Work items are subdivided fairly among threads. */
void execParallel(size_t work_lo, size_t work_hi, size_t nthreads,
  std::function<void(size_t, size_t, size_t)> func);
/// Execute \a func on work items [0; \a nwork[ over \a nthreads threads.
/** The first argument to \a func is the thread number.
 *
 *  Work items are subdivided fairly among threads. */
inline void execParallel(size_t nwork, size_t nthreads,
  std::function<void(size_t, size_t, size_t)> func)
  { execParallel(0, nwork, nthreads, func); }

} // end of namespace detail_threading

using detail_threading::max_threads;
using detail_threading::get_default_nthreads;
using detail_threading::set_default_nthreads;
using detail_threading::Scheduler;
using detail_threading::execSingle;
using detail_threading::execStatic;
using detail_threading::execDynamic;
using detail_threading::execGuided;
using detail_threading::execParallel;

} // end of namespace ducc0

#endif

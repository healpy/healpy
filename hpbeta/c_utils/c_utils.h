/*
 *  This file is part of libc_utils.
 *
 *  libc_utils is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libc_utils is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libc_utils; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libc_utils is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file c_utils.h
 *  Convenience functions
 *
 *  Copyright (C) 2008, 2009, 2010 Max-Planck-Society
 *  \author Martin Reinecke
 *  \note This file should only be included from .c files, NOT from .h files.
 */

#ifndef PLANCK_C_UTILS_H
#define PLANCK_C_UTILS_H

#include <math.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void util_fail_ (const char *file, int line, const char *func, const char *msg);
void util_warn_ (const char *file, int line, const char *func, const char *msg);
void *util_malloc_ (size_t sz);
void util_free_ (void *ptr);

#if defined (__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif

#define UTIL_ASSERT(cond,msg) \
  if(!(cond)) util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
#define UTIL_WARN(cond,msg) \
  if(!(cond)) util_warn_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
#define UTIL_FAIL(msg) \
  util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)

#define ALLOC(ptr,type,num) \
  do { (ptr)=(type *)util_malloc_((num)*sizeof(type)); } while (0)
#define RALLOC(type,num) \
  ((type *)util_malloc_((num)*sizeof(type)))
#define DEALLOC(ptr) \
  do { util_free_(ptr); (ptr)=NULL; } while(0)
#define RESIZE(ptr,type,num) \
  do { util_free_(ptr); ALLOC(ptr,type,num); } while(0)
#define REALLOC(ptr,type,num) \
  do { \
    ptr = (type *)realloc(ptr,(num)*sizeof(type)); \
    UTIL_ASSERT(ptr,"realloc() failed"); \
  } while(0)
#define GROW(ptr,type,sz_old,sz_new) \
  do { \
    if ((sz_new)>(sz_old)) \
      { RESIZE(ptr,type,2*(sz_new));sz_old=2*(sz_new); } \
  } while(0)
#define SET_ARRAY(ptr,i1,i2,val) \
  do { \
    ptrdiff_t cnt_; \
    for (cnt_=(i1);cnt_<(i2);++cnt_) (ptr)[cnt_]=(val); \
    } while(0)
#define COPY_ARRAY(src,dest,i1,i2) \
  do { \
    ptrdiff_t cnt_; \
    for (cnt_=(i1);cnt_<(i2);++cnt_) (dest)[cnt_]=(src)[cnt_]; \
    } while(0)

#define ALLOC2D(ptr,type,num1,num2) \
  do { \
    size_t cnt_, num1_=(num1), num2_=(num2); \
    ALLOC(ptr,type *,num1_); \
    ALLOC(ptr[0],type,num1_*num2_); \
    for (cnt_=1; cnt_<num1_; ++cnt_) \
      ptr[cnt_]=ptr[cnt_-1]+num2_; \
    } while(0)
#define DEALLOC2D(ptr) \
  do { if(ptr) DEALLOC((ptr)[0]); DEALLOC(ptr); } while(0)

#define FAPPROX(a,b,eps) \
  (fabs((a)-(b))<((eps)*fabs(b)))
#define ABSAPPROX(a,b,eps) \
  (fabs((a)-(b))<(eps))
#define IMAX(a,b) \
  (((a)>(b)) ? (a) : (b))
#define IMIN(a,b) \
  (((a)<(b)) ? (a) : (b))

#define SWAP(a,b,type) \
  do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)

#define CHECK_STACK_ALIGN(align) \
  do { \
    double foo; \
    UTIL_WARN((((size_t)(&foo))&(align-1))==0, \
      "WARNING: stack not sufficiently aligned!"); \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif

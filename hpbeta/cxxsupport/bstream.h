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

/*! \file bstream.h
 *  Classes for binary I/O with streams
 *
 *  Copyright (C) 2010 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_BSTREAM_H
#define PLANCK_BSTREAM_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include "datatypes.h"

/*! An object of this class can be cast to a \a bool, which indicates whether
    the computer architecture is big-endian (true) or little-endian (false). */
class EndianTest__
  {
  private:
    bool big_end;

  public:
    EndianTest__()
      {
      union { uint16 i16; uint8 i8; } tmp;
      tmp.i16 = 1;
      big_end = (tmp.i8==0);
      }
    operator bool() const { return big_end; }
  };

const EndianTest__ big_endian;

template<size_t size> inline void byteswap_helper__ (char *);
template<> inline void byteswap_helper__<1> (char *)
  {}
template<> inline void byteswap_helper__<2> (char *val)
  {
  using namespace std; 
  swap (val[0],val[1]);
  }
template<> inline void byteswap_helper__<4> (char *val)
  {
  using namespace std;
  swap (val[0],val[3]); swap (val[1],val[2]);
  }
template<> inline void byteswap_helper__<8> (char *val)
  {
  using namespace std;
  swap (val[0],val[7]); swap (val[1],val[6]);
  swap (val[2],val[5]); swap (val[3],val[4]);
  }

/*! Performs an endianness conversion on \a val.
    \note \a T must be a primitive data type! */
template<typename T> inline void byteswap (T &val)
  { byteswap_helper__<sizeof(T)> (reinterpret_cast<char *> (&val)); }

const bool file_is_lsb=big_endian, file_is_msb=!big_endian,
           file_is_natural=false;

/*! Class for writing binary data to a stream. */
class bostream
  {
  private:
    std::ostream &s;
    bool doswap;

  public:
    /*! Creates a new object which is attached to \a s_ and performs
        endianness conversion if \a doswap_==true. */
    bostream (std::ostream &s_, bool doswap_=false)
      : s(s_), doswap(doswap_) {}

    /*! Writes a binary representation of \a num objects of type \a T
        (stored in \a data) to the attached stream. Endianness conversion
        is performed if requested in the constructor.
        \note \a T must be a primitive data type! */
    template<typename T> bostream &put (const T *data, size_t num)
      {
      if ((sizeof(T)>1) && doswap)
        for (size_t m=0; m<num; ++m)
          {
          T tmp=data[m];
          byteswap (tmp);
          s.write (reinterpret_cast<const char *> (&tmp), sizeof(T));
          }
      else
        s.write (reinterpret_cast<const char *> (data), num*sizeof(T));
      return *this;
      }
    /*! Writes a binary representation of \a data to the attached stream.
        Endianness conversion is performed if requested in the constructor.
        \note \a T must be a primitive data type! */
    template<typename T> bostream &operator<< (const T &data)
      {
      put(&data,1);
      return *this;
      }

    bool getSwap() const
      { return doswap; }
    void setSwap(bool newswap)
      { doswap=newswap; }
    void flipSwap()
      { doswap=!doswap; }
  };

/*! Class for reading binary data from a stream. */
class bistream
  {
  private:
    std::istream &s;
    bool doswap;

  public:
    /*! Creates a new object which is attached to \a s_ and performs
        endianness conversion if \a doswap_==true. */
    bistream (std::istream &s_, bool doswap_=false)
      : s(s_), doswap(doswap_) {}

    /*! Reads a binary representation of \a num objects of type \a T
        from the attached stream and stores them in \a data. Endianness
        conversion is performed if requested in the constructor.
        \note \a T must be a primitive data type! */
    template<typename T> bistream &get (T *data, size_t num)
      {
      s.read (reinterpret_cast<char *> (data), num*sizeof(T));
      if ((sizeof(T)>1) && doswap)
        for (size_t m=0; m<num; ++m)
          byteswap (data[m]);
      return *this;
      }
    /*! Reads a binary representation of \a data from the attached stream.
        Endianness conversion is performed if requested in the constructor.
        \note \a T must be a primitive data type! */
    template<typename T> bistream &operator>> (T &data)
      {
      get (&data,1);
      return *this;
      }

    bool getSwap() const
      { return doswap; }
    void setSwap(bool newswap)
      { doswap=newswap; }
    void flipSwap()
      { doswap=!doswap; }
  };

class bofstream: public std::ofstream
  {
  private:
    bool doswap;

  public:
    /*! */
    bofstream (const char *fname, bool doswap_)
      : std::ofstream(fname,std::ios::binary), doswap(doswap_) {}

    template<typename T> bofstream &operator<< (const T &data)
      {
      if (doswap)
        {
        T tmp = data;
        byteswap (tmp);
        write (reinterpret_cast<const char *> (&tmp), sizeof(T));
        }
      else
        write (reinterpret_cast<const char *> (&data), sizeof(T));
      return *this;
      }
    template<typename T> bofstream &put (const T *data, size_t num)
      {
      if (doswap)
        for (size_t m=0; m<num; ++m)
          {
          T tmp=data[m];
          byteswap (tmp);
          write (reinterpret_cast<const char *> (&tmp), sizeof(T));
          }
      else
        write (reinterpret_cast<const char *> (data), num*sizeof(T));
      return *this;
      }

    bool getSwap() const
      { return doswap; }
    void setSwap(bool newswap)
      { doswap=newswap; }
    void flipSwap()
      { doswap=!doswap; }
  };

class bifstream: public std::ifstream
  {
  private:
    bool doswap;

  public:
    /*! */
    bifstream ()
      : doswap(false) {}
    bifstream (const char *fname, bool doswap_)
      : std::ifstream(fname,std::ios::binary), doswap(doswap_) {}

    void open (const char *fname, bool doswap_)
      {
      doswap=doswap_;
      std::ifstream::open(fname,std::ios::binary);
      }

    template<typename T> bifstream &operator>> (T &data)
      {
      read (reinterpret_cast<char *> (&data), sizeof(T));
      if (doswap) byteswap (data);
      return *this;
      }
    template<typename T> bifstream &get (T *data, size_t num)
      {
      read (reinterpret_cast<char *> (data), num*sizeof(T));
      if (doswap)
        for (size_t m=0; m<num; ++m)
          byteswap (data[m]);
      return *this;
      }

    void rewind()
      { seekg(0,std::ios::beg); }
    void skip(std::streamoff nbytes)
      { seekg(nbytes,std::ios::cur); }

    bool getSwap() const
      { return doswap; }
    void setSwap(bool newswap)
      { doswap=newswap; }
    void flipSwap()
      { doswap=!doswap; }
  };

#endif

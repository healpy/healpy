/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file fitshandle.h
 *  Declaration of the FITS I/O helper class used by LevelS
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_FITSHANDLE_H
#define PLANCK_FITSHANDLE_H

#include <string>
#include <vector>
#include "fitsio.h"
#include "arr.h"
#include "datatypes.h"

/*! \defgroup fitsgroup FITS-related functionality */
/*! \{ */

template<typename T> struct FITSUTIL {};

template<> struct FITSUTIL<signed char>
  { enum { DTYPE=TBYTE }; };
template<> struct FITSUTIL<short>
  { enum { DTYPE=TSHORT }; };
template<> struct FITSUTIL<int>
  { enum { DTYPE=TINT }; };
template<> struct FITSUTIL<long>
  { enum { DTYPE=TLONG }; };
template<> struct FITSUTIL<long long>
  { enum { DTYPE=TLONGLONG }; };
template<> struct FITSUTIL<float>
  { enum { DTYPE=TFLOAT }; };
template<> struct FITSUTIL<double>
  { enum { DTYPE=TDOUBLE }; };

/*! Converts a FITS type code (i.e. the data type of a FITS column) to the
    corresponding Planck type code. */
inline int ftc2type (int ftc)
  {
  switch (ftc)
    {
    case TLOGICAL : return PLANCK_BOOL;
    case TBYTE    : return PLANCK_INT8;
    case TSHORT   : return PLANCK_INT16;
    case TINT32BIT: return PLANCK_INT32;
    case TLONGLONG: return PLANCK_INT64;
    case TFLOAT   : return PLANCK_FLOAT32;
    case TDOUBLE  : return PLANCK_FLOAT64;
    case TSTRING  : return PLANCK_STRING;
    default: throw Message_error ("ftc2type: unsupported component type");
    }
  }

/*! Converts a Planck type code to the corresponding FITS type code
    (i.e. the data type of a FITS column). */
inline int type2ftc (int type)
  {
  switch (type)
    {
    case PLANCK_BOOL   : return TLOGICAL;
    case PLANCK_INT8   : return TBYTE;
    case PLANCK_INT16  : return TSHORT;
    case PLANCK_INT32  : return TINT32BIT;
    case PLANCK_INT64  : return TLONGLONG;
    case PLANCK_FLOAT32: return TFLOAT;
    case PLANCK_FLOAT64: return TDOUBLE;
    case PLANCK_STRING : return TSTRING;
    default: throw Message_error ("type2ftc: unsupported component type");
    }
  }

/*! Class containing information about a single column in a FITS table. */
class fitscolumn
  {
  private:
    std::string name_, unit_;
    int64 repcount_;
    int type_;

  public:
    fitscolumn()
      : repcount_(0), type_(-1) {}
    /*! Creates a \a fitscolumn with name \a nm, unit \a un, a repetition
        count of \a rc, and a FITS type code of \a tp. */
    fitscolumn (const std::string &nm, const std::string &un, int64 rc,
      int tp)
      : name_(nm), unit_(un), repcount_(rc), type_(tp) {}

    /*! Returns the column name. */
    const std::string &name() const {return name_;}
    /*! Returns the column unit string. */
    const std::string &unit() const {return unit_;}
    /*! Returns the repetition count of the column. */
    int64 repcount() const {return repcount_;}
    /*! Returns the FITS type code of the column. */
    int type() const {return type_;}
  };

/*! Class for performing I/O from/to FITS files. */
class fitshandle
  {
  private:
    enum { INVALID = -4711 };

    mutable int status;
    fitsfile *fptr;
    int hdutype_, bitpix_;
    std::vector<int64> axes_;
    std::vector<fitscolumn> columns_;
    int64 nrows_;

    void check_errors() const;

    void clean_data();
    void clean_all();

    void assert_connected (const std::string &func) const
      {
      planck_assert (hdutype_!=INVALID,
        func + ": not connected to a HDU");
      }
    void assert_table_hdu (const std::string &func, unsigned int col) const
      {
      planck_assert ((hdutype_==ASCII_TBL) || (hdutype_==BINARY_TBL),
        func + ": HDU is not a table");
      planck_assert (col>0 && col<=columns_.size(),
        func + ": column number out of range");
      }
    void assert_image_hdu (const std::string &func) const
      {
      planck_assert ((hdutype_==IMAGE_HDU), func + ": HDU is not an image");
      }

    void init_image();
    void init_asciitab();
    void init_bintab();
    void init_data();

    void check_key_present(const std::string &name)const ;

    void read_col (int colnum, void *data, int64 ndata, int dtype,
                   int64 offset) const;
    void write_col (int colnum, const void *data, int64 ndata, int dtype,
                   int64 offset);

  public:
    /*! the list of modes in which a \a fitshandle can be opened. */
    typedef enum { CREATE, /*!< the file must not yet exist */
                   OPEN    /*!< the file must already exist */
                 } openmethod;

    /*! \name File-level access and manipulation. */
    /*! \{ */

    /*! Creates an unconnected \a fitshandle. */
    fitshandle ()
      : status(0), fptr(0), hdutype_(INVALID), bitpix_(INVALID), nrows_(0) {}
    /*! Creates a \a fitshandle connected to file \a fname.
        If \a rwmode == READONLY, no writing access is permitted; if it is
        READWRITE, reading and writing can be performed. */
    fitshandle (const std::string &fname, openmethod mode=OPEN,
      int rwmode=READONLY)
      : status(0), fptr(0), hdutype_(INVALID), bitpix_(INVALID), nrows_(0)
      {
      if (mode==OPEN)
        open (fname, rwmode);
      else
        create (fname);
      }
    /*! Creates a \a fitshandle connected to file \a fname and jumps directly
        to the HDU with the number \a hdunum.
        If \a rwmode == READONLY, no writing access is permitted; if it is
        READWRITE, reading and writing can be performed. */
    fitshandle (const std::string &fname, int hdunum, int rwmode=READONLY)
      : status(0), fptr(0), hdutype_(INVALID), bitpix_(INVALID), nrows_(0)
      {
      open (fname, rwmode);
      goto_hdu (hdunum);
      }

    /*! Performs all necessary cleanups. */
    ~fitshandle() { clean_all(); }

    /*! Connects to the file \a fname.
        If \a rwmode == READONLY, no writing access is permitted; if it is
        READWRITE, reading and writing can be performed. */
    void open (const std::string &fname, int rwmode=READONLY);
    /*! Creates the file \a fname and connects to it. */
    void create (const std::string &fname);
    /*! Closes the current file. */
    void close () { clean_all(); }
    /*! Deletes the file with name \a name. */
    static void delete_file (const std::string &name);
    /*! Jumps to the HDU with the absolute number \a hdu. */
    void goto_hdu (int hdu);
    /*! Returns the number of HDUs in the file. */
    int num_hdus () const;
    /*! Asserts that the PDMTYPE of the current HDU is \a pdmtype. */
    void assert_pdmtype (const std::string &pdmtype) const;
    /*! Inserts a binary table described by \a cols.
        The HDU has the name \a extname. */
    void insert_bintab (const std::vector<fitscolumn> &cols,
      const std::string &extname="xtension");
    /*! Inserts an ASCII table described by \a cols. The width of the
        columns is chosen automatically, in a way that avoids truncation.
        The HDU has the name \a extname. */
    void insert_asctab (const std::vector<fitscolumn> &cols,
      const std::string &extname="xtension");
    /*! Inserts a FITS image with the type given by \a btpx and dimensions
        given by \a Axes. */
    void insert_image (int btpx, const std::vector<int64> &Axes);
    /*! Inserts a 2D FITS image with the type given by \a btpx, whose
        contents are given in \a data. */
    template<typename T>
      void insert_image (int btpx, const arr2<T> &data);

    /*! Computes the checksum for the current HDU and writes it into the
        header. */
    void write_checksum();

    /*! \} */

    /*! \name Information about the current HDU */
    /*! \{ */

    /*! If the current HDU is an image, returns the BITPIX parameter of that
        image, else throws an exception. */
    int bitpix() const
      {
      assert_image_hdu ("fitshandle::bitpix()");
      return bitpix_;
      }
    /*! Returns the FITS type code for the current HDU. */
    int hdutype() const {return hdutype_;}
    /*! Returns the dimensions of the current image. */
    const std::vector<int64> &axes() const
      {
      assert_image_hdu ("fitshandle::axes()");
      return axes_;
      }
    /*! Returns the name of column \a #i. */
    const std::string &colname(int i) const
      {
      assert_table_hdu("fitshandle::colname()",i);
      return columns_[i-1].name();
      }
    /*! Returns the unit of column \a #i. */
    const std::string &colunit(int i) const
      {
      assert_table_hdu("fitshandle::colunit()",i);
      return columns_[i-1].unit();
      }
    /*! Returns repetition count of column \a #i. */
    int64 repcount(int i) const
      {
      assert_table_hdu("fitshandle::repcount()",i);
      return columns_[i-1].repcount();
      }
    /*! Returns the FITS type code for column \a #i. */
    int coltype(int i) const
      {
      assert_table_hdu("fitshandle::coltype()",i);
      return columns_[i-1].type();
      }
    /*! Returns the number of columns in the current table. */
    int ncols() const
      {
      assert_table_hdu("fitshandle::ncols()",1);
      return columns_.size();
      }
    /*! Returns the number of rows in the current table. */
    int64 nrows() const
      {
      assert_table_hdu("fitshandle::nrows()",1);
      return nrows_;
      }
    /*! Returns the total number of elements (nrows*repcount)
        in column \a #i. */
    int64 nelems(int i) const
      {
      assert_table_hdu("fitshandle::nelems()",i);
      if (columns_[i-1].type()==TSTRING) return nrows_;
      return nrows_*columns_[i-1].repcount();
      }

    /*! \} */

    /*! \name Keyword-handling methods */
    /*! \{ */

    /*! Copies all header keywords from the current HDU of \a orig to
        the current HDU of \a *this. */
    void copy_header (const fitshandle &orig);
    /*! Copies all header keywords from the current HDU of \a orig to
        the current HDU of \a *this, prepending a HISTORY keyword to
        every line. */
    void copy_historified_header (const fitshandle &orig);

    /*! Returns a list of all user-defined keys in the current HDU
        in \a keys. */
    void get_all_keys (std::vector<std::string> &keys) const;

    /*! Adds a new header line consisting of \a key, \a value and
        \a comment. */
    template<typename T> void add_key (const std::string &name, const T &value,
      const std::string &comment="");
    /*! Updates \a key with \a value and \a comment. */
    template<typename T> void update_key (const std::string &name,
      const T &value, const std::string &comment="");
    /*! Deletes \a key from the header. */
    void delete_key (const std::string &name);
    /*! Adds \a comment as a comment line. */
    void add_comment (const std::string &comment);
    /*! Reads the value belonging to \a key and returns it in \a value. */
    template<typename T> void get_key (const std::string &name, T &value) const;
    /*! Returms the value belonging to \a key. */
    template<typename T> T get_key (const std::string &name) const
      { T tmp; get_key(name, tmp); return tmp; }
    /*! Returns \a true if \a key is present, else \a false. */
    bool key_present (const std::string &name) const;
    /*! Returns the Planck type code for the key \a name. */
    int get_key_type(const std::string &name) const;

    /*! \} */

    /*! \name Methods for table data I/O */
    /*! \{ */

    void read_column_raw_void
      (int colnum, void *data, int type, int64 num, int64 offset=0) const;
    /*! Copies \a num elements from column \a colnum to the memory pointed
        to by \a data, starting at offset \a offset in the column. */
    template<typename T> void read_column_raw
      (int colnum, T *data, int64 num, int64 offset=0) const
      { read_column_raw_void (colnum, data, typehelper<T>::id, num, offset); }
    /*! Fills \a data with elements from column \a colnum,
        starting at offset \a offset in the column. */
    template<typename T> void read_column
      (int colnum, arr<T> &data, int64 offset=0) const
      { read_column_raw (colnum, &(data[0]), data.size(), offset); }
    /*! Reads the element \a #offset from column \a colnum into \a data. */
    template<typename T> void read_column
      (int colnum, T &data, int64 offset=0) const
      { read_column_raw (colnum, &data, 1, offset); }
    /* Reads the whole column \a colnum into \a data (which is resized
       accordingly). */
    template<typename T> void read_entire_column
      (int colnum, arr<T> &data) const
      { data.alloc(nelems(colnum)); read_column (colnum, data); }


    void write_column_raw_void
      (int colnum, const void *data, int type, int64 num, int64 offset=0);
    /*! Copies \a num elements from the memory pointed to by \a data to the
        column \a colnum, starting at offset \a offset in the column. */
    template<typename T> void write_column_raw
      (int colnum, const T *data, int64 num, int64 offset=0)
      { write_column_raw_void (colnum, data, typehelper<T>::id, num, offset); }
    /*! Copies all elements from \a data to the
        column \a colnum, starting at offset \a offset in the column. */
    template<typename T> void write_column
      (int colnum, const arr<T> &data, int64 offset=0)
      { write_column_raw (colnum, &(data[0]), data.size(), offset); }
    /*! Copies \a data to the column \a colnum, at the position \a offset. */
    template<typename T> void write_column
      (int colnum, const T &data, int64 offset=0)
      { write_column_raw (colnum, &data, 1, offset); }

    /*! \} */

    /*! \name Methods for image data I/O */
    /*! \{ */

    /*! Reads the current image into \a data, which is resized accordingly. */
    template<typename T> void read_image (arr2<T> &data) const;
    /*! Reads the current image into \a data, which is resized accordingly. */
    template<typename T> void read_image (arr3<T> &data) const;
    /*! Reads a partial image, whose dimensions are given by the dimensions
        of \a data, into data. The starting pixel indices are given by
        \a xl and \a yl. */
    template<typename T> void read_subimage
      (arr2<T> &data, int xl, int yl) const;
    /*! Fills \a data with values from the image, starting at the offset
        \a offset in the image. The image is treated as a one-dimensional
        array. */
    template<typename T> void read_subimage (arr<T> &data, int64 offset=0)
      const;
    /*! Writes \a data into the current image. \a data must have the same
        dimensions as specified in the HDU. */
    template<typename T> void write_image (const arr2<T> &data);
    /*! Copies \a data to the image, starting at the offset
        \a offset in the image. The image is treated as a one-dimensional
        array. */
    template<typename T> void write_subimage (const arr<T> &data,
      int64 offset=0);

    /*! \} */

    void add_healpix_keys (int datasize);
  };

/*! \} */

// announce the specialisations
template<> void fitshandle::add_key(const std::string &name,
  const bool &value, const std::string &comment);
template<> void fitshandle::add_key (const std::string &name,
  const std::string &value, const std::string &comment);
template<> void fitshandle::update_key(const std::string &name,
  const bool &value, const std::string &comment);
template<> void fitshandle::update_key (const std::string &name,
  const std::string &value, const std::string &comment);
template<> void fitshandle::get_key(const std::string &name,bool &value) const;
template<> void fitshandle::get_key (const std::string &name,
  std::string &value) const;

#endif

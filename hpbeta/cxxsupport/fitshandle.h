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

/*! \file fitshandle.h
 *  Declaration of the FITS I/O helper class used by LevelS
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_FITSHANDLE_H
#define PLANCK_FITSHANDLE_H

#include <string>
#include <vector>
#include "arr.h"
#include "datatypes.h"
#include "safe_cast.h"

/*! \defgroup fitsgroup FITS-related functionality */
/*! \{ */

/*! Class containing information about a single column in a FITS table. */
class fitscolumn
  {
  private:
    std::string name_, unit_;
    int64 repcount_;
    PDT type_;

  public:
    fitscolumn();
    /*! Creates a \a fitscolumn with name \a nm, unit \a un, a repetition
        count of \a rc, and a Planck type of \a tp. */
    fitscolumn (const std::string &nm, const std::string &un, int64 rc, PDT tp);
    ~fitscolumn();

    /*! Returns the column name. */
    const std::string &name() const {return name_;}
    /*! Returns the column unit string. */
    const std::string &unit() const {return unit_;}
    /*! Returns the repetition count of the column. */
    int64 repcount() const {return repcount_;}
    /*! Returns the Planck type of the column. */
    PDT type() const {return type_;}
  };

/*! Class for performing I/O from/to FITS files. */
class fitshandle
  {
  private:
    enum { INVALID = -4711 };

    mutable int status;
    void *fptr;
    int hdutype_, bitpix_;
    std::vector<int64> axes_;
    std::vector<fitscolumn> columns_;
    int64 nrows_;

    void check_errors() const;

    void clean_data();
    void clean_all();

    bool connected() const { return (hdutype_!=INVALID); }
    bool table_hdu (tsize col) const;
    bool image_hdu () const;

    void init_image();
    void init_asciitab();
    void init_bintab();
    void init_data();

    void read_col (int colnum, void *data, int64 ndata, PDT type,
                   int64 offset) const;
    void write_col (int colnum, const void *data, int64 ndata, PDT type,
                   int64 offset);

    void getKeyHelper(const std::string &name) const;

  public:
    /*! the list of modes in which a \a fitshandle can be opened. */
    enum openmethod { CREATE, /*!< the file must not yet exist */
                      OPEN    /*!< the file must already exist */
                    };

    /*! \name File-level access and manipulation. */
    /*! \{ */

    /*! Creates an unconnected \a fitshandle. */
    fitshandle ();
    /*! Performs all necessary cleanups. */
    ~fitshandle();

    /*! Connects to the file \a fname. */
    void open (const std::string &fname);
    /*! Creates the file \a fname and connects to it. */
    void create (const std::string &fname);
    /*! Closes the current file. */
    void close () { clean_all(); }
    /*! Deletes the file with name \a name. */
    static void delete_file (const std::string &name);
    /*! Returns the name of the connected file. */
    std::string fileName() const;
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
    /*! Inserts a FITS image with the type given by \a type and dimensions
        given by \a Axes. */
    void insert_image (PDT type, const std::vector<int64> &Axes);
    /*! Inserts a 2D FITS image with the type given by \a type, whose
        contents are given in \a data. */
    template<typename T>
      void insert_image (PDT type, const arr2<T> &data);

    /*! Computes the checksum for the current HDU and writes it into the
        header. */
    void write_checksum();

    /*! \} */

    /*! \name Information about the current HDU */
    /*! \{ */

    /*! Returns the dimensions of the current image. */
    const std::vector<int64> &axes() const;
    /*! Returns the name of column \a #i. */
    const std::string &colname(int i) const;
    /*! Returns the unit of column \a #i. */
    const std::string &colunit(int i) const;
    /*! Returns repetition count of column \a #i. */
    int64 repcount(int i) const;
    /*! Returns the Planck type for column \a #i. */
    PDT coltype(int i) const;
    /*! Returns the number of columns in the current table. */
    int ncols() const;
    /*! Returns the number of rows in the current table. */
    int64 nrows() const;
    /*! Returns the total number of elements (nrows*repcount)
        in column \a #i. */
    int64 nelems(int i) const;
    /*! Returns the number of elements that should be read/written in a single
        call for optimum performance. This depends on the current HDU. */
    int64 efficientChunkSize(int i) const;

    /*! \} */

    /*! \name Keyword-handling methods */
    /*! \{ */

    /*! Returns a list of all user-defined keys in the current HDU
        in \a keys. */
    void get_all_keys (std::vector<std::string> &keys) const;

    void set_key_void (const std::string &key, const void *value, PDT type,
      const std::string &comment="");
    /*! Updates \a key with \a value and \a comment. */
    template<typename T> void set_key (const std::string &name,
      const T &value, const std::string &comment="")
      { set_key_void (name, &value, planckType<T>(), comment); }
    /*! Deletes \a key from the header. */
    void delete_key (const std::string &name);
    /*! Adds \a comment as a comment line. */
    void add_comment (const std::string &comment);
    void get_key_void (const std::string &key, void *value, PDT type) const;
    /*! Reads the value belonging to \a key and returns it in \a value. */
    template<typename T> void get_key (const std::string &name, T &value) const
      { get_key_void (name,&value,planckType<T>()); }
    /*! Returms the value belonging to \a key. */
    template<typename T> T get_key (const std::string &name) const
      { T tmp; get_key(name, tmp); return tmp; }
    /*! Returns \a true if \a key is present, else \a false. */
    bool key_present (const std::string &name) const;

    /*! \} */

    /*! \name Methods for table data I/O */
    /*! \{ */

    void read_column_raw_void
      (int colnum, void *data, PDT type, int64 num, int64 offset=0) const;
    /*! Copies \a num elements from column \a colnum to the memory pointed
        to by \a data, starting at offset \a offset in the column. */
    template<typename T> void read_column_raw
      (int colnum, T *data, int64 num, int64 offset=0) const
      { read_column_raw_void (colnum, data, planckType<T>(), num, offset); }
    /*! Fills \a data with elements from column \a colnum,
        starting at offset \a offset in the column. */
    template<typename T> void read_column
      (int colnum, arr<T> &data, int64 offset=0) const
      { read_column_raw (colnum, &(data[0]), data.size(), offset); }
    /*! Reads the element \a #offset from column \a colnum into \a data. */
    template<typename T> void read_column
      (int colnum, T &data, int64 offset=0) const
      { read_column_raw (colnum, &data, 1, offset); }
    /*! Reads the whole column \a colnum into \a data (which is resized
       accordingly). */
    template<typename T> void read_entire_column
      (int colnum, arr<T> &data) const
      {
      data.alloc(safe_cast<tsize>(nelems(colnum)));
      read_column (colnum, data);
      }

    void write_column_raw_void
      (int colnum, const void *data, PDT type, int64 num, int64 offset=0);
    /*! Copies \a num elements from the memory pointed to by \a data to the
        column \a colnum, starting at offset \a offset in the column. */
    template<typename T> void write_column_raw
      (int colnum, const T *data, int64 num, int64 offset=0)
      { write_column_raw_void (colnum, data, planckType<T>(), num, offset); }
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
    void read_subimage_void (void *data, PDT type, tsize ndata, int64 offset=0)
      const;
    /*! Fills \a data with values from the image, starting at the offset
        \a offset in the image. The image is treated as a one-dimensional
        array. */
    template<typename T> void read_subimage (arr<T> &data, int64 offset=0)
      const
      { read_subimage_void (&data[0],planckType<T>(),data.size(),offset); }

    void write_image2D_void (const void *data, PDT type, tsize s1,
      tsize s2);
    /*! Writes \a data into the current image. \a data must have the same
        dimensions as specified in the HDU. */
    template<typename T> void write_image (const arr2<T> &data)
      {
      write_image2D_void (&data[0][0],planckType<T>(),data.size1(),
                          data.size2());
      }

    void write_subimage_void (const void *data, PDT type, tsize sz,
      int64 offset);
    /*! Copies \a data to the image, starting at the offset
        \a offset in the image. The image is treated as a one-dimensional
        array. */
    template<typename T> void write_subimage (const arr<T> &data,
      int64 offset=0)
      { write_subimage_void(&data[0],planckType<T>(),data.size(),offset); }

    /*! \} */
  };

/*! \} */

#endif

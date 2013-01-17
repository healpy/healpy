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

/*
 *  This file contains the implementation of the FITS I/O helper class
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <sstream>
#include <cstring>
#include <vector>
#include "fitsio.h"
#include "fitshandle.h"
#include "string_utils.h"

#define FPTR (static_cast<fitsfile *> (fptr))
#define OFPTR (static_cast<fitsfile *> (orig.fptr))

using namespace std;

namespace {

template<typename T> inline int fitsType();
template<> inline int fitsType<float> () { return TFLOAT;  }
template<> inline int fitsType<double>() { return TDOUBLE; }

int type2bitpix (PDT type)
  {
  switch (type)
    {
    case PLANCK_FLOAT32: return FLOAT_IMG;
    case PLANCK_FLOAT64: return DOUBLE_IMG;
    default: planck_fail ("unsupported component type");
    }
  }

/*! Converts a FITS type code (i.e. the data type of a FITS column) to the
    corresponding Planck type code. */
PDT ftc2type (int ftc)
  {
  switch (ftc)
    {
    case TLOGICAL : return PLANCK_BOOL;
    case TBYTE    : return PLANCK_INT8;
    case TSHORT   : return PLANCK_INT16;
    case TINT     :
    case TINT32BIT: return PLANCK_INT32;
    case TLONGLONG: return PLANCK_INT64;
    case TFLOAT   : return PLANCK_FLOAT32;
    case TDOUBLE  : return PLANCK_FLOAT64;
    case TSTRING  : return PLANCK_STRING;
    default: planck_fail ("unsupported component type");
    }
  }

/*! Converts a Planck type code to the corresponding FITS type code
    (i.e. the data type of a FITS column). */
int type2ftc (PDT type)
  {
  switch (type)
    {
    case PLANCK_BOOL   : return TLOGICAL;
    case PLANCK_INT8   :
    case PLANCK_UINT8  : return TBYTE;
    case PLANCK_INT16  : return TSHORT;
    case PLANCK_INT32  : return TINT;
    case PLANCK_INT64  : return TLONGLONG;
    case PLANCK_FLOAT32: return TFLOAT;
    case PLANCK_FLOAT64: return TDOUBLE;
    case PLANCK_STRING : return TSTRING;
    default: planck_fail ("unsupported component type");
    }
  }

const char *type2fitschar (PDT type)
  {
  switch (type)
    {
    case PLANCK_BOOL   : return "L";
    case PLANCK_FLOAT32: return "E";
    case PLANCK_FLOAT64: return "D";
    case PLANCK_INT8   :
    case PLANCK_UINT8  : return "B";
    case PLANCK_INT16  : return "I";
    case PLANCK_INT32  : return "J";
    case PLANCK_INT64  : return "K";
    case PLANCK_STRING : return "A";
    default:
      planck_fail(string("unknown data type ")+type2string(type));
    }
  }

const char *type2asciiform (PDT type)
  {
  switch (type)
    {
    case PLANCK_FLOAT32: return "E14.7";
    case PLANCK_FLOAT64: return "D23.15";
    case PLANCK_UINT8  : return "I3";
    case PLANCK_INT8   : return "I4";
    case PLANCK_INT16  : return "I6";
    case PLANCK_INT32  : return "I11";
    case PLANCK_INT64  : return "I22";
    default:
      planck_fail(string("unknown data type ")+type2string(type));
    }
  }

string fixkey (const string &key)
  {
  for (tsize m=0; m<key.size(); ++m)
    if (islower(key[m])) return string("HIERARCH "+key);

  return key;
  }
} // unnamed namespace

fitscolumn::fitscolumn()
  : repcount_(0), type_(PLANCK_INVALID) {}

fitscolumn::fitscolumn (const string &nm, const string &un, int64 rc, PDT tp)
  : name_(nm), unit_(un), repcount_(rc), type_(tp) {}

fitscolumn::~fitscolumn () {}

void fitshandle::check_errors() const
  {
  char msg[81];
  if (status==0)
    {
    while (fits_read_errmsg(msg))
      cerr << "STALE FITS ERROR MESSAGE: " << msg << endl;
    fits_clear_errmsg();
    return;
    }
  fits_get_errstatus (status, msg);
  cerr << msg << endl;
  while (fits_read_errmsg(msg)) cerr << msg << endl;
  fits_clear_errmsg();
  status=0;
  planck_fail("FITS error");
  }

void fitshandle::clean_data()
  {
  if (!fptr) return;
  axes_.clear();
  columns_.clear();
  hdutype_=INVALID;
  bitpix_=INVALID;
  nrows_=0;
  }

void fitshandle::clean_all()
  {
  if (!fptr) return;
  clean_data();
  fits_close_file (FPTR, &status);
  check_errors();
  fptr=0;
  }

bool fitshandle::table_hdu (tsize col) const
  {
  if ((hdutype_!=ASCII_TBL) && (hdutype_!=BINARY_TBL)) return false;
  if ((col<=0) || (col>columns_.size())) return false;
  return true;
  }
bool fitshandle::image_hdu () const
  { return hdutype_==IMAGE_HDU; }

void fitshandle::init_image()
  {
  int naxis;
  fits_get_img_type (FPTR, &bitpix_, &status);
  fits_get_img_dim (FPTR, &naxis, &status);
  check_errors();
  arr<LONGLONG> naxes(naxis);
  fits_get_img_sizell (FPTR, naxis, &naxes[0], &status);
  for (long m=0; m<naxis; ++m) axes_.push_back(naxes[naxis-m-1]);
  check_errors();
  }

void fitshandle::init_asciitab()
  {
  char ttype[81], tunit[81], tform[81];
  int ncol, typecode;
  fits_get_num_cols (FPTR, &ncol, &status);
  { LONGLONG tmp; fits_get_num_rowsll (FPTR, &tmp, &status); nrows_=tmp; }
  check_errors();
  for (int m=1; m<=ncol; ++m)
    {
    fits_get_acolparms (FPTR, m, ttype, 0, tunit, tform,
                        0, 0, 0, 0, &status);
    fits_ascii_tform (tform, &typecode, 0,0, &status);
    check_errors();
    columns_.push_back (fitscolumn (ttype,tunit,1,ftc2type(typecode)));
    }
  }

void fitshandle::init_bintab()
  {
  char ttype[81], tunit[81], tform[81];
  LONGLONG repc;
  int ncol, typecode;
  fits_get_num_cols (FPTR, &ncol, &status);
  { LONGLONG tmp; fits_get_num_rowsll (FPTR, &tmp, &status); nrows_=tmp; }
  check_errors();
  for (int m=1; m<=ncol; ++m)
    {
    fits_get_bcolparmsll (FPTR, m, ttype, tunit, tform, &repc,
                        0, 0, 0, 0, &status);
    fits_binary_tform (tform, &typecode, 0,0, &status);
    check_errors();
    columns_.push_back (fitscolumn (ttype,tunit,repc,ftc2type(typecode)));
    }
  }

void fitshandle::init_data()
  {
  clean_data();
  fits_get_hdu_type (FPTR, &hdutype_, &status);
  check_errors();
  switch(hdutype_)
    {
    case IMAGE_HDU:
      init_image(); break;
    case ASCII_TBL:
      init_asciitab(); break;
    case BINARY_TBL:
      init_bintab(); break;
    default:
      planck_fail("init_data(): unsupported HDU type"); break;
    }
  }

void fitshandle::read_col (int colnum, void *data, int64 ndata, PDT type,
  int64 offset) const
  {
  planck_assert(table_hdu(colnum),"incorrect FITS table access");
  int64 repc = columns_[colnum-1].repcount();
  planck_assert (ndata<=(repc*nrows_-offset),"read_column(): array too large");
  int64 frow = offset/repc+1;
  int64 felem = offset%repc+1;
  fits_read_col (FPTR, type2ftc(type), colnum, frow, felem, ndata, 0, data, 0,
    &status);
  check_errors();
  }
void fitshandle::write_col (int colnum, const void *data, int64 ndata,
  PDT type, int64 offset)
  {
  planck_assert(table_hdu(colnum),"incorrect FITS table access");
  int64 repc = columns_[colnum-1].repcount();
  int64 frow = offset/repc+1;
  int64 felem = offset%repc+1;
  fits_write_col (FPTR, type2ftc(type), colnum, frow, felem, ndata,
    const_cast<void *>(data), &status);
  nrows_ = max(nrows_,offset+ndata);
  check_errors();
  }

void fitshandle::getKeyHelper(const string &name) const
  {
  if (status==KEY_NO_EXIST)
    {
    fits_clear_errmsg();
    status=0;
    planck_fail("fitshandle::get_key(): key '"+name+"' not found");
    }
  check_errors();
  }

fitshandle::fitshandle ()
  : status(0), fptr(0), hdutype_(INVALID), bitpix_(INVALID), nrows_(0) {}

fitshandle::~fitshandle()
  { clean_all(); }

void fitshandle::open (const string &fname)
  {
  clean_all();
  fitsfile *ptr;
  fits_open_file(&ptr, fname.c_str(), READONLY, &status);
  fptr=ptr;
  check_errors();
  init_data();
  }

void fitshandle::create (const string &fname)
  {
  clean_all();
  fitsfile *ptr;
  fits_create_file(&ptr, fname.c_str(), &status);
  fptr=ptr;
  fits_write_imghdr(FPTR, 8, 0, 0, &status); // insert empty PHDU
  fits_write_date(FPTR, &status);
  check_errors();
  init_data();
  }

// static
void fitshandle::delete_file (const string &name)
  {
  fitsfile *ptr;
  int stat = 0;
  fits_open_file(&ptr, name.c_str(), READWRITE, &stat);
  fits_delete_file(ptr, &stat);
  if (stat==0) return;

  char msg[81];
  fits_get_errstatus (stat, msg);
  cerr << msg << endl;
  while (fits_read_errmsg(msg)) cerr << msg << endl;
  planck_fail("FITS error");
  }

string fitshandle::fileName() const
  {
  planck_assert(connected(),"handle not connected to a file");
  char *fname = new char[2048];
  fits_file_name(FPTR, fname, &status);
  check_errors();
  string result(fname);
  delete[] fname;
  return result;
  }

void fitshandle::goto_hdu (int hdu)
  {
  int curhdu;
  fits_get_hdu_num(FPTR,&curhdu);
  if (curhdu!=hdu)
    {
    fits_movabs_hdu(FPTR, hdu, &hdutype_, &status);
    check_errors();
    init_data();
    }
  }

int fitshandle::num_hdus () const
  {
  int result;
  fits_get_num_hdus (FPTR, &result, &status);
  check_errors();
  return result;
  }

void fitshandle::insert_bintab (const vector<fitscolumn> &cols,
  const string &extname)
  {
  clean_data();
  int ncol=cols.size();
  arr2b<char> ttype(ncol,81), tform(ncol,81), tunit(ncol,81);

  for (long m=0; m<ncol; ++m)
    {
    strcpy (ttype[m], cols[m].name().c_str());
    strcpy (tunit[m], cols[m].unit().c_str());
    ostringstream x;
    x << cols[m].repcount() << type2fitschar(cols[m].type());
    strcpy (tform[m], x.str().c_str());
    }
  fits_insert_btbl (FPTR, nrows_, ncol, ttype.p0(), tform.p0(), tunit.p0(),
    const_cast<char *>(extname.c_str()), 0, &status);
  check_errors();
  init_data();
  }

void fitshandle::insert_asctab (const vector<fitscolumn> &cols,
  const string &extname)
  {
  clean_data();
  int ncol=cols.size();
  arr2b<char> ttype(ncol,81), tform(ncol,81), tunit(ncol,81);

  for (long m=0; m<ncol; ++m)
    {
    strcpy (ttype[m], cols[m].name().c_str());
    strcpy (tunit[m], cols[m].unit().c_str());
    ostringstream x;
    if (cols[m].type()!=TSTRING)
      {
      planck_assert (cols[m].repcount()==1,"bad repcount for ASCII table");
      x << type2asciiform(cols[m].type());
      }
    else
      {
      x << "A" << dataToString(cols[m].repcount());
      }
    strcpy (tform[m], x.str().c_str());
    }
  fits_insert_atbl (FPTR, 0, nrows_, ncol, ttype.p0(), 0, tform.p0(),
    tunit.p0(), const_cast<char *>(extname.c_str()), &status);
  check_errors();
  init_data();
  }

void fitshandle::insert_image (PDT type, const vector<int64> &Axes)
  {
  clean_data();
  arr<LONGLONG> tmpax(Axes.size());
  for (long m=0; m<long(Axes.size()); m++) tmpax[m]=Axes[Axes.size()-1-m];
  fits_insert_imgll (FPTR, type2bitpix(type), Axes.size(), &tmpax[0], &status);
  check_errors();
  init_data();
  }

template<typename T>
  void fitshandle::insert_image (PDT type, const arr2<T> &data)
  {
  clean_data();
  arr<LONGLONG> tmpax(2);
  tmpax[0] = data.size2(); tmpax[1] = data.size1();
  fits_insert_imgll (FPTR, type2bitpix(type), 2, &tmpax[0], &status);
  arr2<T> &tmparr = const_cast<arr2<T> &> (data);
  fits_write_img (FPTR, fitsType<T>(), 1, tmpax[0]*tmpax[1],
    &tmparr[0][0], &status);
  check_errors();
  init_data();
  }

template void fitshandle::insert_image (PDT type, const arr2<double> &data);
template void fitshandle::insert_image (PDT type, const arr2<float> &data);

void fitshandle::write_checksum()
  {
  planck_assert(connected(),"handle not connected to a file");
  fits_write_chksum (FPTR, &status);
  check_errors();
  }

const vector<int64> &fitshandle::axes() const
  {
  planck_assert(image_hdu(),"not connected to an image");
  return axes_;
  }
const string &fitshandle::colname(int i) const
  {
  planck_assert(table_hdu(i),"incorrect FITS table access");
  return columns_[i-1].name();
  }
const string &fitshandle::colunit(int i) const
  {
  planck_assert(table_hdu(i),"incorrect FITS table access");
  return columns_[i-1].unit();
  }
int64 fitshandle::repcount(int i) const
  {
  planck_assert(table_hdu(i),"incorrect FITS table access");
  return columns_[i-1].repcount();
  }
PDT fitshandle::coltype(int i) const
  {
  planck_assert(table_hdu(i),"incorrect FITS table access");
  return columns_[i-1].type();
  }
int fitshandle::ncols() const
  {
  planck_assert(table_hdu(1),"incorrect FITS table access");
  return columns_.size();
  }
int64 fitshandle::nrows() const
  {
  planck_assert(table_hdu(1),"incorrect FITS table access");
  return nrows_;
  }
int64 fitshandle::nelems(int i) const
  {
  planck_assert(table_hdu(i),"incorrect FITS table access");
  if (columns_[i-1].type()==PLANCK_STRING) return nrows_;
  return nrows_*columns_[i-1].repcount();
  }
int64 fitshandle::efficientChunkSize(int i) const
  {
  planck_assert(table_hdu(1),"incorrect FITS table access");
  long int res;
  fits_get_rowsize(FPTR, &res, &status);
  planck_assert(res>=1,"bad recommended FITS chunk size");
  check_errors();
  return res*columns_[i-1].repcount();
  }

void fitshandle::get_all_keys(vector<string> &keys) const
  {
  keys.clear();
  char card[81];
  const char *inclist[] = { "*" };
  planck_assert(connected(),"handle not connected to a file");
  fits_read_record (FPTR, 0, card, &status);
  check_errors();
  while (true)
    {
    fits_find_nextkey (FPTR, const_cast<char **>(inclist), 1,
      0, 0, card, &status);
    if (status!=0) break;
    if (fits_get_keyclass(card)==TYP_USER_KEY)
      {
      char keyname[80];
      int dummy;
      fits_get_keyname(card, keyname, &dummy, &status);
      check_errors();
      keys.push_back(trim(keyname));
      }
    check_errors();
    }
  if (status==KEY_NO_EXIST) { fits_clear_errmsg(); status=0; }
  check_errors();
  }

void fitshandle::set_key_void (const string &key, const void *value,
  PDT type, const string &comment)
  {
  planck_assert(connected(),"handle not connected to a file");
  string key2 = fixkey(key);
  switch (type)
    {
    case PLANCK_INT8:
    case PLANCK_UINT8:
    case PLANCK_INT16:
    case PLANCK_INT32:
    case PLANCK_INT64:
    case PLANCK_FLOAT32:
    case PLANCK_FLOAT64:
      fits_update_key (FPTR, type2ftc(type), const_cast<char *>(key2.c_str()),
        const_cast<void *>(value), const_cast<char *>(comment.c_str()),
        &status);
      break;
    case PLANCK_BOOL:
      {
      int val = *(static_cast<const bool *>(value));
      fits_update_key (FPTR, TLOGICAL, const_cast<char *>(key2.c_str()),
        &val, const_cast<char *>(comment.c_str()), &status);
      break;
      }
    case PLANCK_STRING:
      {
      string val = *(static_cast<const string *>(value));
      fits_update_key_longstr (FPTR, const_cast<char *>(key2.c_str()),
        const_cast<char *>(val.c_str()), const_cast<char *>(comment.c_str()),
        &status);
      break;
      }
    default:
      planck_fail ("unsupported data type in set_key_void()");
    }
  check_errors();
  }

void fitshandle::get_key_void (const string &name, void *value, PDT type) const
  {
  planck_assert(connected(),"handle not connected to a file");
  switch (type)
    {
    case PLANCK_INT8:
    case PLANCK_UINT8:
    case PLANCK_INT16:
    case PLANCK_INT32:
    case PLANCK_INT64:
    case PLANCK_FLOAT32:
    case PLANCK_FLOAT64:
      fits_read_key (FPTR, type2ftc(type), const_cast<char *>(name.c_str()),
        value, 0, &status);
      getKeyHelper(name);
      break;
    case PLANCK_BOOL:
      {
      int val;
      fits_read_key (FPTR, TLOGICAL, const_cast<char *>(name.c_str()), &val, 0,
        &status);
      getKeyHelper(name);
      *(static_cast<bool *>(value))=val;
      break;
      }
    case PLANCK_STRING:
      {
      char *tmp=0;
      fits_read_key_longstr (FPTR, const_cast<char *>(name.c_str()), &tmp, 0,
        &status);
      getKeyHelper(name);
      *(static_cast<string *>(value))=tmp;
      if (tmp) free(tmp);
      break;
      }
    default:
      planck_fail ("unsupported data type in get_key_void()");
    }
  check_errors();
  }

void fitshandle::delete_key (const string &name)
  {
  planck_assert(connected(),"handle not connected to a file");
  fits_delete_key (FPTR, const_cast<char *>(name.c_str()), &status);
  check_errors();
  }

void fitshandle::add_comment(const string &comment)
  {
  planck_assert(connected(),"handle not connected to a file");
  fits_write_comment(FPTR,const_cast<char *>(comment.c_str()),&status);
  check_errors();
  }

bool fitshandle::key_present(const string &name) const
  {
  char card[81];
  planck_assert(connected(),"handle not connected to a file");
  fits_read_card(FPTR, const_cast<char *>(name.c_str()), card, &status);
  if (status==KEY_NO_EXIST)
    { fits_clear_errmsg(); status=0; return false; }
  check_errors();
  return true;
  }

void fitshandle::assert_pdmtype (const string &pdmtype) const
  {
  string type;
  get_key("PDMTYPE",type);
  if (pdmtype==type) return;
  cerr << "PDMTYPE " << pdmtype << " expected, but found " << type << endl;
  }

void fitshandle::read_column_raw_void
  (int colnum, void *data, PDT type, int64 num, int64 offset) const
  {
  switch (type)
    {
    case PLANCK_INT8:
    case PLANCK_UINT8:
    case PLANCK_INT16:
    case PLANCK_INT32:
    case PLANCK_INT64:
    case PLANCK_FLOAT32:
    case PLANCK_FLOAT64:
    case PLANCK_BOOL:
      read_col (colnum, data, num, type, offset); break;
    case PLANCK_STRING:
      {
      string *data2 = static_cast<string *> (data);
      planck_assert(table_hdu(colnum),"incorrect FITS table access");
      planck_assert (num<=(nrows_-offset),
        "read_column(): array too large");
      arr2b<char> tdata(safe_cast<tsize>(num),
                        safe_cast<tsize>(columns_[colnum-1].repcount()+1));
      int dispwidth;
      fits_get_col_display_width(FPTR, colnum, &dispwidth, &status);
      planck_assert(dispwidth<=columns_[colnum-1].repcount(),"column too wide");
      fits_read_col (FPTR, TSTRING, colnum, offset+1, 1, num,
        0, tdata.p0(), 0, &status);
      check_errors();
      for (long m=0;m<num;++m) data2[m]=tdata[m];
      break;
      }
    default:
      planck_fail ("unsupported data type in read_column_raw_void()");
    }
  }

void fitshandle::write_column_raw_void
  (int colnum, const void *data, PDT type, int64 num, int64 offset)
  {
  switch (type)
    {
    case PLANCK_INT8:
    case PLANCK_UINT8:
    case PLANCK_INT16:
    case PLANCK_INT32:
    case PLANCK_INT64:
    case PLANCK_FLOAT32:
    case PLANCK_FLOAT64:
    case PLANCK_BOOL:
      write_col (colnum, data, num, type, offset); break;
    case PLANCK_STRING:
      {
      const string *data2 = static_cast<const string *> (data);
      planck_assert(table_hdu(colnum),"incorrect FITS table access");
      tsize stringlen = safe_cast<tsize>(columns_[colnum-1].repcount()+1);
      arr2b<char> tdata(safe_cast<tsize>(num), stringlen);
      for (long m=0;m<num;++m)
        {
        strncpy(tdata[m],data2[m].c_str(),stringlen-1);
        tdata[m][stringlen-1] = '\0';
        }
      fits_write_col (FPTR, TSTRING, colnum, offset+1, 1, num,
        tdata.p0(), &status);
      nrows_ = max(nrows_,offset+num);
      check_errors();
      break;
      }
    default:
      planck_fail ("unsupported data type in write_column_raw_void()");
    }
  }

void fitshandle::write_image2D_void (const void *data, PDT type, tsize s1,
  tsize s2)
  {
  planck_assert(image_hdu(),"not connected to an image");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  planck_assert (axes_[0]==int64(s1), "wrong size of dimension 1");
  planck_assert (axes_[1]==int64(s2), "wrong size of dimension 2");

  fits_write_img (FPTR, type2ftc(type), 1, axes_[0]*axes_[1],
    const_cast<void *>(data), &status);
  check_errors();
  }

void fitshandle::write_subimage_void (const void *data, PDT type, tsize sz,
  int64 offset)
  {
  planck_assert(image_hdu(),"not connected to an image");
  fits_write_img (FPTR, type2ftc(type), 1+offset, sz, const_cast<void *>(data),
    &status);
  check_errors();
  }

template<typename T> void fitshandle::read_image (arr2<T> &data) const
  {
  planck_assert(image_hdu(),"not connected to an image");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  data.alloc(safe_cast<tsize>(axes_[0]), safe_cast<tsize>(axes_[1]));
  fits_read_img (FPTR, fitsType<T>(), 1, axes_[0]*axes_[1], 0, &data[0][0], 0,
    &status);
  check_errors();
  }

template void fitshandle::read_image (arr2<float> &data) const;
template void fitshandle::read_image (arr2<double> &data) const;

template<typename T> void fitshandle::read_image (arr3<T> &data) const
  {
  planck_assert(image_hdu(),"not connected to an image");
  planck_assert (axes_.size()==3, "wrong number of dimensions");
  data.alloc(safe_cast<tsize>(axes_[0]), safe_cast<tsize>(axes_[1]),
    safe_cast<tsize>(axes_[2]));
  fits_read_img (FPTR, fitsType<T>(), 1, axes_[0]*axes_[1]*axes_[2],
    0, &data(0,0,0), 0, &status);
  check_errors();
  }

template void fitshandle::read_image (arr3<float> &data) const;
template void fitshandle::read_image (arr3<double> &data) const;

template<typename T> void fitshandle::read_subimage
  (arr2<T> &data, int xl, int yl) const
  {
  planck_assert(image_hdu(),"not connected to an image");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  for (tsize m=0; m<data.size1(); ++m)
    fits_read_img (FPTR, fitsType<T>(), (xl+m)*axes_[1]+yl+1,
      data.size2(), 0, &data[m][0], 0, &status);
  check_errors();
  }

template void fitshandle::read_subimage
  (arr2<float> &data, int xl, int yl) const;
template void fitshandle::read_subimage
  (arr2<double> &data, int xl, int yl) const;

void fitshandle::read_subimage_void (void *data, PDT type, tsize ndata,
  int64 offset) const
  {
  planck_assert(image_hdu(),"not connected to an image");
  fits_read_img (FPTR, type2ftc(type), 1+offset, ndata, 0, data, 0, &status);
  check_errors();
  }

namespace {

class cfitsio_checker
  {
  public:
    cfitsio_checker()
      {
      float fitsversion;
      planck_assert(fits_get_version(&fitsversion),
        "error calling fits_get_version()");
      int v_header  = nearest<int>(1000.*CFITSIO_VERSION),
          v_library = nearest<int>(1000.*fitsversion);
      if (v_header!=v_library)
        cerr << endl << "WARNING: version mismatch between CFITSIO header (v"
             << dataToString(v_header*0.001) << ") and linked library (v"
             << dataToString(v_library*0.001) << ")." << endl << endl;
      }
  };

cfitsio_checker Cfitsio_Checker;

} // unnamed namespace

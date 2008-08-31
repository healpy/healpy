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

/*
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <string>
#include "alm_fitsio.h"
#include "alm.h"
#include "fitshandle.h"
#include "xcomplex.h"

using namespace std;

void get_almsize(fitshandle &inp, int &lmax, int &mmax)
  {
  if (inp.key_present("MAX-LPOL") && inp.key_present("MAX-MPOL"))
    {
    inp.get_key ("MAX-LPOL",lmax);
    inp.get_key ("MAX-MPOL",mmax);
    return;
    }

  int n_alms = inp.nelems(1);
  arr<int> index;
  const int chunksize=1024*256;
  int offset=0;
  lmax=-1;
  mmax=-1;
  while (offset<n_alms)
    {
    int ppix=min(chunksize,n_alms-offset);
    index.alloc(ppix);
    inp.read_column(1,index,offset);

    for (int i=0; i<ppix; ++i)
      {
      int l = isqrt(index[i]-1);
      int m = index[i] - l*l - l - 1;
      if (l>lmax) lmax=l;
      if (m>mmax) mmax=m;
      }
    offset+=chunksize;
    }
  }

void get_almsize(const string &filename, int &lmax, int &mmax, int hdunum)
  {
  fitshandle inp;
  inp.open (filename);
  inp.goto_hdu(hdunum);
  get_almsize (inp, lmax, mmax);
  }

void get_almsize_pol(const string &filename, int &lmax, int &mmax)
  {
  int tlmax, tmmax;
  fitshandle inp;
  inp.open (filename);
  lmax=mmax=0;
  for (int hdu=2; hdu<=4; ++hdu)
    {
    inp.goto_hdu(hdu);
    get_almsize (inp,tlmax,tmmax);
    if (tlmax>lmax) lmax=tlmax;
    if (tmmax>mmax) mmax=tmmax;
    }
  }

template<typename T> void read_Alm_from_fits
  (fitshandle &inp, Alm<xcomplex<T> >&alms, int lmax, int mmax)
  {
  int n_alms = inp.nelems(1);
  arr<int> index;
  arr<T> re, im;

  alms.Set(lmax, mmax);
  alms.SetToZero();
  int max_index = lmax*lmax + lmax + mmax + 1;
  const int chunksize=1024*256;
  int offset=0;
  while (offset<n_alms)
    {
    int ppix=min(chunksize,n_alms-offset);
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    inp.read_column(1,index,offset);
    inp.read_column(2,re,offset);
    inp.read_column(3,im,offset);

    for (int i=0; i<ppix; ++i)
      {
      if (index[i]>max_index) return;

      int l = isqrt(index[i]-1);
      int m = index[i] - l*l - l - 1;
      planck_assert(m>=0,"negative m encountered");
      planck_assert(l>=m, "wrong l,m combination");
      if ((l<=lmax) && (m<=mmax))
        alms(l,m).Set (re[i], im[i]);
      }
    offset+=chunksize;
    }
  }

template void read_Alm_from_fits (fitshandle &inp,
  Alm<xcomplex<double> > &alms, int lmax, int mmax);
template void read_Alm_from_fits (fitshandle &inp,
  Alm<xcomplex<float> > &alms, int lmax, int mmax);


template<typename T> void read_Alm_from_fits
  (const string &filename, Alm<xcomplex<T> >&alms, int lmax, int mmax,
  int hdunum)
  {
  fitshandle inp;
  inp.open (filename);
  inp.goto_hdu(hdunum);
  read_Alm_from_fits(inp,alms,lmax,mmax);
  }

template void read_Alm_from_fits (const string &filename,
  Alm<xcomplex<double> > &alms, int lmax, int mmax, int hdunum);
template void read_Alm_from_fits (const string &filename,
  Alm<xcomplex<float> > &alms, int lmax, int mmax, int hdunum);


template<typename T> void write_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<T> > &alms, int lmax, int mmax,
  int datatype)
  {
  vector<fitscolumn> cols;
  cols.push_back (fitscolumn("index","l*l+l+m+1",1,TINT32BIT));
  cols.push_back (fitscolumn("real","unknown",1,datatype));
  cols.push_back (fitscolumn("imag","unknown",1,datatype));
  out.insert_bintab(cols);
  arr<int> index;
  arr<double> re, im;

  int lm=alms.Lmax(), mm=alms.Mmax();
  int n_alms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);

  int l=0, m=0;
  const int chunksize=1024*256;
  int offset=0;
  while (offset<n_alms)
    {
    int ppix=min(chunksize,n_alms-offset);
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    for (int i=0; i<ppix; ++i)
      {
      index[i] = l*l + l + m + 1;
      if ((l<=lm) && (m<=mm))
        { re[i] = alms(l,m).re; im[i] = alms(l,m).im; }
      else
        { re[i] = 0; im[i] = 0; }
      ++m;
      if ((m>l) || (m>mmax)) { ++l; m=0; }
      }
    out.write_column(1,index,offset);
    out.write_column(2,re,offset);
    out.write_column(3,im,offset);

    offset+=chunksize;
    }
  out.add_key("MAX-LPOL",lmax,"highest l in the table");
  out.add_key("MAX-MPOL",mmax,"highest m in the table");
  }

template void write_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<double> > &alms, int lmax,
   int mmax, int datatype);
template void write_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<float> > &alms, int lmax,
   int mmax, int datatype);


template<typename T> void write_compressed_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<T> > &alms, int lmax, int mmax,
  int datatype)
  {
  vector<fitscolumn> cols;
  cols.push_back (fitscolumn("index","l*l+l+m+1",1,TINT32BIT));
  cols.push_back (fitscolumn("real","unknown",1,datatype));
  cols.push_back (fitscolumn("imag","unknown",1,datatype));
  out.insert_bintab(cols);
  arr<int> index;
  arr<double> re, im;

  int n_alms = 0;
  for (int m=0; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      if (alms(l,m).norm()>0) ++n_alms;

  int l=0, m=0;
  const int chunksize=1024*256;
  int real_lmax=0, real_mmax=0;
  int offset=0;
  while (offset<n_alms)
    {
    int ppix=min(chunksize,n_alms-offset);
    index.alloc(ppix);
    re.alloc(ppix); im.alloc(ppix);
    for (int i=0; i<ppix; ++i)
      {
      while (alms(l,m).norm()==0)
        {
        ++m;
        if ((m>l) || (m>mmax)) { ++l; m=0; }
        }
      index[i] = l*l + l + m + 1;
      re[i] = alms(l,m).re;
      im[i] = alms(l,m).im;
      if (l>real_lmax) real_lmax=l;
      if (m>real_mmax) real_mmax=m;
      ++m;
      if ((m>l) || (m>mmax)) { ++l; m=0; }
      }
    out.write_column(1,index,offset);
    out.write_column(2,re,offset);
    out.write_column(3,im,offset);

    offset+=chunksize;
    }
  out.add_key("MAX-LPOL",real_lmax,"highest l in the table");
  out.add_key("MAX-MPOL",real_mmax,"highest m in the table");
  }

template void write_compressed_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<double> > &alms, int lmax,
   int mmax, int datatype);
template void write_compressed_Alm_to_fits
  (fitshandle &out, const Alm<xcomplex<float> > &alms, int lmax,
   int mmax, int datatype);

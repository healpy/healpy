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
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2020 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ducc0/healpix/healpix_base.h"
#include "ducc0/math/geom_utils.h"
#include "ducc0/math/constants.h"
#include "ducc0/infra/mav.h"
#include "ducc0/math/space_filling.h"

namespace ducc0 {

namespace detail_healpix {

using namespace std;

namespace {

template<typename T> inline T spread(int v);
template<> inline int spread(int v)
  { return spread_bits_2D_32(v); }
template<> inline int64_t spread(int v)
  { return spread_bits_2D_64(v); }

template<typename T> inline T interleave(int x, int y);
template<> inline int interleave(int x, int y)
  { return coord2morton2D_32({uint32_t(x), uint32_t(y)}); }
template<> inline int64_t interleave(int x, int y)
  { return coord2morton2D_64({uint32_t(x), uint32_t(y)}); }

template<typename T> inline void deinterleave(T v, int &x, int &y);
template<> inline void deinterleave(int v, int &x, int &y)
  {
  auto res = morton2coord2D_32(v);
  x = res[0];
  y = res[1];
  }
template<> inline void deinterleave(int64_t v, int &x, int &y)
  {
  auto res = morton2coord2D_64(v);
  x = res[0];
  y = res[1];
  }

}

template<typename I> int T_Healpix_Base<I>::nside2order (I nside)
  {
  MR_assert (nside>I(0), "invalid value for Nside");
  return ((nside)&(nside-1)) ? -1 : ilog2(nside);
  }
template<typename I> I T_Healpix_Base<I>::npix2nside (I npix)
  {
  I res=isqrt(npix/I(12));
  MR_assert (npix==res*res*I(12), "invalid value for npix");
  return res;
  }

template<typename I> I T_Healpix_Base<I>::ring_above (double z) const
  {
  double az=abs(z);
  if (az<=twothird) // equatorial region
    return I(nside_*(2-1.5*z));
  I iring = I(nside_*sqrt(3*(1-az)));
  return (z>0) ? iring : 4*nside_-iring-1;
  }

namespace {

/* Short note on the "zone":
   zone = 0: pixel lies completely outside the queried shape
          1: pixel may overlap with the shape, pixel center is outside
          2: pixel center is inside the shape, but maybe not the complete pixel
          3: pixel lies completely inside the shape */

template<typename I, typename I2> inline void check_pixel (size_t o, size_t order_,
  size_t omax, size_t zone, rangeset<I2> &pixset, I pix, vector<pair<I,size_t> > &stk,
  bool inclusive, size_t &stacktop)
  {
  if (zone==0) return;

  if (o<order_)
    {
    if (zone>=3)
      {
      int sdist=2*(order_-o); // the "bit-shift distance" between map orders
      pixset.append(pix<<sdist,(pix+1)<<sdist); // output all subpixels
      }
    else // (1<=zone<=2)
      for (size_t i=0; i<4; ++i)
        stk.push_back(make_pair(4*pix+3-i,o+1)); // add children
    }
  else if (o>order_) // this implies that inclusive==true
    {
    if (zone>=2) // pixel center in shape
      {
      pixset.append(pix>>(2*(o-order_))); // output the parent pixel at order_
      stk.resize(stacktop); // unwind the stack
      }
    else // (zone==1): pixel center in safety range
      {
      if (o<omax) // check sublevels
        for (int i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
      else // at resolution limit
        {
        pixset.append(pix>>(2*(o-order_))); // output the parent pixel at order_
        stk.resize(stacktop); // unwind the stack
        }
      }
    }
  else // o==order_
    {
    if (zone>=2)
      pixset.append(pix);
    else if (inclusive) // and (zone>=1)
      {
      if (order_<omax) // check sublevels
        {
        stacktop=stk.size(); // remember current stack position
        for (size_t i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
        }
      else // at resolution limit
        pixset.append(pix); // output the pixel
      }
    }
  }

template<typename I> bool check_pixel_ring (const T_Healpix_Base<I> &b1,
  const T_Healpix_Base<I> &b2, I pix, I nr, I ipix1, int fct,
  double cz, double cphi, double cosrp2, I cpix)
  {
  if (pix>=nr) pix-=nr;
  if (pix<0) pix+=nr;
  pix+=ipix1;
  if (pix==cpix) return false; // disk center in pixel => overlap
  int px,py,pf;
  b1.pix2xyf(pix,px,py,pf);
  for (int i=0; i<fct-1; ++i) // go along the 4 edges
    {
    I ox=fct*px, oy=fct*py;
    double pz,pphi;
    b2.pix2zphi(b2.xyf2pix(ox+i,oy,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox+fct-1,oy+i,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox+fct-1-i,oy+fct-1,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    b2.pix2zphi(b2.xyf2pix(ox,oy+fct-1-i,pf),pz,pphi);
    if (cosdist_zphi(pz,pphi,cz,cphi)>cosrp2) // overlap
      return false;
    }
  return true;
  }

} // unnamed namespace

template<typename I> template<typename I2>
  void T_Healpix_Base<I>::query_disc_internal
  (pointing ptg, double radius, int fact, rangeset<I2> &pixset) const
  {
  bool inclusive = (fact!=0);
  pixset.clear();
  ptg.normalize();

  if (scheme_==RING)
    {
    I fct=1;
    if (inclusive)
      {
      MR_assert (((I(1)<<order_max)/nside_)>=fact,
        "invalid oversampling factor");
      fct = fact;
      }
    T_Healpix_Base b2;
    double rsmall, rbig;
    if (fct>1)
      {
      b2.SetNside(fct*nside_,RING);
      rsmall = radius+b2.max_pixrad();
      rbig = radius+max_pixrad();
      }
    else
      rsmall = rbig = inclusive ? radius+max_pixrad() : radius;

    if (rsmall>=pi)
      { pixset.append(0,npix_); return; }

    rbig = min(pi,rbig);

    double cosrsmall = cos(rsmall);
    double cosrbig = cos(rbig);

    double z0 = cos(ptg.theta);
    double xa = 1./sqrt((1-z0)*(1+z0));

    I cpix=zphi2pix(z0,ptg.phi);

    double rlat1 = ptg.theta - rsmall;
    double zmax = cos(rlat1);
    I irmin = ring_above (zmax)+1;

    if ((rlat1<=0) && (irmin>1)) // north pole in the disk
      {
      I sp,rp; bool dummy;
      get_ring_info_small(irmin-1,sp,rp,dummy);
      pixset.append(0,sp+rp);
      }

    if ((fct>1) && (rlat1>0)) irmin=max(I(1),irmin-1);

    double rlat2 = ptg.theta + rsmall;
    double zmin = cos(rlat2);
    I irmax = ring_above (zmin);

    if ((fct>1) && (rlat2<pi)) irmax=min(4*nside_-1,irmax+1);

    for (I iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      double x = (cosrbig-z*z0)*xa;
      double ysq = 1-z*z-x*x;
      double dphi=-1;
      if (ysq<=0) // no intersection, ring completely inside or outside
        dphi = (fct==1) ? 0: pi-1e-15;
      else
        dphi = atan2(sqrt(ysq),x);
      if (dphi>0)
        {
        I nr, ipix1;
        bool shifted;
        get_ring_info_small(iz,ipix1,nr,shifted);
        double shift = shifted ? 0.5 : 0.;

        I ipix2 = ipix1 + nr - 1; // highest pixel number in the ring

        I ip_lo = ifloor<I>(nr*inv_twopi*(ptg.phi-dphi) - shift)+1;
        I ip_hi = ifloor<I>(nr*inv_twopi*(ptg.phi+dphi) - shift);

        if (fct>1)
          {
          while ((ip_lo<=ip_hi) && check_pixel_ring
                (*this,b2,ip_lo,nr,ipix1,fct,z0,ptg.phi,cosrsmall,cpix))
            ++ip_lo;
          while ((ip_hi>ip_lo) && check_pixel_ring
                (*this,b2,ip_hi,nr,ipix1,fct,z0,ptg.phi,cosrsmall,cpix))
            --ip_hi;
          }

        if (ip_lo<=ip_hi)
          {
          if (ip_hi>=nr)
            { ip_lo-=nr; ip_hi-=nr; }
          if (ip_lo<0)
            {
            pixset.append(ipix1,ipix1+ip_hi+1);
            pixset.append(ipix1+ip_lo+nr,ipix2+1);
            }
          else
            pixset.append(ipix1+ip_lo,ipix1+ip_hi+1);
          }
        }
      }
    if ((rlat2>=pi) && (irmax+1<4*nside_)) // south pole in the disk
      {
      I sp,rp; bool dummy;
      get_ring_info_small(irmax+1,sp,rp,dummy);
      pixset.append(sp,npix_);
      }
    }
  else // scheme_==NEST
    {
    if (radius>=pi) // disk covers the whole sphere
      { pixset.append(0,npix_); return; }

    int oplus = 0;
    if (inclusive)
      {
      MR_assert ((I(1)<<(order_max-order_))>=fact,
        "invalid oversampling factor");
      MR_assert ((fact&(fact-1))==0,
        "oversampling factor must be a power of 2");
      oplus=ilog2(fact);
      }
    int omax=order_+oplus; // the order up to which we test

    vec3 vptg(ptg);
    vector<T_Healpix_Base<I> > base(omax+1);
    vector<double> crpdr(omax+1), crmdr(omax+1);
    double cosrad=cos(radius);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      crpdr[o] = (radius+dr>pi) ? -1. : cos(radius+dr);
      crmdr[o] = (radius-dr<0.) ?  1. : cos(radius-dr);
      }
    vector<pair<I,size_t> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(I(11-i),0));

    size_t stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      I pix=stk.back().first;
      int o=stk.back().second;
      stk.pop_back();

      double z,phi;
      base[o].pix2zphi(pix,z,phi);
      // cosine of angular distance between pixel center and disk center
      double cangdist=cosdist_zphi(vptg.z,ptg.phi,z,phi);

      if (cangdist>crpdr[o])
        {
        size_t zone = (cangdist<cosrad) ? 1 : ((cangdist<=crmdr[o]) ? 2 : 3);

        check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
          stacktop);
        }
      }
    }
  }

template<typename I> void T_Healpix_Base<I>::query_disc
  (pointing ptg, double radius, rangeset<I> &pixset) const
  {
  query_disc_internal (ptg, radius, 0, pixset);
  }

template<typename I> void T_Healpix_Base<I>::query_disc_inclusive
  (pointing ptg, double radius, rangeset<I> &pixset, int fact) const
  {
  MR_assert(fact>0,"fact must be a positive integer");
  if ((sizeof(I)<8) && (((I(1)<<order_max)/nside_)<fact))
    {
    T_Healpix_Base<int64_t> base2(nside_,scheme_,SET_NSIDE);
    base2.query_disc_internal(ptg,radius,fact,pixset);
    return;
    }
  query_disc_internal (ptg, radius, fact, pixset);
  }

template<typename I> template<typename I2>
  void T_Healpix_Base<I>::query_multidisc (const vector<vec3> &norm,
  const vector<double> &rad, int fact, rangeset<I2> &pixset) const
  {
  bool inclusive = (fact!=0);
  size_t nv=norm.size();
  MR_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    I fct=1;
    if (inclusive)
      {
      MR_assert (((I(1)<<order_max)/nside_)>=fact,
        "invalid oversampling factor");
      fct = fact;
      }
    T_Healpix_Base b2;
    double rpsmall, rpbig;
    if (fct>1)
      {
      b2.SetNside(fct*nside_,RING);
      rpsmall = b2.max_pixrad();
      rpbig = max_pixrad();
      }
    else
      rpsmall = rpbig = inclusive ? max_pixrad() : 0;

    I irmin=1, irmax=4*nside_-1;
    vector<double> z0,xa,cosrsmall,cosrbig;
    vector<pointing> ptg;
    vector<I> cpix;
    for (size_t i=0; i<nv; ++i)
      {
      double rsmall=rad[i]+rpsmall;
      if (rsmall<pi)
        {
        double rbig=min(pi,rad[i]+rpbig);
        pointing pnt=pointing(norm[i]);
        cosrsmall.push_back(cos(rsmall));
        cosrbig.push_back(cos(rbig));
        double cth=cos(pnt.theta);
        z0.push_back(cth);
        if (fct>1) cpix.push_back(zphi2pix(cth,pnt.phi));
        xa.push_back(1./sqrt((1-cth)*(1+cth)));
        ptg.push_back(pnt);

        double rlat1 = pnt.theta - rsmall;
        double zmax = cos(rlat1);
        I irmin_t = (rlat1<=0) ? 1 : ring_above (zmax)+1;

        if ((fct>1) && (rlat1>0)) irmin_t=max(I(1),irmin_t-1);

        double rlat2 = pnt.theta + rsmall;
        double zmin = cos(rlat2);
        I irmax_t = (rlat2>=pi) ? 4*nside_-1 : ring_above (zmin);

        if ((fct>1) && (rlat2<pi)) irmax_t=min(4*nside_-1,irmax_t+1);

        if (irmax_t < irmax) irmax=irmax_t;
        if (irmin_t > irmin) irmin=irmin_t;
        }
      }

    for (I iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      I ipix1,nr;
      bool shifted;
      get_ring_info_small(iz,ipix1,nr,shifted);
      double shift = shifted ? 0.5 : 0.;
      rangeset<I2> tr;
      tr.append(ipix1,ipix1+nr);
      for (size_t j=0; j<z0.size(); ++j)
        {
        double x = (cosrbig[j]-z*z0[j])*xa[j];
        double ysq = 1.-z*z-x*x;
        double dphi = (ysq<=0) ? pi-1e-15 : atan2(sqrt(ysq),x);
        I ip_lo = ifloor<I>(nr*inv_twopi*(ptg[j].phi-dphi) - shift)+1;
        I ip_hi = ifloor<I>(nr*inv_twopi*(ptg[j].phi+dphi) - shift);
        if (fct>1)
          {
          while ((ip_lo<=ip_hi) && check_pixel_ring
            (*this,b2,ip_lo,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
            ++ip_lo;
          while ((ip_hi>ip_lo) && check_pixel_ring
            (*this,b2,ip_hi,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
            --ip_hi;
          }
        if (ip_hi>=nr)
          { ip_lo-=nr; ip_hi-=nr;}
        if (ip_lo<0)
          tr.remove(ipix1+ip_hi+1,ipix1+ip_lo+nr);
        else
          tr.intersect(ipix1+ip_lo,ipix1+ip_hi+1);
        }
      pixset.append(tr);
      }
    }
  else // scheme_ == NEST
    {
    size_t oplus = 0;
    if (inclusive)
      {
      MR_assert ((I(1)<<(order_max-order_))>=fact,
        "invalid oversampling factor");
      MR_assert ((fact&(fact-1))==0,
        "oversampling factor must be a power of 2");
      oplus=ilog2(fact);
      }
    size_t omax=size_t(order_)+oplus; // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    vector<T_Healpix_Base<I> > base(omax+1);
    mav<double,3> crlimit({size_t(omax)+1,nv,3});
    for (size_t o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (size_t i=0; i<nv; ++i)
        {
        crlimit.v(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit.v(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit.v(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<I,size_t> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(I(11-i),0));

    size_t stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      I pix=stk.back().first;
      size_t o=stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      size_t zone=3;
      for (size_t i=0; i<nv; ++i)
        {
        double crad=dotprod(pv,norm[i]);
        for (size_t iz=0; iz<zone; ++iz)
          if (crad<crlimit(o,i,iz))
            if ((zone=iz)==0) goto bailout;
        }

      check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
        stacktop);
      bailout:;
      }
    }
  }

template<typename I> void T_Healpix_Base<I>::query_multidisc_general
  (const vector<vec3> &norm, const vector<double> &rad, bool inclusive,
  const vector<int> &cmds, rangeset<I> &pixset) const
  {
  size_t nv=norm.size();
  MR_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    MR_fail ("not yet implemented");
    }
  else // scheme_ == NEST
    {
    size_t oplus=inclusive ? 2 : 0;
    size_t omax=min<size_t>(order_max,size_t(order_)+oplus); // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    vector<T_Healpix_Base<I> > base(omax+1);
    mav<double,3> crlimit({size_t(omax+1),nv,3});
    for (size_t o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (size_t i=0; i<nv; ++i)
        {
        crlimit.v(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit.v(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit.v(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<I,size_t> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(I(11-i),0));

    size_t stacktop=0; // a place to save a stack position
    vector<size_t> zone(nv);

    vector<size_t> zstk; zstk.reserve(cmds.size());

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      I pix=stk.back().first;
      size_t o=stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      for (size_t i=0; i<nv; ++i)
        {
        zone[i]=3;
        double crad=dotprod(pv,norm[i]);
        for (size_t iz=0; iz<zone[i]; ++iz)
          if (crad<crlimit(o,i,iz))
            zone[i]=iz;
        }

      for (size_t i=0; i<cmds.size(); ++i)
        {
        size_t tmp;
        switch (cmds[i])
          {
          case -1: // union
            tmp=zstk.back(); zstk.pop_back();
            zstk.back() = max(zstk.back(),tmp);
            break;
          case -2: // intersection
            tmp=zstk.back(); zstk.pop_back();
            zstk.back() = min(zstk.back(),tmp);
            break;
          default: // add value
            zstk.push_back(zone[cmds[i]]);
          }
        }
      MR_assert(zstk.size()==1,"inconsistent commands");
      size_t zn=zstk[0]; zstk.pop_back();

      check_pixel (o, order_, omax, zn, pixset, pix, stk, inclusive,
        stacktop);
      }
    }
  }

template<typename I> inline void T_Healpix_Base<I>::nest2xyf (I pix, int &ix,
  int &iy, int &face_num) const
  {
  face_num = pix>>(2*order_);
  pix &= (npface_-1);
  deinterleave<I>(pix, ix, iy);
  }

template<typename I> inline I T_Healpix_Base<I>::xyf2nest (int ix, int iy,
  int face_num) const
  { return (I(face_num)<<(2*order_)) + interleave<I>(ix, iy); }

namespace {

// low-level hack to accelerate divisions with a result of [0;3]
template<typename I> inline I special_div(I a, I b)
  {
  I t=(a>=(b<<1));
  a-=t*(b<<1);
  return (t<<1)+(a>=b);
  }

} // unnamed namespace

template<typename I> void T_Healpix_Base<I>::ring2xyf (I pix, int &ix, int &iy,
  int &face_num) const
  {
  I iring, iphi, kshift, nr;
  I nl2 = 2*nside_;

  if (pix<ncap_) // North Polar cap
    {
    iring = (1+isqrt(1+2*pix))>>1; //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    face_num=special_div(iphi-1,nr);
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    I ip = pix - ncap_;
    I tmp = (order_>=0) ? ip>>(order_+2) : ip/(4*nside_);
    iring = tmp+nside_;
    iphi = ip-tmp*4*nside_ + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
    I ire = tmp+1,
      irm = nl2+1-tmp;
    I ifm = iphi - (ire>>1) + nside_ -1,
      ifp = iphi - (irm>>1) + nside_ -1;
    if (order_>=0)
      { ifm >>= order_; ifp >>= order_; }
    else
      { ifm /= nside_; ifp /= nside_; }
    face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else // South Polar cap
    {
    I ip = npix_ - pix;
    iring = (1+isqrt(2*ip-1))>>1; //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    face_num=special_div(iphi-1,nr)+8;
    }

  I irt = iring - ((2+(face_num>>2))*nside_) + 1;
  I ipt = 2*iphi- jpll[face_num]*nr - kshift -1;
//  I ipt = 2*iphi- (((face_num&3)<<1)+1-((face_num>>2)&1))*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  ix =  (ipt-irt) >>1;
  iy = (-ipt-irt) >>1;
  }

template<typename I> I T_Healpix_Base<I>::xyf2ring (int ix, int iy,
  int face_num) const
  {
  I nl4 = 4*nside_;
  I jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  I nr, kshift, n_before;

  bool shifted;
  get_ring_info_small(jr,n_before,nr,shifted);
  nr>>=2;
  kshift=1-shifted;
  I jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  MR_assert(jp<=4*nr,"must not happen");
  if (jp<1) jp+=nl4; // assumption: if this triggers, then nl4==4*nr

  return n_before + jp - 1;
  }

template<typename I> T_Healpix_Base<I>::T_Healpix_Base ()
  : order_(-1), nside_(0), npface_(0), ncap_(0), npix_(0),
    fact1_(0), fact2_(0), scheme_(RING) {}

template<typename I> void T_Healpix_Base<I>::Set (int order,
  Ordering_Scheme scheme)
  {
  MR_assert ((order>=0)&&(order<=order_max), "bad order");
  order_  = order;
  nside_  = I(1)<<order;
  npface_ = nside_<<order_;
  ncap_   = (npface_-nside_)<<1;
  npix_   = 12*npface_;
  fact2_  = 4./npix_;
  fact1_  = (nside_<<1)*fact2_;
  scheme_ = scheme;
  }
template<typename I> void T_Healpix_Base<I>::SetNside (I nside,
  Ordering_Scheme scheme)
  {
  order_  = nside2order(nside);
  MR_assert ((scheme!=NEST) || (order_>=0),
    "SetNside: nside must be power of 2 for nested maps");
  nside_  = nside;
  npface_ = nside_*nside_;
  ncap_   = (npface_-nside_)<<1;
  npix_   = 12*npface_;
  fact2_  = 4./npix_;
  fact1_  = (nside_<<1)*fact2_;
  scheme_ = scheme;
  }

template<typename I> double T_Healpix_Base<I>::ring2z (I ring) const
  {
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
  }

template<typename I> I T_Healpix_Base<I>::pix2ring (I pix) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      return (1+I(isqrt(1+2*pix)))>>1; // counted from North pole
    else if (pix<(npix_-ncap_)) // Equatorial region
      return (pix-ncap_)/(4*nside_) + nside_; // counted from North pole
    else // South Polar cap
      return 4*nside_-((1+I(isqrt(2*(npix_-pix)-1)))>>1);
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);
    return (I(jrll[face_num])<<order_) - ix - iy - 1;
    }
  }

template<typename I> I T_Healpix_Base<I>::nest2ring (I pix) const
  {
  MR_assert(order_>=0, "hierarchical map required");
  int ix, iy, face_num;
  nest2xyf (pix, ix, iy, face_num);
  return xyf2ring (ix, iy, face_num);
  }

template<typename I> I T_Healpix_Base<I>::ring2nest (I pix) const
  {
  MR_assert(order_>=0, "hierarchical map required");
  int ix, iy, face_num;
  ring2xyf (pix, ix, iy, face_num);
  return xyf2nest (ix, iy, face_num);
  }

template<typename I> inline I T_Healpix_Base<I>::nest_peano_helper
  (I pix, int dir) const
  {
  int face = (pix>>(2*order_));
  I result = 0;
  int state = ((peano_face2path[dir][face]<<4))|(dir<<7);
  int shift=2*order_-4;
  for (; shift>=0; shift-=4)
    {
    state=peano_arr2[(state&0xF0) | ((pix>>shift)&0xF)];
    result = (result<<4) | (state&0xF);
    }
  if (shift==-2)
    {
    state=peano_arr[((state>>2)&0xFC) | (pix&0x3)];
    result = (result<<2) | (state&0x3);
    }

  return result + (I(peano_face2face[dir][face])<<(2*order_));
  }

template<typename I> I T_Healpix_Base<I>::nest2peano (I pix) const
  { return nest_peano_helper(pix,0); }

template<typename I> I T_Healpix_Base<I>::peano2nest (I pix) const
  { return nest_peano_helper(pix,1); }

template<typename I> I T_Healpix_Base<I>::loc2pix (double z, double phi,
  double sth, bool have_sth) const
  {
  double za = abs(z);
  double tt = fmodulo(phi*inv_halfpi,4.0); // in [0,4)

  if (scheme_==RING)
    {
    if (za<=twothird) // Equatorial region
      {
      I nl4 = 4*nside_;
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*z*0.75;
      I jp = I(temp1-temp2); // index of  ascending edge line
      I jm = I(temp1+temp2); // index of descending edge line

      // ring number counted from z=2/3
      I ir = nside_ + 1 + jp - jm; // in {1,2n+1}
      I kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

      I t1 = jp+jm-nside_+kshift+1+nl4+nl4;
      I ip = (order_>0) ?
        (t1>>1)&(nl4-1) : ((t1>>1)%nl4); // in {0,4n-1}

      return ncap_ + (ir-1)*nl4 + ip;
      }
    else  // North & South polar caps
      {
      double tp = tt-I(tt);
      double tmp = ((za<0.99)||(!have_sth)) ?
                   nside_*sqrt(3*(1-za)) :
                   nside_*sth/sqrt((1.+za)/3.);

      I jp = I(tp*tmp); // increasing edge line index
      I jm = I((1.0-tp)*tmp); // decreasing edge line index

      I ir = jp+jm+1; // ring number counted from the closest pole
      I ip = I(tt*ir); // in {0,4*ir-1}
      MR_assert((ip>=0)&&(ip<4*ir),"must not happen");
      //ip %= 4*ir;

      return (z>0)  ?  2*ir*(ir-1) + ip  :  npix_ - 2*ir*(ir+1) + ip;
      }
    }
  else // scheme_ == NEST
    {
    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*(z*0.75);
      I jp = I(temp1-temp2); // index of  ascending edge line
      I jm = I(temp1+temp2); // index of descending edge line
      I ifp = jp >> order_;  // in {0,4}
      I ifm = jm >> order_;
      int face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));

      int ix = jm & (nside_-1),
          iy = nside_ - (jp & (nside_-1)) - 1;
      return xyf2nest(ix,iy,face_num);
      }
    else // polar region, za > 2/3
      {
      int ntt = min(3,int(tt));
      double tp = tt-ntt;
      double tmp = ((za<0.99)||(!have_sth)) ?
                   nside_*sqrt(3*(1-za)) :
                   nside_*sth/sqrt((1.+za)/3.);

      I jp = I(tp*tmp); // increasing edge line index
      I jm = I((1.0-tp)*tmp); // decreasing edge line index
      jp=min(jp,nside_-1); // for points too close to the boundary
      jm=min(jm,nside_-1);
      return (z>=0) ?
        xyf2nest(nside_-jm -1,nside_-jp-1,ntt) : xyf2nest(jp,jm,ntt+8);
      }
    }
  }

template<typename I> void T_Healpix_Base<I>::pix2loc (I pix, double &z,
  double &phi, double &sth, bool &have_sth) const
  {
  have_sth=false;
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      {
      I iring = (1+I(isqrt(1+2*pix)))>>1; // counted from North pole
      I iphi  = (pix+1) - 2*iring*(iring-1);

      double tmp=(iring*iring)*fact2_;
      z = 1.0 - tmp;
      if (z>0.99) { sth=sqrt(tmp*(2.-tmp)); have_sth=true; }
      phi = (iphi-0.5) * halfpi/iring;
      }
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      I nl4 = 4*nside_;
      I ip  = pix - ncap_;
      I tmp = (order_>=0) ? ip>>(order_+2) : ip/nl4;
      I iring = tmp + nside_,
        iphi = ip-nl4*tmp+1;
      // 1 if iring+nside is odd, 1/2 otherwise
      double fodd = ((iring+nside_)&1) ? 1 : 0.5;

      z = (2*nside_-iring)*fact1_;
      phi = (iphi-fodd) * pi*0.75*fact1_;
      }
    else // South Polar cap
      {
      I ip = npix_ - pix;
      I iring = (1+I(isqrt(2*ip-1)))>>1; // counted from South pole
      I iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

      double tmp=(iring*iring)*fact2_;
      z = tmp - 1.0;
      if (z<-0.99) { sth=sqrt(tmp*(2.-tmp)); have_sth=true; }
      phi = (iphi-0.5) * halfpi/iring;
      }
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);

    I jr = (I(jrll[face_num])<<order_) - ix - iy - 1;

    I nr;
    if (jr<nside_)
      {
      nr = jr;
      double tmp=(nr*nr)*fact2_;
      z = 1 - tmp;
      if (z>0.99) { sth=sqrt(tmp*(2.-tmp)); have_sth=true; }
      }
    else if (jr > 3*nside_)
      {
      nr = nside_*4-jr;
      double tmp=(nr*nr)*fact2_;
      z = tmp - 1;
      if (z<-0.99) { sth=sqrt(tmp*(2.-tmp)); have_sth=true; }
      }
    else
      {
      nr = nside_;
      z = (2*nside_-jr)*fact1_;
      }

    I tmp=I(jpll[face_num])*nr+ix-iy;
    MR_assert(tmp<8*nr,"must not happen");
    if (tmp<0) tmp+=8*nr;
    phi = (nr==nside_) ? 0.75*halfpi*tmp*fact1_ :
                         (0.5*halfpi*tmp)/nr;
    }
  }

template<typename I> template<typename I2>
  void T_Healpix_Base<I>::query_polygon_internal
  (const vector<pointing> &vertex, int fact, rangeset<I2> &pixset) const
  {
  bool inclusive = (fact!=0);
  size_t nv=vertex.size();
  size_t ncirc = inclusive ? nv+1 : nv;
  MR_assert(nv>=3,"not enough vertices in polygon");
  vector<vec3> vv(nv);
  for (size_t i=0; i<nv; ++i)
    vv[i]=vertex[i].to_vec3();
  vector<vec3> normal(ncirc);
  int flip=0;
  for (size_t i=0; i<nv; ++i)
    {
    normal[i]=crossprod(vv[i],vv[(i+1)%nv]).Norm();
    double hnd=dotprod(normal[i],vv[(i+2)%nv]);
    MR_assert(abs(hnd)>1e-10,"degenerate corner");
    if (i==0)
      flip = (hnd<0.) ? -1 : 1;
    else
      MR_assert(flip*hnd>0,"polygon is not convex");
    normal[i]*=flip;
    }
  vector<double> rad(ncirc,halfpi);
  if (inclusive)
    {
    double cosrad;
    find_enclosing_circle (vv, normal[nv], cosrad);
    rad[nv]=acos(cosrad);
    }
  query_multidisc(normal,rad,fact,pixset);
  }

template<typename I> void T_Healpix_Base<I>::query_polygon
  (const vector<pointing> &vertex, rangeset<I> &pixset) const
  {
  query_polygon_internal(vertex, 0, pixset);
  }

template<typename I> void T_Healpix_Base<I>::query_polygon_inclusive
  (const vector<pointing> &vertex, rangeset<I> &pixset, int fact) const
  {
  MR_assert(fact>0,"fact must be a positive integer");
  if ((sizeof(I)<8) && (((I(1)<<order_max)/nside_)<fact))
    {
    T_Healpix_Base<int64_t> base2(nside_,scheme_,SET_NSIDE);
    base2.query_polygon_internal(vertex,fact,pixset);
    return;
    }
  query_polygon_internal(vertex, fact, pixset);
  }

template<typename I> void T_Healpix_Base<I>::query_strip_internal
  (double theta1, double theta2, bool inclusive, rangeset<I> &pixset) const
  {
  if (scheme_==RING)
    {
    I ring1 = max(I(1),1+ring_above(cos(theta1))),
      ring2 = min(4*nside_-1,ring_above(cos(theta2)));
    if (inclusive)
      {
      ring1 = max(I(1),ring1-1);
      ring2 = min(4*nside_-1,ring2+1);
      }

    I sp1,rp1,sp2,rp2;
    bool dummy;
    get_ring_info_small(ring1,sp1,rp1,dummy);
    get_ring_info_small(ring2,sp2,rp2,dummy);
    I pix1 = sp1,
      pix2 = sp2+rp2;
    if (pix1<=pix2) pixset.append(pix1,pix2);
    }
  else
    MR_fail("query_strip not yet implemented for NESTED");
  }

template<typename I> void T_Healpix_Base<I>::query_strip (double theta1,
  double theta2, bool inclusive, rangeset<I> &pixset) const
  {
  pixset.clear();

  if (theta1<theta2)
    query_strip_internal(theta1,theta2,inclusive,pixset);
  else
    {
    query_strip_internal(0.,theta2,inclusive,pixset);
    rangeset<I> ps2;
    query_strip_internal(theta1,pi,inclusive,ps2);
    pixset.append(ps2);
    }
  }

template<typename I> inline void T_Healpix_Base<I>::get_ring_info_small
  (I ring, I &startpix, I &ringpix, bool &shifted) const
  {
  if (ring < nside_)
    {
    shifted = true;
    ringpix = 4*ring;
    startpix = 2*ring*(ring-1);
    }
  else if (ring < 3*nside_)
    {
    shifted = ((ring-nside_) & 1) == 0;
    ringpix = 4*nside_;
    startpix = ncap_ + (ring-nside_)*ringpix;
    }
  else
    {
    shifted = true;
    I nr= 4*nside_-ring;
    ringpix = 4*nr;
    startpix = npix_-2*nr*(nr+1);
    }
  }

template<typename I> void T_Healpix_Base<I>::get_ring_info (I ring, I &startpix,
  I &ringpix, double &costheta, double &sintheta, bool &shifted) const
  {
  I northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    costheta = 1 - tmp;
    sintheta = sqrt(tmp*(2-tmp));
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    costheta = (2*nside_-northring)*fact1_;
    sintheta = sqrt((1+costheta)*(1-costheta));
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    costheta = -costheta;
    startpix = npix_ - startpix - ringpix;
    }
  }
template<typename I> void T_Healpix_Base<I>::get_ring_info2 (I ring,
  I &startpix, I &ringpix, double &theta, bool &shifted) const
  {
  I northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    double costheta = 1 - tmp;
    double sintheta = sqrt(tmp*(2-tmp));
    theta = atan2(sintheta,costheta);
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    theta = acos((2*nside_-northring)*fact1_);
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    theta = pi-theta;
    startpix = npix_ - startpix - ringpix;
    }
  }

template<typename I> void T_Healpix_Base<I>::neighbors (I pix,
  array<I,8> &result) const
  {
  int ix, iy, face_num;
  (scheme_==RING) ?
    ring2xyf(pix,ix,iy,face_num) : nest2xyf(pix,ix,iy,face_num);

  const I nsm1 = nside_-1;
  if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
    {
    if (scheme_==RING)
      for (size_t m=0; m<8; ++m)
        result[m] = xyf2ring(ix+nb_xoffset[m],iy+nb_yoffset[m],face_num);
    else
      {
      I fpix = I(face_num)<<(2*order_),
        px0=spread<I>(ix  ), py0=spread<I>(iy  )<<1,
        pxp=spread<I>(ix+1), pyp=spread<I>(iy+1)<<1,
        pxm=spread<I>(ix-1), pym=spread<I>(iy-1)<<1;

      result[0] = fpix+pxm+py0; result[1] = fpix+pxm+pyp;
      result[2] = fpix+px0+pyp; result[3] = fpix+pxp+pyp;
      result[4] = fpix+pxp+py0; result[5] = fpix+pxp+pym;
      result[6] = fpix+px0+pym; result[7] = fpix+pxm+pym;
      }
    }
  else
    {
    for (size_t i=0; i<8; ++i)
      {
      int x=ix+nb_xoffset[i], y=iy+nb_yoffset[i];
      int nbnum=4;
      if (x<0)
        { x+=nside_; nbnum-=1; }
      else if (x>=nside_)
        { x-=nside_; nbnum+=1; }
      if (y<0)
        { y+=nside_; nbnum-=3; }
      else if (y>=nside_)
        { y-=nside_; nbnum+=3; }

      int f = nb_facearray[nbnum][face_num];
      if (f>=0)
        {
        int bits = nb_swaparray[nbnum][face_num>>2];
        if (bits&1) x=nside_-x-1;
        if (bits&2) y=nside_-y-1;
        if (bits&4) std::swap(x,y);
        result[i] = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
        }
      else
        result[i] = -1;
      }
    }
  }

template<typename I> void T_Healpix_Base<I>::get_interpol (const pointing &ptg,
  array<I,4> &pix, array<double,4> &wgt) const
  {
  MR_assert((ptg.theta>=0)&&(ptg.theta<=pi),"invalid theta value");
  double z = cos (ptg.theta);
  I ir1 = ring_above(z);
  I ir2 = ir1+1;
  double theta1, theta2, w1, tmp, dphi;
  I sp,nr;
  bool shift;
  I i1,i2;
  if (ir1>0)
    {
    get_ring_info2 (ir1, sp, nr, theta1, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? I(tmp)-1 : I(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[0] = sp+i1; pix[1] = sp+i2;
    wgt[0] = 1-w1; wgt[1] = w1;
    }
  if (ir2<(4*nside_))
    {
    get_ring_info2 (ir2, sp, nr, theta2, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? I(tmp)-1 : I(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[2] = sp+i1; pix[3] = sp+i2;
    wgt[2] = 1-w1; wgt[3] = w1;
    }

  if (ir1==0)
    {
    double wtheta = ptg.theta/theta2;
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    double fac = (1-wtheta)*0.25;
    wgt[0] = fac; wgt[1] = fac; wgt[2] += fac; wgt[3] +=fac;
    pix[0] = (pix[2]+2)&3;
    pix[1] = (pix[3]+2)&3;
    }
  else if (ir2==4*nside_)
    {
    double wtheta = (ptg.theta-theta1)/(pi-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    double fac = wtheta*0.25;
    wgt[0] += fac; wgt[1] += fac; wgt[2] = fac; wgt[3] =fac;
    pix[2] = ((pix[0]+2)&3)+npix_-4;
    pix[3] = ((pix[1]+2)&3)+npix_-4;
    }
  else
    {
    double wtheta = (ptg.theta-theta1)/(theta2-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    }

  if (scheme_==NEST)
    for (size_t m=0; m<pix.size(); ++m)
      pix[m] = ring2nest(pix[m]);
  }

template<typename I> void T_Healpix_Base<I>::swap (T_Healpix_Base &other)
  {
  std::swap(order_,other.order_);
  std::swap(nside_,other.nside_);
  std::swap(npface_,other.npface_);
  std::swap(ncap_,other.ncap_);
  std::swap(npix_,other.npix_);
  std::swap(fact1_,other.fact1_);
  std::swap(fact2_,other.fact2_);
  std::swap(scheme_,other.scheme_);
  }

template<typename I> double T_Healpix_Base<I>::max_pixrad() const
  {
  vec3 va,vb;
  va.set_z_phi (2./3., pi/(4*nside_));
  double t1 = 1.-1./nside_;
  t1*=t1;
  vb.set_z_phi (1-t1/3, 0);
  return v_angle(va,vb);
  }

template<typename I> double T_Healpix_Base<I>::max_pixrad(I ring) const
  {
  if (ring>=2*nside_) ring=4*nside_-ring;
  double z=ring2z(ring), z_up=ring2z(ring-1);
  vec3 mypos, uppos;
  uppos.set_z_phi(z_up,0);
  if (ring<=nside_)
    {
    mypos.set_z_phi(z,pi/(4*ring));
    double v1=v_angle(mypos,uppos);
    if (ring!=1) return v1;
    uppos.set_z_phi(ring2z(ring+1),pi/(4*(min(nside_,ring+1))));
    return max(v1,v_angle(mypos,uppos));
    }
  mypos.set_z_phi(z,0);
  double vdist=v_angle(mypos,uppos);
  double hdist=sqrt(1.-z*z)*pi/(4*nside_);
  return max(hdist,vdist);
  }

template<typename I> void T_Healpix_Base<I>::xyf2loc (double x, double y,
  int face, double &z, double &phi, double &sth, bool &have_sth) const
  {
  have_sth = false;
  double jr = jrll[face] - x - y;
  double nr;
  if (jr<1)
    {
    nr = jr;
    double tmp = nr*nr/3.;
    z = 1 - tmp;
    if (z > 0.99)
      {
      sth = std::sqrt(tmp*(2.0-tmp));
      have_sth = true;
      }
    }
  else if (jr>3)
    {
    nr = 4-jr;
    double tmp = nr*nr/3.;
    z = tmp - 1;
    if (z<-0.99)
      {
      sth = std::sqrt(tmp*(2.-tmp));
      have_sth = true;
      }
    }
  else
    {
    nr = 1;
    z = (2-jr)*2./3.;
    }

  double tmp=jpll[face]*nr+x-y;
  if (tmp<0) tmp+=8;
  if (tmp>=8) tmp-=8;
  phi = (nr<1e-15) ? 0 : (0.5*halfpi*tmp)/nr;
  }

namespace {

vec3 locToVec3 (double z, double phi, double sth, bool have_sth)
  {
  if (have_sth)
    return vec3(sth*cos(phi),sth*sin(phi),z);
  else
    {
    vec3 res;
    res.set_z_phi (z, phi);
    return res;
    }
  }

} // unnamed namespace

template<typename I> void T_Healpix_Base<I>::boundaries(I pix, size_t step,
  vector<vec3> &out) const
  {
  out.resize(4*step);
  int ix, iy, face;
  pix2xyf(pix, ix, iy, face);
  double dc = 0.5 / nside_;
  double xc = (ix + 0.5)/nside_, yc = (iy + 0.5)/nside_;
  double d = 1.0/(step*nside_);
  for (size_t i=0; i<step; ++i)
    {
    double z, phi, sth;
    bool have_sth;
    xyf2loc(xc+dc-i*d, yc+dc, face, z, phi, sth, have_sth);
    out[i] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc-dc, yc+dc-i*d, face, z, phi, sth, have_sth);
    out[i+step] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc-dc+i*d, yc-dc, face, z, phi, sth, have_sth);
    out[i+2*step] = locToVec3(z, phi, sth, have_sth);
    xyf2loc(xc+dc, yc-dc+i*d, face, z, phi, sth, have_sth);
    out[i+3*step] = locToVec3(z, phi, sth, have_sth);
    }
  }

template<typename I> vector<int> T_Healpix_Base<I>::swap_cycles() const
  {
  MR_assert(order_>=0, "need hierarchical map");
  MR_assert(order_<=13, "map too large");
  vector<int> result(swap_clen[order_]);
  size_t ofs=0;
  for (int m=0; m<order_;++m) ofs+=swap_clen[m];
  for (size_t m=0; m<result.size();++m) result[m]=swap_cycle[m+ofs];
  return result;
  }

template class T_Healpix_Base<int>;
template class T_Healpix_Base<int64_t>;

}}

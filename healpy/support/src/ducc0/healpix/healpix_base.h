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

/*! \file healpix_base.h
 *
 *  \copyright Copyright (C) 2003-2020 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_BASE_H
#define HEALPIX_BASE_H

#include <array>
#include <cmath>
#include <vector>
#include "ducc0/infra/error_handling.h"
#include "ducc0/math/math_utils.h"
#include "ducc0/math/vec3.h"
#include "ducc0/healpix/healpix_tables.h"
#include "ducc0/math/pointing.h"
#include "ducc0/math/rangeset.h"

namespace ducc0 {

namespace detail_healpix {

template<typename I> struct Orderhelper__ {};
template<> struct Orderhelper__<int> {enum{omax=13};};
template<> struct Orderhelper__<int64_t> {enum{omax=29};};

/*! Functionality related to the HEALPix pixelisation. */
template<typename I> class T_Healpix_Base: public Healpix_Tables
  {
  protected:
    /*! The order of the map; -1 for nonhierarchical map. */
    int order_;
    /*! The N_side parameter of the map; 0 if not allocated. */
    I nside_;
    I npface_, ncap_, npix_;
    double fact1_, fact2_;
    /*! The map's ordering scheme. */
    Ordering_Scheme scheme_;

    /*! Returns the number of the next ring to the north of \a z=cos(theta).
        It may return 0; in this case \a z lies north of all rings. */
    inline I ring_above (double z) const;
    void in_ring (I iz, double phi0, double dphi, rangeset<I> &pixset) const;

    template<typename I2> void query_multidisc (const std::vector<vec3> &norm,
      const std::vector<double> &rad, int fact, rangeset<I2> &pixset) const;

    void query_multidisc_general (const std::vector<vec3> &norm, const std::vector<double> &rad,
      bool inclusive, const std::vector<int> &cmds, rangeset<I> &pixset) const;

    void query_strip_internal (double theta1, double theta2, bool inclusive,
      rangeset<I> &pixset) const;

    I xyf2nest(int ix, int iy, int face_num) const;
    void nest2xyf(I pix, int &ix, int &iy, int &face_num) const;
    I xyf2ring(int ix, int iy, int face_num) const;
    void ring2xyf(I pix, int &ix, int &iy, int &face_num) const;

    I loc2pix (double z, double phi, double sth, bool have_sth) const;
    void pix2loc (I pix, double &z, double &phi, double &sth, bool &have_sth)
      const;

    void xyf2loc(double x, double y, int face, double &z, double &ph,
      double &sth, bool &have_sth) const;

    I nest_peano_helper (I pix, int dir) const;

    typedef I (T_Healpix_Base::*swapfunc)(I pix) const;
//    using swapfunc = I(I) const;

  public:
    enum {order_max=Orderhelper__<I>::omax};

    /*! Calculates the map order from its \a N_side parameter.
        Returns -1 if \a nside is not a power of 2.
        \param nside the \a N_side parameter */
    static int nside2order (I nside);
    /*! Calculates the \a N_side parameter from the number of pixels.
        \param npix the number of pixels */
    static I npix2nside (I npix);
    /*! Constructs an unallocated object. */
    T_Healpix_Base ();
    /*! Constructs an object with a given \a order and the ordering
        scheme \a scheme. */
    T_Healpix_Base (int order, Ordering_Scheme scheme)
      { Set (order, scheme); }
    /*! Constructs an object with a given \a nside and the ordering
        scheme \a scheme. The \a nside_dummy parameter must be set to
        SET_NSIDE. */
    T_Healpix_Base (I nside, Ordering_Scheme scheme, const nside_dummy)
      { SetNside (nside, scheme); }

    /*! Adjusts the object to \a order and \a scheme. */
    void Set (int order, Ordering_Scheme scheme);
    /*! Adjusts the object to \a nside and \a scheme. */
    void SetNside (I nside, Ordering_Scheme scheme);

    /*! Returns the z-coordinate of the ring \a ring. This also works
        for the (not really existing) rings 0 and 4*nside. */
    double ring2z (I ring) const;
    /*! Returns the number of the ring in which \a pix lies. */
    I pix2ring (I pix) const;

    I xyf2pix(int ix, int iy, int face_num) const
      {
      return (scheme_==RING) ?
        xyf2ring(ix,iy,face_num) : xyf2nest(ix,iy,face_num);
      }
    void pix2xyf(I pix, int &ix, int &iy, int &face_num) const
      {
      (scheme_==RING) ?
        ring2xyf(pix,ix,iy,face_num) : nest2xyf(pix,ix,iy,face_num);
      }

    /*! Translates a pixel number from NEST to RING. */
    I nest2ring (I pix) const;
    /*! Translates a pixel number from RING to NEST. */
    I ring2nest (I pix) const;
    /*! Translates a pixel number from NEST to its Peano index. */
    I nest2peano (I pix) const;
    /*! Translates a pixel number from its Peano index to NEST. */
    I peano2nest (I pix) const;

    /*! Returns the number of the pixel which contains the angular coordinates
        (\a z:=cos(theta), \a phi).
        \note This method is inaccurate near the poles at high resolutions. */
    I zphi2pix (double z, double phi) const
      { return loc2pix(z,phi,0.,false); }

    /*! Returns the number of the pixel which contains the angular coordinates
        \a ang. */
    I ang2pix (const pointing &ang) const
      {
      const double pi_=3.141592653589793238462643383279502884197;
      MR_assert((ang.theta>=0)&&(ang.theta<=pi_),"invalid theta value");
      return ((ang.theta<0.01) || (ang.theta > 3.14159-0.01)) ?
        loc2pix(cos(ang.theta),ang.phi,sin(ang.theta),true) :
        loc2pix(cos(ang.theta),ang.phi,0.,false);
      }
    /*! Returns the number of the pixel which contains the vector \a vec
        (\a vec is normalized if necessary). */
    I vec2pix (const vec3 &vec) const
      {
      double xl = 1./vec.Length();
      double phi = safe_atan2(vec.y,vec.x);
      double nz = vec.z*xl;
      if (std::abs(nz)>0.99)
        return loc2pix (nz,phi,std::sqrt(vec.x*vec.x+vec.y*vec.y)*xl,true);
      else
        return loc2pix (nz,phi,0,false);
      }

    /*! Returns the angular coordinates (\a z:=cos(theta), \a phi) of the center
        of the pixel with number \a pix.
        \note This method is inaccurate near the poles at high resolutions. */
    void pix2zphi (I pix, double &z, double &phi) const
      {
      bool dum_b;
      double dum_d;
      pix2loc(pix,z,phi,dum_d,dum_b);
      }

    /*! Returns the angular coordinates of the center of the pixel with
        number \a pix. */
    pointing pix2ang (I pix) const
      {
      double z, phi, sth;
      bool have_sth;
      pix2loc (pix,z,phi,sth,have_sth);
      return have_sth ? pointing(atan2(sth,z),phi) : pointing(acos(z),phi);
      }
    /*! Returns the vector to the center of the pixel with number \a pix. */
    vec3 pix2vec (I pix) const
      {
      double z, phi, sth;
      bool have_sth;
      pix2loc (pix,z,phi,sth,have_sth);
      if (have_sth)
        return vec3(sth*cos(phi),sth*sin(phi),z);
      else
        {
        vec3 res;
        res.set_z_phi (z, phi);
        return res;
        }
      }
    /*! Returns the pixel number for this T_Healpix_Base corresponding to the
        pixel number \a pix in \a b.
        \note \a b.Nside()\%Nside() must be 0. */
    I pixel_import (I pix, const T_Healpix_Base &b) const
      {
      I ratio = b.nside_/nside_;
      MR_assert(nside_*ratio==b.nside_,"bad nside ratio");
      int x, y, f;
      b.pix2xyf(pix, x, y, f);
      x/=ratio; y/=ratio;
      return xyf2pix(x, y, f);
      }

    template<typename I2> void query_disc_internal (pointing ptg, double radius,
      int fact, rangeset<I2> &pixset) const;

    /*! Returns the range set of all pixels whose centers lie within the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param pixset a \a rangeset object containing the indices of all pixels
           whose centers lie inside the disk
        \note This method is more efficient in the RING scheme. */
    void query_disc (pointing ptg, double radius, rangeset<I> &pixset) const;
    /*! Returns the range set of all pixels whose centers lie within the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \note This method is more efficient in the RING scheme. */
    rangeset<I> query_disc (pointing ptg, double radius) const
      {
      rangeset<I> res;
      query_disc(ptg, radius, res);
      return res;
      }
    /*! Returns the range set of all pixels which overlap with the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param pixset a \a rangeset object containing the indices of all pixels
           overlapping with the disk.
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the disk at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    void query_disc_inclusive (pointing ptg, double radius, rangeset<I> &pixset,
      int fact=1) const;
    /*! Returns the range set of all pixels which overlap with the disk
        defined by \a dir and \a radius.
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the disk at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    rangeset<I> query_disc_inclusive (pointing ptg, double radius,
      int fact=1) const
      {
      rangeset<I> res;
      query_disc_inclusive(ptg, radius, res, fact);
      return res;
      }

    /*! \deprecated Please use the version based on \a rangeset */
    void query_disc (const pointing &dir, double radius,
      std::vector<I> &listpix) const
      {
      rangeset<I> pixset;
      query_disc(dir,radius,pixset);
      pixset.toVector(listpix);
      }
    /*! \deprecated Please use the version based on \a rangeset */
    void query_disc_inclusive (const pointing &dir, double radius,
      std::vector<I> &listpix, int fact=1) const
      {
      rangeset<I> pixset;
      query_disc_inclusive(dir,radius,pixset,fact);
      pixset.toVector(listpix);
      }

    template<typename I2> void query_polygon_internal
      (const std::vector<pointing> &vertex, int fact,
      rangeset<I2> &pixset) const;

    /*! Returns a range set of pixels whose centers lie within the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \param pixset a \a rangeset object containing the indices of all pixels
           whose centers lie inside the polygon
        \note This method is more efficient in the RING scheme. */
    void query_polygon (const std::vector<pointing> &vertex,
      rangeset<I> &pixset) const;
    /*! Returns a range set of pixels whose centers lie within the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \note This method is more efficient in the RING scheme. */
    rangeset<I> query_polygon (const std::vector<pointing> &vertex) const
      {
      rangeset<I> res;
      query_polygon(vertex, res);
      return res;
      }

    /*! Returns a range set of pixels which overlap with the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \param pixset a \a rangeset object containing the indices of all pixels
           overlapping with the polygon.
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the polygon at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    void query_polygon_inclusive (const std::vector<pointing> &vertex,
      rangeset<I> &pixset, int fact=1) const;
    /*! Returns a range set of pixels which overlap with the convex
        polygon defined by the \a vertex array.
        \param vertex array containing the vertices of the polygon.
        \param fact The overlapping test will be done at the resolution
           \a fact*nside. For NESTED ordering, \a fact must be a power of 2,
           else it can be any positive integer. A typical choice would be 4.
        \note This method may return some pixels which don't overlap with
           the polygon at all. The higher \a fact is chosen, the fewer false
           positives are returned, at the cost of increased run time.
        \note This method is more efficient in the RING scheme. */
    rangeset<I> query_polygon_inclusive (const std::vector<pointing> &vertex,
      int fact=1) const
      {
      rangeset<I> res;
      query_polygon_inclusive(vertex, res, fact);
      return res;
      }

    /*! Returns a range set of pixels whose centers lie within the colatitude
        range defined by \a theta1 and \a theta2 (if \a inclusive==false), or
        which overlap with this region (if \a inclusive==true). If
        \a theta1<theta2, the region between both angles is considered,
        otherwise the regions \a 0<theta<theta2 and \a theta1<theta<pi.
        \param theta1 first colatitude
        \param theta2 second colatitude
        \param inclusive if \a false, return the exact set of pixels whose
           pixels centers lie within the region; if \a true, return all pixels
           that overlap with the region. */
    void query_strip (double theta1, double theta2, bool inclusive,
      rangeset<I> &pixset) const;
    /*! Returns a range set of pixels whose centers lie within the colatitude
        range defined by \a theta1 and \a theta2 (if \a inclusive==false), or
        which overlap with this region (if \a inclusive==true). If
        \a theta1<theta2, the region between both angles is considered,
        otherwise the regions \a 0<theta<theta2 and \a theta1<theta<pi.
        \param theta1 first colatitude
        \param theta2 second colatitude
        \param inclusive if \a false, return the exact set of pixels whose
           pixels centers lie within the region; if \a true, return all pixels
           that overlap with the region. */
    rangeset<I> query_strip (double theta1, double theta2, bool inclusive) const
      {
      rangeset<I> res;
      query_strip(theta1, theta2, inclusive, res);
      return res;
      }

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param costheta the cosine of the colatitude of the ring
        \param sintheta the sine of the colatitude of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info (I ring, I &startpix, I &ringpix,
      double &costheta, double &sintheta, bool &shifted) const;
    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param theta the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info2 (I ring, I &startpix, I &ringpix,
      double &theta, bool &shifted) const;
    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info_small (I ring, I &startpix, I &ringpix,
        bool &shifted) const;
    /*! Returns the neighboring pixels of \a pix in \a result.
        On exit, \a result contains (in this order)
        the pixel numbers of the SW, W, NW, N, NE, E, SE and S neighbor
        of \a pix. If a neighbor does not exist (this can only be the case
        for the W, N, E and S neighbors), its entry is set to -1.

        \note This method works in both RING and NEST schemes, but is
          considerably faster in the NEST scheme. */
    void neighbors (I pix, std::array<I,8> &result) const;
    /*! Returns interpolation information for the direction \a ptg.
        The surrounding pixels are returned in \a pix, their corresponding
        weights in \a wgt.
        \note This method works in both RING and NEST schemes, but is
          considerably faster in the RING scheme. */
    void get_interpol (const pointing &ptg, std::array<I,4> &pix,
                       std::array<double,4> &wgt) const;

    /*! Returns the order parameter of the object. */
    int Order() const { return order_; }
    /*! Returns the \a N_side parameter of the object. */
    I Nside() const { return nside_; }
    /*! Returns the number of pixels of the object. */
    I Npix() const { return npix_; }
    /*! Returns the ordering scheme of the object. */
    Ordering_Scheme Scheme() const { return scheme_; }

    /*! Returns \a true, if both objects have the same nside and scheme,
        else \a false. */
    bool conformable (const T_Healpix_Base &other) const
      { return ((nside_==other.nside_) && (scheme_==other.scheme_)); }

    /*! Swaps the contents of two Healpix_Base objects. */
    void swap (T_Healpix_Base &other);

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners. */
    double max_pixrad() const;

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners in a given ring. */
    double max_pixrad(I ring) const;

    /*! Returns a set of points along the boundary of the given pixel.
        \a step=1 gives 4 points on the corners. The first point corresponds
        to the northernmost corner, the subsequent points follow the pixel
        boundary through west, south and east corners.
        \param pix pixel index number
        \param step the number of returned points is 4*step. */
    void boundaries (I pix, size_t step, std::vector<vec3> &out) const;

    std::vector<int> swap_cycles() const;
  };

/*! T_Healpix_Base for Nside up to 2^13. */
typedef T_Healpix_Base<int> Healpix_Base;
/*! T_Healpix_Base for Nside up to 2^29. */
typedef T_Healpix_Base<int64_t> Healpix_Base2;

}

using detail_healpix::Healpix_Base;
using detail_healpix::Healpix_Base2;

}

#endif

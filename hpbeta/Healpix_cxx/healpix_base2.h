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

/*! \file healpix_base2.h
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_BASE2_H
#define HEALPIX_BASE2_H

#include "healpix_base.h"
#include "datatypes.h"

/*! Functionality related to the HEALPix pixelisation. Supports resolutions up to
    N_side = 2^29. */
class Healpix_Base2
  {
  protected:
    enum { order_max=29 };

    class Tablefiller
      {
      public:
        Tablefiller();
      };
    static Tablefiller Filler;
    friend class Tablefiller;

    static short ctab[0x100], utab[0x100];

    static const int jrll[];
    static const int jpll[];

    /*! The order of the map; -1 for nonhierarchical map. */
    int order_;
    /*! The N_side parameter of the map; 0 if not allocated. */
    int64 nside_;
    int64 npface_, ncap_, npix_;
    double fact1_, fact2_;
    /*! The map's ordering scheme. */
    Healpix_Ordering_Scheme scheme_;

    inline int64 ring_above (double z) const;
    void in_ring (int64 iz, double phi0, double dphi,
      std::vector<int64> &listir) const;

    int64 xyf2nest(int ix, int iy, int face_num) const;
    void nest2xyf(int64 pix, int &ix, int &iy, int &face_num) const;
    int64 xyf2ring(int ix, int iy, int face_num) const;
    void ring2xyf(int64 pix, int &ix, int &iy, int &face_num) const;

    typedef int64 (Healpix_Base2::*swapfunc)(int64 pix) const;
    typedef void (Healpix_Base2::*pix2xyf)
                 (int64 pix, int &x, int &y, int &f) const;
    typedef int64 (Healpix_Base2::*xyf2pix) (int x, int y, int f) const;

  public:
    /*! Calculates the map order from its \a N_side parameter.
        Returns -1 if \a nside is not a power of 2.
        \param nside the \a N_side parameter */
    static int nside2order (int64 nside)
      {
      planck_assert (nside>0, "invalid value for Nside");
      if ((nside)&(nside-1)) return -1;
      return ilog2(nside);
      }
    /*! Calculates the \a N_side parameter from the number of pixels.
        \param npix the number of pixels */
    static int64 npix2nside (int64 npix);
    /*! Constructs an unallocated object. */
    Healpix_Base2 ()
      : order_(-1), nside_(0), npface_(0), ncap_(0), npix_(0),
        fact1_(0), fact2_(0), scheme_(RING) {}
    /*! Constructs an object with a given \a order and the ordering
        scheme \a scheme. */
    Healpix_Base2 (int order, Healpix_Ordering_Scheme scheme)
      { Set (order, scheme); }
    /*! Constructs an object with a given \a nside and the ordering
        scheme \a scheme. The \a nside_dummy parameter must be set to
        SET_NSIDE. */
    Healpix_Base2 (int64 nside, Healpix_Ordering_Scheme scheme,
        const nside_dummy)
      { SetNside (nside, scheme); }

    /* Adjusts the object to \a order and \a scheme. */
    void Set (int order, Healpix_Ordering_Scheme scheme)
      {
      planck_assert ((order>=0)&&(order<=order_max), "bad order");
      order_  = order;
      nside_  = int64(1)<<order;
      npface_ = nside_<<order_;
      ncap_   = (npface_-nside_)<<1;
      npix_   = 12*npface_;
      fact2_  = 4./npix_;
      fact1_  = (nside_<<1)*fact2_;
      scheme_ = scheme;
      }
    /* Adjusts the object to \a nside and \a scheme. */
    void SetNside (int64 nside, Healpix_Ordering_Scheme scheme)
      {
      order_  = nside2order(nside);
      planck_assert ((scheme!=NEST) || (order_>=0),
        "SetNside: nside must be power of 2 for nested maps");
      nside_  = nside;
      npface_ = nside_*nside_;
      ncap_   = (npface_-nside_)<<1;
      npix_   = 12*npface_;
      fact2_  = 4./npix_;
      fact1_  = (nside_<<1)*fact2_;
      scheme_ = scheme;
      }

    /*! Returns the z-coordinate of the ring \a ring. This also works
        for the (not really existing) rings 0 and 4*nside. */
    double ring2z (int64 ring) const;
    /*! Returns the number of the ring in which \a pix lies. */
    int64 pix2ring (int64 pix) const;

    /*! Translates a pixel number from NEST to RING. */
    int64 nest2ring (int64 pix) const;
    /*! Translates a pixel number from RING to NEST. */
    int64 ring2nest (int64 pix) const;
    /*! Translates a pixel number from NEST to its Peano index. */
    int64 nest2peano (int64 pix) const;
    /*! Translates a pixel number from its Peano index to NEST. */
    int64 peano2nest (int64 pix) const;

    int64 ang2pix_z_phi (double z, double phi) const;

    /*! Returns the number of the pixel which contains the angular coordinates
        \a ang. */
    int64 ang2pix (const pointing &ang) const
      { return ang2pix_z_phi (cos(ang.theta), ang.phi); }
    /*! Returns the number of the pixel which contains the vector \a vec
        (\a vec is normalized if necessary). */
    int64 vec2pix (const vec3 &vec) const
      { return ang2pix_z_phi (vec.z/vec.Length(), safe_atan2(vec.y,vec.x)); }

    void pix2ang_z_phi (int64 pix, double &z, double &phi) const;

    /*! Returns the angular coordinates of the center of the pixel with
        number \a pix. */
    pointing pix2ang (int64 pix) const
      {
      double z, phi;
      pix2ang_z_phi (pix,z,phi);
      return pointing(acos(z),phi);
      }
    /*! Returns the vector to the center of the pixel with number \a pix. */
    vec3 pix2vec (int64 pix) const
      {
      double z, phi;
      pix2ang_z_phi (pix,z,phi);
      vec3 res;
      res.set_z_phi (z, phi);
      return res;
      }

    /*! Returns the numbers of all pixels whose centers lie within \a radius
        of \a dir in \a listpix.
        \param dir the angular coordinates of the disc center
        \param radius the radius (in radians) of the disc
        \param listpix a vector containing the numbers of all pixels within
               the disc
        \note This method is more efficient in the RING scheme. */
    void query_disc (const pointing &dir, double radius,
      std::vector<int64> &listpix) const;
    /*! Returns the numbers of all pixels that lie at least partially within
        \a radius of \a dir in \a listpix. It may also return a few pixels
        which do not lie in the disk at all.
        \param dir the angular coordinates of the disc center
        \param radius the radius (in radians) of the disc
        \param listpix a vector containing the numbers of all pixels within
               the disc
        \note This method works in both RING and NEST schemes, but is
          considerably faster in the RING scheme. */
//FIXME: factor 1.362 not OK!
    void query_disc_inclusive (const pointing &dir, double radius,
      std::vector<int64> &listpix) const
        { query_disc (dir,radius+1.362*pi/(4*nside_),listpix); }

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
        \param ringpix the number of pixels in the ring
        \param costheta the cosine of the colatitude (in radians) of the ring
        \param sintheta the sine of the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info (int64 ring, int64 &startpix, int64 &ringpix,
      double &costheta, double &sintheta, bool &shifted) const;
    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
        \param ringpix the number of pixels in the ring
        \param theta the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info2 (int64 ring, int64 &startpix, int64 &ringpix,
      double &theta, bool &shifted) const;

    /*! Returns the neighboring pixels of \a pix in \a result.
        On exit, \a result contains (in this order)
        the pixel numbers of the SW, W, NW, N, NE, E, SE and S neighbor
        of \a pix. If a neighbor does not exist (this can only be the case
        for the W, N, E and S neighbors), its entry is set to -1.

        \note This method works in both RING and NEST schemes, but is
          considerably faster in the NEST scheme. */
    void neighbors (int64 pix, fix_arr<int64,8> &result) const;
    /*! Returns interpolation information for the direction \a ptg.
        The surrounding pixels are returned in \a pix, their corresponding
        weights in \a wgt.
        \note This method works in both RING and NEST schemes, but is
          considerably faster in the RING scheme. */
    void get_interpol (const pointing &ptg, fix_arr<int64,4> &pix,
                       fix_arr<double,4> &wgt) const;

    /*! Returns the order parameter of the object. */
    int Order() const { return order_; }
    /*! Returns the \a N_side parameter of the object. */
    int64 Nside() const { return nside_; }
    /*! Returns the number of pixels of the object. */
    int64 Npix() const { return npix_; }
    /*! Returns the ordering scheme of the object. */
    Healpix_Ordering_Scheme Scheme() const { return scheme_; }

    /*! Returns \a true, if both objects have the same nside and scheme,
        else  \a false. */
    bool conformable (const Healpix_Base2 &other) const
      { return ((nside_==other.nside_) && (scheme_==other.scheme_)); }

    /*! Swaps the contents of two Healpix_Base objects. */
    void swap (Healpix_Base2 &other)
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

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners. */
    double max_pixrad() const;
  };

#endif

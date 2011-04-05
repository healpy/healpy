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

/*! \file healpix_base.h
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_BASE_H
#define HEALPIX_BASE_H

#include <vector>
#include "cxxutils.h"
#include "pointing.h"
#include "arr.h"
#include "rangeset.h"

/*! The two possible ordering schemes of a HEALPix map. */
enum Healpix_Ordering_Scheme { RING, /*!< RING scheme */
                               NEST  /*!< NESTED scheme */
                             };

Healpix_Ordering_Scheme string2HealpixScheme (const std::string &inp);

class nside_dummy {};
extern const nside_dummy SET_NSIDE;

/*! Functionality related to the HEALPix pixelisation. */
class Healpix_Base
  {
  protected:
    enum { order_max=13 };

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
    int nside_;
    int npface_, ncap_, npix_;
    double fact1_, fact2_;
    /*! The map's ordering scheme. */
    Healpix_Ordering_Scheme scheme_;

    /*! Returns the number of the next ring to the north of \a z=cos(theta).
        It may return 0; in this case \a z lies north of all rings. */
    inline int ring_above (double z) const;
    void in_ring (int iz, double phi0, double dphi, rangeset<int> &pixset)
      const;

    void query_multidisc (const arr<vec3> &norm, const arr<double> &rad,
      bool inclusive, rangeset<int> &pixset) const;

    int xyf2nest(int ix, int iy, int face_num) const;
    void nest2xyf(int pix, int &ix, int &iy, int &face_num) const;
    int xyf2ring(int ix, int iy, int face_num) const;
    void ring2xyf(int pix, int &ix, int &iy, int &face_num) const;

    typedef int (Healpix_Base::*swapfunc)(int pix) const;
    typedef void (Healpix_Base::*pix2xyf)
                 (int pix, int &x, int &y, int &f) const;
    typedef int (Healpix_Base::*xyf2pix) (int x, int y, int f) const;

  public:
    /*! Calculates the map order from its \a N_side parameter.
        Returns -1 if \a nside is not a power of 2.
        \param nside the \a N_side parameter */
    static int nside2order (int nside);
    /*! Calculates the \a N_side parameter from the number of pixels.
        \param npix the number of pixels */
    static int npix2nside (int npix);
    /*! Constructs an unallocated object. */
    Healpix_Base ();
    /*! Constructs an object with a given \a order and the ordering
        scheme \a scheme. */
    Healpix_Base (int order, Healpix_Ordering_Scheme scheme)
      { Set (order, scheme); }
    /*! Constructs an object with a given \a nside and the ordering
        scheme \a scheme. The \a nside_dummy parameter must be set to
        SET_NSIDE. */
    Healpix_Base (int nside, Healpix_Ordering_Scheme scheme, const nside_dummy)
      { SetNside (nside, scheme); }

    /*! Adjusts the object to \a order and \a scheme. */
    void Set (int order, Healpix_Ordering_Scheme scheme);
    /*! Adjusts the object to \a nside and \a scheme. */
    void SetNside (int nside, Healpix_Ordering_Scheme scheme);

    /*! Returns the z-coordinate of the ring \a ring. This also works
        for the (not really existing) rings 0 and 4*nside. */
    double ring2z (int ring) const;
    /*! Returns the number of the ring in which \a pix lies. */
    int pix2ring (int pix) const;

    /*! Translates a pixel number from NEST to RING. */
    int nest2ring (int pix) const;
    /*! Translates a pixel number from RING to NEST. */
    int ring2nest (int pix) const;
    /*! Translates a pixel number from NEST to its Peano index. */
    int nest2peano (int pix) const;
    /*! Translates a pixel number from its Peano index to NEST. */
    int peano2nest (int pix) const;

    /*! Returns the number of the pixel which contains the angular coordinates
        (\a z:=cos(theta), \a phi). */
    int zphi2pix (double z, double phi) const;

    /*! Returns the number of the pixel which contains the angular coordinates
        \a ang. */
    int ang2pix (const pointing &ang) const
      { return zphi2pix (cos(ang.theta), ang.phi); }
    /*! Returns the number of the pixel which contains the vector \a vec
        (\a vec is normalized if necessary). */
    int vec2pix (const vec3 &vec) const
      { return zphi2pix (vec.z/vec.Length(), safe_atan2(vec.y,vec.x)); }

    /*! Returns the angular coordinates (\a z:=cos(theta), \a phi) of the center
        of the pixel with number \a pix. */
    void pix2zphi (int pix, double &z, double &phi) const;

    /*! Returns the angular coordinates of the center of the pixel with
        number \a pix. */
    pointing pix2ang (int pix) const
      {
      double z, phi;
      pix2zphi (pix,z,phi);
      return pointing(acos(z),phi);
      }
    /*! Returns the vector to the center of the pixel with number \a pix. */
    vec3 pix2vec (int pix) const
      {
      double z, phi;
      pix2zphi (pix,z,phi);
      vec3 res;
      res.set_z_phi (z, phi);
      return res;
      }

    /*! Returns a range set of pixels whose centers lie within the disk
        defined by \a dir and \a radius (if \a inclusive==false), or which
        overlap with this disk (if \a inclusive==true).
        \param dir the angular coordinates of the disk center
        \param radius the radius (in radians) of the disk
        \param inclusive if \a false, return the exact set of pixels whose
           pixels centers lie within the disk; if \a true, return all pixels
           that overlap with the disk, and maybe a few more.
        \param pixset a \a rangeset object containing the indices of all pixels
           within the disk
        \note This method is more efficient in the RING scheme, but the
           algorithm used for \a inclusive==true returns fewer false positives
           in the NEST scheme. */
    void query_disc (pointing ptg, double radius, bool inclusive,
      rangeset<int> &pixset) const;

    /*! \deprecated Please use the version based on \a rangeset */
    void query_disc (const pointing &dir, double radius,
      std::vector<int> &listpix) const
      {
      rangeset<int> pixset;
      query_disc(dir,radius,false,pixset);
      pixset.toVector(listpix);
      }
    /*! \deprecated Please use the version based on \a rangeset */
    void query_disc_inclusive (const pointing &dir, double radius,
      std::vector<int> &listpix) const
      {
      rangeset<int> pixset;
      query_disc(dir,radius,true,pixset);
      pixset.toVector(listpix);
      }

    /*! Returns a range set of pixels whose centers lie within the convex
        polygon defined by the \a vertex array (if \a inclusive==false), or
        which overlap with this polygon (if \a inclusive==true).
        \param vertex array containing the vertices of the polygon.
        \param inclusive if \a false, return the exact set of pixels whose
           pixels centers lie within the polygon; if \a true, return all pixels
           that overlap with the polygon, and maybe a few more.
        \note This method is currently only implemented in the NEST scheme. */
    void query_polygon (const std::vector<pointing> &vertex, bool inclusive,
      rangeset<int> &pixset) const;

    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param costheta the cosine of the colatitude of the ring
        \param sintheta the sine of the colatitude of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info (int ring, int &startpix, int &ringpix,
      double &costheta, double &sintheta, bool &shifted) const;
    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring
        \param theta the colatitude (in radians) of the ring
        \param shifted if \a true, the center of the first pixel is not at
               \a phi=0 */
    void get_ring_info2 (int ring, int &startpix, int &ringpix,
      double &theta, bool &shifted) const;
    /*! Returns useful information about a given ring of the map.
        \param ring the ring number (the number of the first ring is 1)
        \param startpix the number of the first pixel in the ring
               (NOTE: this is always given in the RING numbering scheme!)
        \param ringpix the number of pixels in the ring */
    void get_ring_info_small (int ring, int &startpix, int &ringpix) const;
    /*! Returns the neighboring pixels of \a pix in \a result.
        On exit, \a result contains (in this order)
        the pixel numbers of the SW, W, NW, N, NE, E, SE and S neighbor
        of \a pix. If a neighbor does not exist (this can only be the case
        for the W, N, E and S neighbors), its entry is set to -1.

        \note This method works in both RING and NEST schemes, but is
          considerably faster in the NEST scheme. */
    void neighbors (int pix, fix_arr<int,8> &result) const;
    /*! Returns interpolation information for the direction \a ptg.
        The surrounding pixels are returned in \a pix, their corresponding
        weights in \a wgt.
        \note This method works in both RING and NEST schemes, but is
          considerably faster in the RING scheme. */
    void get_interpol (const pointing &ptg, fix_arr<int,4> &pix,
                       fix_arr<double,4> &wgt) const;

    /*! Returns the order parameter of the object. */
    int Order() const { return order_; }
    /*! Returns the \a N_side parameter of the object. */
    int Nside() const { return nside_; }
    /*! Returns the number of pixels of the object. */
    int Npix() const { return npix_; }
    /*! Returns the ordering scheme of the object. */
    Healpix_Ordering_Scheme Scheme() const { return scheme_; }

    /*! Returns \a true, if both objects have the same nside and scheme,
        else  \a false. */
    bool conformable (const Healpix_Base &other) const
      { return ((nside_==other.nside_) && (scheme_==other.scheme_)); }

    /*! Swaps the contents of two Healpix_Base objects. */
    void swap (Healpix_Base &other);

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners. */
    double max_pixrad() const;

    /*! Returns the maximum angular distance (in radian) between any pixel
        center and its corners in a given ring. */
    double max_pixrad(int ring) const;

    arr<int> swap_cycles() const;
  };

#endif

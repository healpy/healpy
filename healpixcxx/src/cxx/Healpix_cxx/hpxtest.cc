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
 *  Copyright (C) 2004-2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

/*

Candidates for testing the validity of the Healpix routines:

- done: ang2pix(pix2ang(i)) = i for all Healpix_Bases
- done: pix2ang(ang2pix(ptg)) dot ptg > 1-delta for all Healpix_Bases
- done: ring2nest(nest2ring(i)) = i for all hierarchical Healpix_Bases
- done: downgrade(upgrade(map)) = map for all maps
- done: map and downgraded map should have same average
- done: alm2map(map2alm(map)) approx map (same for pol)
- partly done: neighbor tests
- powspec -> alm -> powspec (should produce similar powspecs, also for pol)
- done: two swap_schemes() should produce original map
- done: query_disc tests (dot products etc.)
- a_lms: test Set(), Scale(), Add(), alm(l,m) = alm.mstart(m)[l], etc.

*/

#include <iostream>
#include "healpix_base.h"
#include "healpix_map.h"
#include "arr.h"
#include "planck_rng.h"
#include "lsconstants.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "geom_utils.h"
#include "walltimer.h"
#include "announce.h"

using namespace std;

const int nsamples = 1000000;

planck_rng rng;

namespace {

void random_dir (pointing &ptg)
  {
  ptg.theta = acos(rng.rand_uni()*2-1);
  ptg.phi = rng.rand_uni()*twopi;
  }
void random_zphi (double &z, double &phi)
  {
  z = rng.rand_uni()*2-1;
  phi = rng.rand_uni()*twopi;
  }

template<typename I> string bname()
  { return string("(basetype: ")+type2typename<I>()+")"; }

template<typename I> void check_ringnestring()
  {
  cout << "testing ring2nest(nest2ring(m))==m " << bname<I>() << endl;
  for (int order=0; order<=T_Healpix_Base<I>::order_max; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.ring2nest(base.nest2ring(pix))!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      }
    }
  }

template<typename I> void check_nestpeanonest()
  {
  cout << "testing peano2nest(nest2peano(m))==m " << bname<I>() << endl;
  for (int order=0; order<=T_Healpix_Base<I>::order_max; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.peano2nest(base.nest2peano(pix))!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      }
    }
  }

template<typename I> void check_pixzphipix()
  {
  cout << "testing zphi2pix(pix2zphi(m))==m " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,RING), base2 (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi;
      I pix = I(rng.rand_uni()*base1.Npix());
      base1.pix2zphi(pix,z,phi);
      if (base1.zphi2pix(z,phi)!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      base2.pix2zphi(pix,z,phi);
      if (base2.zphi2pix(z,phi)!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi;
      I pix = I(rng.rand_uni()*base.Npix());
      base.pix2zphi(pix,z,phi);
      if (base.zphi2pix(z,phi)!=pix)
        cout << "  PROBLEM: nside = " << nside << ", pixel = " << pix << endl;
      }
    }
  }

template<typename I> void check_zphipixzphi()
  {
  cout << "testing pix2zphi(zphi2pix(ptg)) approx zphi " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,NEST), base2 (order,RING);
    double mincos = min (cos(base1.max_pixrad()),0.999999999999999);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi,z2,phi2;
      random_zphi (z,phi);
      base1.pix2zphi(base1.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        cout << "  PROBLEM: order = " << order
             << ", zphi = " << z << ", " << phi << endl;
      base2.pix2zphi(base2.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        cout << "  PROBLEM: order = " << order
             << ", zphi = " << z << ", " << phi << endl;
      }
    }
  for (int nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    double mincos = min (cos(base.max_pixrad()),0.999999999999999);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi,z2,phi2;
      random_zphi (z,phi);
      base.pix2zphi(base.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        cout << "  PROBLEM: nside = " << nside
             << ", zphi = " << z << ", " << phi << endl;
      }
    }
  }

template<typename I> void check_pixangpix()
  {
  cout << "testing ang2pix(pix2ang(m))==m " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,RING), base2 (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base1.Npix());
      if (base1.ang2pix(base1.pix2ang(pix))!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      if (base2.ang2pix(base2.pix2ang(pix))!=pix)
        cout << "  PROBLEM: order = " << order << ", pixel = " << pix << endl;
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.ang2pix(base.pix2ang(pix))!=pix)
        cout << "  PROBLEM: nside = " << nside << ", pixel = " << pix << endl;
      }
    }
  }

template<typename I> void check_neighbors()
  {
  cout << "testing neighbor function " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST), base2(order,RING);
    double maxang = 2.01*base.max_pixrad();
    for (int m=0; m<nsamples/10; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      fix_arr<I,8> nb,nb2;
      vec3 pixpt = base.pix2vec(pix);
      base.neighbors(pix,nb);
      base2.neighbors(base.nest2ring(pix),nb2);
      for (int n=0; n<8; ++n)
        if (nb[n]<0)
          planck_assert(nb2[n]<0,"neighbor inconsistency");
        else
          planck_assert(base.nest2ring(nb[n])==nb2[n],"neighbor inconsistency");
      sort(&nb[0],&nb[0]+8);
      int check=0;
      for (int n=0; n<8; ++n)
        {
        if (nb[n]<0)
          ++check;
        else
          {
          if (v_angle(base.pix2vec(nb[n]),pixpt)>maxang)
            cout << " PROBLEM: order = " << order << ", pix = " << pix << endl;
          if ((n>0) && (nb[n]==nb[n-1]))
            cout << " PROBLEM: order = " << order << ", pix = " << pix << endl;
          }
        }
      planck_assert((check<=1)||((order==0)&&(check<=2)),"too few neighbors");
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    double maxang = 2.01*base.max_pixrad();
    for (int m=0; m<nsamples/10; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      fix_arr<I,8> nb;
      vec3 pixpt = base.pix2vec(pix);
      base.neighbors(pix,nb);
      for (int n=0; n<8; ++n)
        if ((nb[n]>=0) && (v_angle(base.pix2vec(nb[n]),pixpt)>maxang))
          cout << "  PROBLEM: nside = " << nside << ", pix = " << pix << endl;
      }
    }
  }

void check_swap_scheme()
  {
  cout << "testing whether double swap_scheme() returns the original map"
       << endl << "(for orders 0 to 10)." << endl;
  for (int order=0; order<=10; ++order)
    {
    Healpix_Map<uint8> map(order,NEST);
    for (int m=0; m<map.Npix(); ++m) map[m]=uint8(m&0xFF);
    map.swap_scheme();
    map.swap_scheme();
    for (int m=0; m<map.Npix(); ++m)
      if (map[m]!=(m&0xFF))
        cout << "  PROBLEM: order = " << order << ", pix = " << m << endl;
    }
  }

void check_query_disc_strict (Healpix_Ordering_Scheme scheme)
  {
  cout << "testing whether all pixels found by query_disc() really" << endl
       << "lie inside the disk (and vice versa)" << endl;
  cout << "Ordering scheme: " << (scheme==RING ? "RING" : "NEST") << endl;
  for (int order=0; order<=5; ++order)
    {
    Healpix_Map<bool> map (order,scheme);
    map.fill(false);
    Healpix_Map<vec3> vmap(order,scheme);
    for (int m=0; m<vmap.Npix(); ++m)
      vmap[m]=vmap.pix2vec(m);
    rangeset<int> pixset;
    for (int m=0; m<100000; ++m)
      {
      pointing ptg;
      random_dir (ptg);
      double rad = pi/1 * rng.rand_uni();
      map.query_disc(ptg,rad,pixset);
      vec3 vptg=ptg;
      double cosrad=cos(rad);
      for (tsize j=0; j<pixset.size(); ++j)
        for (int i=pixset.ivbegin(j); i<pixset.ivend(j); ++i)
          map[i] = true;
      for (int i=0; i<map.Npix(); ++i)
        {
        bool inside = dotprod(vmap[i],vptg)>cosrad;
        if (inside^map[i])
          cout << "  PROBLEM: order = " << order << ", ptg = " << ptg << endl;
        }
      for (tsize j=0; j<pixset.size(); ++j)
        for (int i=pixset.ivbegin(j); i<pixset.ivend(j); ++i)
          map[i] = false;
      }
    }
  }

template<typename I>void check_query_disc()
  {
  cout << "checking query_disc() " << bname<I>() << endl;
  int omax=min(20,T_Healpix_Base<I>::order_max);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> rbase (order,RING), nbase (order,NEST);
    rangeset<I> pixset;
    int niter=max(1,min(1000,100000>>order));
    for (int m=0; m<niter; ++m)
      {
      pointing ptg;
      random_dir (ptg);
      double rad = pi/1 * rng.rand_uni();
      rbase.query_disc(ptg,rad,pixset);
      rangeset<I> pslast=pixset;
      for (tsize fct=5; fct>0; --fct)
        {
        rangeset<I> psi;
        rbase.query_disc_inclusive(ptg,rad,psi,fct);
        if (!psi.contains(pslast))
          cout << "  PROBLEM: pixel sets inconsistent" << endl;
        swap(pslast,psi);
        }
      I nval = pixset.nval();
      nbase.query_disc(ptg,rad,pixset);
      pslast=pixset;
      for (tsize fct=8; fct>0; fct>>=1)
        {
        rangeset<I> psi;
        nbase.query_disc_inclusive(ptg,rad,psi,fct);
        if (!psi.contains(pslast))
          cout << "  PROBLEM: pixel sets inconsistent" << endl;
        swap(pslast,psi);
        }
      if (nval!=pixset.nval())
        cout << "  PROBLEM: number of pixels different: "
             << nval << " vs. " << pixset.nval() << endl;
      }
    }
  }
template<typename I>void check_query_polygon()
  {
  cout << "checking query_polygon() " << bname<I>() << endl;
  int omax=min(20,T_Healpix_Base<I>::order_max);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> rbase (order,RING), nbase (order,NEST);
    rangeset<I> pixset;
    int niter=max(1,min(1000,100000>>order));
    for (int m=0; m<niter; ++m)
      {
      vector<pointing> corner(3);
      random_dir(corner[0]); random_dir(corner[1]); random_dir(corner[2]);
      rbase.query_polygon(corner,pixset);
      I nval = pixset.nval();
      nbase.query_polygon(corner,pixset);
      if (nval!=pixset.nval())
        cout << "  PROBLEM: number of pixels different: "
             << nval << " vs. " << pixset.nval() << endl;
      rbase.query_polygon_inclusive(corner,pixset,4);
      I nv1=pixset.nval();
      nbase.query_polygon_inclusive(corner,pixset,4);
      I nv2=pixset.nval();
      if (nv1<nv2)
        cout << "  PROBLEM: inclusive(RING)<inclusive(NEST): "
             << nv1 << " vs. " << nv2 << endl;
      if (nv2<nval)
        cout << "  PROBLEM: inclusive(NEST)<non-inclusive: "
             << nv2 << " vs. " << nval << endl;
      }
    }
  }

void helper_oop (int order)
  {
  Healpix_Map<double> map (order,RING), map2 (order,NEST), map3 (order,RING);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      cout << "PROBLEM: order = " << order << endl;
  }
void helper_udgrade (int order, Healpix_Ordering_Scheme s1,
  Healpix_Ordering_Scheme s2)
  {
  Healpix_Map<double> map (order,s1), map2 (order+2,s2), map3 (order, s1);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      cout << "PROBLEM: order = " << order << endl;
  }
void helper_udgrade2 (int nside)
  {
  Healpix_Map<double> map (nside,RING,SET_NSIDE), map2 (nside*3,RING,SET_NSIDE),
    map3 (nside, RING,SET_NSIDE);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      cout << "PROBLEM: nside = " << nside << endl;
  }

void check_import()
  {
  cout << "testing out-of-place swapping" << endl;
  for (int order=0; order<=7; ++order)
    helper_oop(order);
  cout << "testing downgrade(upgrade(map)) == map" << endl;
  for (int order=0; order<=7; ++order)
    {
    helper_udgrade(order,RING,RING);
    helper_udgrade(order,RING,NEST);
    helper_udgrade(order,NEST,NEST);
    helper_udgrade(order,NEST,RING);
    }
  for (int nside=3; nside<500; nside+=nside/2+1)
    helper_udgrade2(nside);
  }

void check_average()
  {
  cout << "testing whether average(map) == average(downgraded map)" << endl;
  for (int order=1; order<=10; ++order)
    {
    Healpix_Map<double> map (order,RING), map2(1,RING);
    for (int m=0; m<map.Npix(); ++m)
      map[m] = rng.rand_uni()+0.01;
    map2.Import(map);
    double avg=map.average(), avg2=map2.average();
    if (!approx(avg,avg2,1e-10))
      cout << "PROBLEM: order = " << order << " " << avg/avg2-1 << endl;
    }
  for (int nside=3; nside<1000; nside += nside/2+1)
    {
    Healpix_Map<double> map (nside,RING,SET_NSIDE), map2(1,RING,SET_NSIDE);
    for (int m=0; m<map.Npix(); ++m)
      map[m] = rng.rand_uni()+0.01;
    map2.Import(map);
    double avg=map.average(), avg2=map2.average();
    if (!approx(avg,avg2,1e-10))
      cout << "PROBLEM: nside = " << nside << " " << avg/avg2-1 << endl;
    }
  }

void random_alm (Alm<xcomplex<double> >&almT, Alm<xcomplex<double> >&almG,
  Alm<xcomplex<double> >&almC, int lmax, int mmax)
  {
  almT.Set(lmax,mmax); almG.Set(lmax,mmax); almC.Set(lmax,mmax);

  for (int l=0; l<=lmax; ++l)
    {
    almT(l,0).re=rng.rand_gauss(); almT(l,0).im=0.;
    almG(l,0).re=rng.rand_gauss(); almG(l,0).im=0.;
    almC(l,0).re=rng.rand_gauss(); almC(l,0).im=0.;
    }

  for (int m=1; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      {
      almT(l,m).re=rng.rand_gauss(); almT(l,m).im=rng.rand_gauss();
      almG(l,m).re=rng.rand_gauss(); almG(l,m).im=rng.rand_gauss();
      almC(l,m).re=rng.rand_gauss(); almC(l,m).im=rng.rand_gauss();
      }
  almG(0,0)=almC(0,0)=almG(1,0)=almC(1,0)=almG(1,1)=almC(1,1)=0;
  }

void random_alm (Alm<xcomplex<double> >&alm, int lmax, int mmax)
  {
  alm.Set(lmax,mmax);

  for (int l=0; l<=lmax; ++l)
    { alm(l,0).re=rng.rand_gauss(); alm(l,0).im=0.; }

  for (int m=1; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      { alm(l,m).re=rng.rand_gauss(); alm(l,m).im=rng.rand_gauss(); }
  }

void check_alm (const Alm<xcomplex<double> >&oalm,
  const Alm<xcomplex<double> >&alm, double epsilon)
  {
  for (int m=0; m<=alm.Mmax(); ++m)
    for (int l=m; l<=alm.Lmax(); ++l)
      {
      if (!abs_approx(oalm(l,m).re,alm(l,m).re,epsilon))
        cout << "Problemr " << l << " " << m << endl;
      if (!abs_approx(oalm(l,m).im,alm(l,m).im,epsilon))
        cout << "Problemi " << l << " " << m <<  endl;
      }
  }

void check_alm2map2alm (int lmax, int mmax, int nside)
  {
  cout << "testing whether a_lm->map->a_lm returns original a_lm" << endl;
  cout << "lmax=" << lmax <<", mmax=" << mmax << ", nside=" << nside << endl;
  const double epsilon = 1e-8;
  Alm<xcomplex<double> > oalmT(lmax,mmax),almT(lmax,mmax),
    oalmG(lmax,mmax),almG(lmax,mmax),oalmC(lmax,mmax),almC(lmax,mmax);
  Healpix_Map<double> mapT(nside,RING,SET_NSIDE), mapQ(nside,RING,SET_NSIDE),
    mapU(nside,RING,SET_NSIDE);

  random_alm(oalmT,oalmG,oalmC,lmax,mmax);
  alm2map(oalmT,mapT);
  map2alm_iter2(mapT,almT,1e-12,1e-12);
  check_alm (oalmT, almT, epsilon);

  alm2map_spin(oalmG,oalmC,mapQ,mapU,1);
  map2alm_spin_iter2(mapQ,mapU,almG,almC,1,1e-12,1e-12);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);

  alm2map_pol(oalmT,oalmG,oalmC,mapT,mapQ,mapU);
  map2alm_pol_iter2(mapT,mapQ,mapU,almT,almG,almC,1e-12,1e-12);
  check_alm (oalmT, almT, epsilon);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);
  }

void check_smooth_alm ()
  {
  cout << "testing whether unsmooth(smooth(a_lm)) returns a_lm" << endl;
  const double epsilon = 1e-14;
  const double fwhm = 100.*arcmin2rad;
  const int lmax=300, mmax=300;
  Alm<xcomplex<double> > oalmT(lmax,mmax),almT(lmax,mmax),
    oalmG(lmax,mmax),almG(lmax,mmax),oalmC(lmax,mmax),almC(lmax,mmax);

  random_alm(oalmT,oalmG,oalmC,lmax,mmax);

  almT=oalmT; almG=oalmG; almC=oalmC;
  smoothWithGauss (almT, fwhm);
  smoothWithGauss (almT, -fwhm);
  check_alm (oalmT, almT, epsilon);
  almT=oalmT;
  smoothWithGauss (almT, almG, almC, fwhm);
  smoothWithGauss (almT, almG, almC, -fwhm);
  check_alm (oalmT, almT, epsilon);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);
  }

void check_rot_alm ()
  {
  cout << "testing whether rot^-1(rot(a_lm)) returns a_lm" << endl;
  const double epsilon = 2e-13;
  const int lmax=300;
  Alm<xcomplex<double> > oalm(lmax,lmax),alm(lmax,lmax);

  random_alm(oalm,lmax,lmax);

  alm=oalm;
  rotate_alm (alm,3,4,5);
  rotate_alm (alm,-5,-4,-3);
  check_alm (oalm, alm, epsilon);
  }

void check_isqrt()
  {
  cout << "testing whether isqrt() works reliably" << endl;
  uint64 val=uint64(0xF234)<<16, valsq=val*val;
  if (isqrt(valsq)!=val) cout << "PROBLEM1" << endl;
  if (isqrt(valsq-1)!=val-1) cout << "PROBLEM2" << endl;
  }

void check_pix2ang_acc()
  {
  cout << "testing accuracy of pix2ang at the poles" << endl;
  for (int m=0; m<=29;++m)
    {
    Healpix_Base2 base(m,RING);
    if (base.pix2ang(1).theta==0.)
      cout << "PROBLEM: order " << m << endl;
    }
  }

const int nsteps=1000000;

template<typename I>void perf_neighbors(const string &name,
  Healpix_Ordering_Scheme scheme)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    fix_arr<I,8> nres;
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      base.neighbors(pix,nres);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2ang(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      pointing p(base.pix2ang(pix));
      dummy+=p.theta+p.phi;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2vec(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      vec3 v(base.pix2vec(pix));
      dummy+=v.x+v.y+v.z;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2zphi(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    double z,phi;
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      base.pix2zphi(pix,z,phi);
      dummy+=z+phi;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_zphi2pix(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  double dz=2./sqrt(nsteps);
  double dph=twopi/sqrt(nsteps);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    for (double z=-1; z<1; z+=dz)
      for (double phi=0; phi<twopi; phi+=dph)
        {
        dummy+=base.zphi2pix(z,phi);
        ++cnt;
        }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_ang2pix(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  double dth=pi/sqrt(nsteps);
  double dph=twopi/sqrt(nsteps);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    for (double theta=0; theta<pi; theta+=dth)
      for (double phi=0; phi<twopi; phi+=dph)
        {
        dummy+=base.ang2pix(pointing(theta+1e-15*phi,phi));
        ++cnt;
        }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_ring2nest(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.ring2nest(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_nest2ring(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.nest2ring(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_peano2nest(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.peano2nest(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_nest2peano(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.nest2peano(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_query_disc(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    rangeset<I> pix;
    base.query_disc(vec3(1,0,0),halfpi/9,pix);
    dummy+=pix.size();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }
template<typename I>void perf_query_triangle(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  vector<pointing> corner;
  corner.push_back(vec3(1,0.01,0.01));
  corner.push_back(vec3(0.01,1,0.01));
  corner.push_back(vec3(0.01,0.01,1));
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    rangeset<I> pix;
    base.query_polygon(corner,pix);
    dummy+=pix.size();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }
template<typename I>void perf_query_polygon(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  vector<pointing> corner;
  corner.push_back(vec3(1,0.01,0.01));
  corner.push_back(vec3(1,1,-0.3));
  corner.push_back(vec3(0.01,1,0.01));
  corner.push_back(vec3(0.01,0.01,1));
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    rangeset<I> pix;
    base.query_polygon(corner,pix);
    dummy+=pix.size();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }

void perftest()
  {
  double dummy=0;
  cout << "Measuring performance of Healpix_Base methods." << endl;
  perf_pix2zphi<int>   ("pix2zphi (RING):int  ",RING,dummy);
  perf_pix2zphi<int>   ("pix2zphi (NEST):int  ",NEST,dummy);
  perf_pix2zphi<int64> ("pix2zphi (RING):int64",RING,dummy);
  perf_pix2zphi<int64> ("pix2zphi (NEST):int64",NEST,dummy);
  perf_zphi2pix<int>   ("zphi2pix (RING):int  ",RING,dummy);
  perf_zphi2pix<int>   ("zphi2pix (NEST):int  ",NEST,dummy);
  perf_zphi2pix<int64> ("zphi2pix (RING):int64",RING,dummy);
  perf_zphi2pix<int64> ("zphi2pix (NEST):int64",NEST,dummy);
  perf_pix2ang<int>    ("pix2ang  (RING):int  ",RING,dummy);
  perf_pix2ang<int>    ("pix2ang  (NEST):int  ",NEST,dummy);
  perf_pix2ang<int64>  ("pix2ang  (RING):int64",RING,dummy);
  perf_pix2ang<int64>  ("pix2ang  (NEST):int64",NEST,dummy);
  perf_ang2pix<int>    ("ang2pix  (RING):int  ",RING,dummy);
  perf_ang2pix<int>    ("ang2pix  (NEST):int  ",NEST,dummy);
  perf_ang2pix<int64>  ("ang2pix  (RING):int64",RING,dummy);
  perf_ang2pix<int64>  ("ang2pix  (NEST):int64",NEST,dummy);
  perf_pix2vec<int>    ("pix2vec  (RING):int  ",RING,dummy);
  perf_pix2vec<int>    ("pix2vec  (NEST):int  ",NEST,dummy);
  perf_pix2vec<int64>  ("pix2vec  (RING):int64",RING,dummy);
  perf_pix2vec<int64>  ("pix2vec  (NEST):int64",NEST,dummy);
  perf_neighbors<int>  ("neighbors(NEST):int  ",NEST);
  perf_neighbors<int>  ("neighbors(RING):int  ",RING);
  perf_neighbors<int64>("neighbors(NEST):int64",NEST);
  perf_neighbors<int64>("neighbors(RING):int64",RING);
  perf_ring2nest<int>  ("ring2nest      :int  ",dummy);
  perf_ring2nest<int64>("ring2nest      :int64",dummy);
  perf_nest2ring<int>  ("nest2ring      :int  ",dummy);
  perf_nest2ring<int64>("nest2ring      :int64",dummy);
  perf_peano2nest<int>  ("peano2nest     :int  ",dummy);
  perf_peano2nest<int64>("peano2nest     :int64",dummy);
  perf_nest2peano<int>  ("nest2peano     :int  ",dummy);
  perf_nest2peano<int64>("nest2peano     :int64",dummy);
  perf_query_disc<int>      ("query_disc    (RING):int  ",RING,dummy);
  perf_query_disc<int>      ("query_disc    (NEST):int  ",NEST,dummy);
  perf_query_disc<int64>    ("query_disc    (RING):int64",RING,dummy);
  perf_query_disc<int64>    ("query_disc    (NEST):int64",NEST,dummy);
  perf_query_triangle<int>  ("query_triangle(RING):int  ",RING,dummy);
  perf_query_triangle<int>  ("query_triangle(NEST):int  ",NEST,dummy);
  perf_query_triangle<int64>("query_triangle(RING):int64",RING,dummy);
  perf_query_triangle<int64>("query_triangle(NEST):int64",NEST,dummy);
  perf_query_polygon<int>   ("query_polygon (RING):int  ",RING,dummy);
  perf_query_polygon<int>   ("query_polygon (NEST):int  ",NEST,dummy);
  perf_query_polygon<int64> ("query_polygon (RING):int64",RING,dummy);
  perf_query_polygon<int64> ("query_polygon (NEST):int64",NEST,dummy);

  if (dummy<0) cout << dummy << endl;
  }
} // unnamed namespace

int main(int argc, const char **argv)
  {
  module_startup ("hpxtest",argc,argv,1,"");
  perftest();
  check_isqrt();
  check_pix2ang_acc();
  check_smooth_alm();
  check_rot_alm();
  check_alm2map2alm(620,620,256);
  check_alm2map2alm(620,2,256);
  check_average();
  check_import();
  check_ringnestring<int>();
  check_ringnestring<int64>();
  check_nestpeanonest<int>();
  check_nestpeanonest<int64>();
  check_pixzphipix<int>();
  check_pixzphipix<int64>();
  check_zphipixzphi<int>();
  check_zphipixzphi<int64>();
  check_pixangpix<int>();
  check_pixangpix<int64>();
  check_neighbors<int>();
  check_neighbors<int64>();
  check_swap_scheme();
  check_query_disc_strict(RING);
  check_query_disc_strict(NEST);
  check_query_disc<int>();
  check_query_disc<int64>();
  check_query_polygon<int>();
  check_query_polygon<int64>();
  }

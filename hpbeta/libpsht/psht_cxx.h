/*
 *  This file is part of libpsht.
 *
 *  libpsht is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libpsht is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libpsht; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libpsht is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file psht_cxx.h
 *  Spherical transform library
 *
 *  Copyright (C) 2008, 2009 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_PSHT_CXX_H
#define PLANCK_PSHT_CXX_H

#include "psht.h"
#include "psht_geomhelpers.h"
#include "psht_almhelpers.h"
#include "xcomplex.h"

class psht_base
  {
  protected:
    psht_geom_info *ginfo;
    psht_alm_info *ainfo;

  public:
    psht_base()
      : ginfo(0), ainfo(0) {}
    ~psht_base()
      {
      psht_destroy_geom_info(ginfo);
      psht_destroy_alm_info(ainfo);
      }

    void set_general_geometry (int nrings, const int *nph, const ptrdiff_t *ofs,
      const int *stride, const double *phi0, const double *theta,
      const double *weight)
      {
      psht_make_geom_info (nrings, nph, ofs, stride, phi0, theta, weight,
        &ginfo);
      }

    void set_ECP_geometry (int nrings, int nphi)
      { psht_make_ecp_geom_info (nrings, nphi, 0., 1, &ginfo); }

    void set_Healpix_geometry (int nside)
      { psht_make_healpix_geom_info (nside, 1, &ginfo); }

    void set_weighted_Healpix_geometry (int nside, const double *weight)
      { psht_make_weighted_healpix_geom_info (nside, 1, weight, &ginfo); }

    void set_triangular_alm_info (int lmax, int mmax)
      { psht_make_triangular_alm_info (lmax, mmax, 1, &ainfo); }
  };

template<typename T> class psht_joblist
  {};

template<> class psht_joblist<float>: public psht_base
  {
  private:
    pshts_joblist *jobs;

    static pshts_cmplx *conv (xcomplex<float> *ptr)
      { return reinterpret_cast<pshts_cmplx *>(ptr); }
    static const pshts_cmplx *conv (const xcomplex<float> *ptr)
      { return reinterpret_cast<const pshts_cmplx *>(ptr); }

  public:
    psht_joblist()
      { pshts_make_joblist (&jobs); }
    ~psht_joblist()
      { pshts_destroy_joblist (jobs); }

    void clear_jobs()
      { pshts_clear_joblist (jobs); }

    void add_alm2map (const xcomplex<float> *alm, float *map, bool add)
      { pshts_add_job_alm2map (jobs, conv(alm), map, add); }
    void add_map2alm (const float *map, xcomplex<float> *alm, bool add)
      { pshts_add_job_map2alm (jobs, map, conv(alm), add); }
    void add_alm2map_pol (const xcomplex<float> *almT,
      const xcomplex<float> *almG, const xcomplex<float> *almC,
      float *mapT, float *mapQ, float *mapU, bool add)
      {
      pshts_add_job_alm2map_pol (jobs, conv(almT), conv(almG), conv(almC),
        mapT, mapQ, mapU, add);
      }
    void add_map2alm_pol (const float *mapT, const float *mapQ,
      const float *mapU, xcomplex<float> *almT, xcomplex<float> *almG,
      xcomplex<float> *almC, bool add)
      {
      pshts_add_job_map2alm_pol (jobs, mapT, mapQ, mapU, conv(almT),
        conv(almG), conv(almC), add);
      }
    void add_alm2map_spin (const xcomplex<float> *alm1,
      const xcomplex<float> *alm2, float *map1, float *map2, int spin, bool add)
      {
      pshts_add_job_alm2map_spin (jobs, conv(alm1), conv(alm2), map1, map2,
        spin, add);
      }
    void add_map2alm_spin (const float *map1, const float *map2,
      xcomplex<float> *alm1, xcomplex<float> *alm2, int spin, bool add)
      {
      pshts_add_job_map2alm_spin (jobs, map1, map2, conv(alm1), conv(alm2),
        spin, add);
      }
    void add_alm2map_der1 (const xcomplex<float> *alm, float *mapdth,
      float *mapdph, bool add)
      {
      pshts_add_job_alm2map_deriv1 (jobs, conv(alm), mapdth, mapdph, add);
      }

    void execute()
      { pshts_execute_jobs (jobs,ginfo,ainfo); }
  };

template<> class psht_joblist<double>: public psht_base
  {
  private:
    pshtd_joblist *jobs;

    static pshtd_cmplx *conv (xcomplex<double> *ptr)
      { return reinterpret_cast<pshtd_cmplx *>(ptr); }
    static const pshtd_cmplx *conv (const xcomplex<double> *ptr)
      { return reinterpret_cast<const pshtd_cmplx *>(ptr); }

  public:
    psht_joblist()
      { pshtd_make_joblist (&jobs); }
    ~psht_joblist()
      { pshtd_destroy_joblist (jobs); }

    void clear_jobs()
      { pshtd_clear_joblist (jobs); }

    void add_alm2map (const xcomplex<double> *alm, double *map, bool add)
      { pshtd_add_job_alm2map (jobs, conv(alm), map, add); }
    void add_map2alm (const double *map, xcomplex<double> *alm, bool add)
      { pshtd_add_job_map2alm (jobs, map, conv(alm), add); }
    void add_alm2map_pol (const xcomplex<double> *almT,
      const xcomplex<double> *almG, const xcomplex<double> *almC,
      double *mapT, double *mapQ, double *mapU, bool add)
      {
      pshtd_add_job_alm2map_pol (jobs, conv(almT), conv(almG), conv(almC),
        mapT, mapQ, mapU, add);
      }
    void add_map2alm_pol (const double *mapT, const double *mapQ,
      const double *mapU, xcomplex<double> *almT, xcomplex<double> *almG,
      xcomplex<double> *almC, bool add)
      {
      pshtd_add_job_map2alm_pol (jobs, mapT, mapQ, mapU, conv(almT),
        conv(almG), conv(almC), add);
      }
    void add_alm2map_spin (const xcomplex<double> *alm1,
      const xcomplex<double> *alm2, double *map1, double *map2, int spin,
      bool add)
      {
      pshtd_add_job_alm2map_spin (jobs, conv(alm1), conv(alm2), map1, map2,
        spin, add);
      }
    void add_map2alm_spin (const double *map1, const double *map2,
      xcomplex<double> *alm1, xcomplex<double> *alm2, int spin, bool add)
      {
      pshtd_add_job_map2alm_spin (jobs, map1, map2, conv(alm1), conv(alm2),
        spin, add);
      }
    void add_alm2map_der1 (const xcomplex<double> *alm, double *mapdth,
      double *mapdph, bool add)
      {
      pshtd_add_job_alm2map_deriv1 (jobs, conv(alm), mapdth, mapdph, add);
      }

    void execute()
      { pshtd_execute_jobs (jobs,ginfo,ainfo); }
  };

#endif

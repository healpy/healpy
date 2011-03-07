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

/*! \file psht.h
 *  Interface for the spherical transform library.
 *
 *  Copyright (C) 2006-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_PSHT_H
#define PLANCK_PSHT_H

#include <stddef.h>
#include "sse_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! PSHT type for storing single-precision complex values */
typedef struct
  {
  float re,im;
  } pshts_cmplx;

/*! PSHT type for storing double-precision complex values */
typedef struct
  {
  double re,im;
  } pshtd_cmplx;

extern const pshts_cmplx pshts_cmplx_null;
extern const pshtd_cmplx pshtd_cmplx_null;

/*! Helper type containing information about a single ring.
    \note No user serviceable parts inside! */
typedef struct
  {
  double theta, phi0, weight, cth, sth;
  ptrdiff_t ofs;
  int nph, stride;
  } psht_ringinfo;

/*! Helper type containing information about a pair of rings with colatitudes
    symmetric around the equator.
    \note No user serviceable parts inside! */
typedef struct
  {
  psht_ringinfo r1,r2;
  } psht_ringpair;

/*! Enumeration of PSHT job types.
    \note No user serviceable parts inside! */
typedef enum { MAP2ALM, ALM2MAP, ALM2MAP_DERIV1 } psht_jobtype;

/*! Type holding all required information about a map geometry.
    \note No user serviceable parts inside! */
typedef struct
  {
  psht_ringpair *pair;
  int npairs;
  } psht_geom_info;

/*! Type holding all required information about a double precision SHT.
    \note No user serviceable parts inside! */
typedef struct
  {
  psht_jobtype type;
  int spin;
  int add_output;
  int nmaps, nalm;
  double *map[3];
  pshtd_cmplx *alm[3];
  pshtd_cmplx *phas1[3], *phas2[3];
  double *norm_l;
  union {
#ifdef PLANCK_HAVE_SSE2
    v2df *v[3];
    v2df2 *v2[3];
#else
    pshtd_cmplx *c[3];
#endif
    } alm_tmp;
  } pshtd_job;

enum { psht_maxjobs=10 };

/*! Type holding a list of simultaneous double precision SHT jobs.
    \note No user serviceable parts inside! */
typedef struct
  {
  pshtd_job job[psht_maxjobs];
  int njobs;
  } pshtd_joblist;

/*! Type holding all required information about a single precision SHT.
    \note No user serviceable parts inside! */
typedef struct
  {
  psht_jobtype type;
  int spin;
  int add_output;
  int nmaps, nalm;
  float *map[3];
  pshts_cmplx *alm[3];
  pshtd_cmplx *phas1[3], *phas2[3];
  double *norm_l;
  union {
#ifdef PLANCK_HAVE_SSE2
    v2df *v[3];
    v2df2 *v2[3];
#else
    pshtd_cmplx *c[3];
#endif
    } alm_tmp;
  } pshts_job;

/*! Type holding a list of simultaneous single precision SHT jobs.
    \note No user serviceable parts inside! */
typedef struct
  {
  pshts_job job[psht_maxjobs];
  int njobs;
  } pshts_joblist;

/*! \defgroup almgroup Helpers for calculation of a_lm indices */
/*! \{ */

/*! Helper type for index calculation in a_lm arrays. */
typedef struct
  {
  /*! Maximum \a l index of the array */
  int lmax;
  /*! Number of different \a m values in this object */
  int nm;
  /*! Array with \a nm entries containing the individual m values */
  int *mval;
  /*! Array with \a nm entries containing the (hypothetical) indices of
      the coefficients with quantum numbers 0,\a mval[i] */
  ptrdiff_t *mvstart;
  /*! Stride between a_lm and a_(l+1),m */
  ptrdiff_t stride;
  } psht_alm_info;

/*! Creates an Alm data structure information from the following parameters:
    \param lmax maximum \a l quantum number (>=0)
    \param mmax maximum \a m quantum number (0<= \a mmax <= \a lmax)
    \param stride the stride between consecutive a_lm entries
    \param mstart the index of the (hypothetical) coefficient with the
      quantum numbers 0,\a m. Must have \a mmax+1 entries.
    \param alm_info will hold a pointer to the newly created data structure
 */
void psht_make_alm_info (int lmax, int mmax, int stride,
  const ptrdiff_t *mstart, psht_alm_info **alm_info);
/*! Creates an Alm data structure information from the following parameters:
    \param lmax maximum \a l quantum number (>=0)
    \param nm number of different \a m (<=\a lmax+1)
    \param stride the stride between consecutive a_lm entries
    \param mval array with \a nm entries containing the individual m values
    \param mvstart array with \a nm entries containing the (hypothetical)
      indices of the coefficients with the quantum numbers 0,\a mval[i]
    \param alm_info will hold a pointer to the newly created data structure
 */
void psht_make_general_alm_info (int lmax, int nm, int stride, const int *mval,
  const ptrdiff_t *mvstart, psht_alm_info **alm_info);
/*! Returns the index of the coefficient with quantum numbers \a l,
    \a mval[mi]. */
ptrdiff_t psht_alm_index (const psht_alm_info *self, int l, int mi);
/*! Deallocates the a_lm info object. */
void psht_destroy_alm_info (psht_alm_info *info);

/* \} */

/*! \defgroup geominfogroup Functions for dealing with geometry information */
/*! \{ */

/*! Creates a geometry information from a set of ring descriptions.
    All arrays passed to this function must have \a nrings elements.
    \param nrings the number of rings in the map
    \param nph the number of pixels in each ring
    \param ofs the index of the first pixel in each ring in the map array
    \param stride the stride between consecutive pixels
    \param phi0 the azimuth (in radians) of the first pixel in each ring
    \param theta the colatitude (in radians) of each ring
    \param weight the pixel weight to be used for the ring
    \param geom_info will hold a pointer to the newly created data structure
 */
void psht_make_geom_info (int nrings, const int *nph, const ptrdiff_t *ofs,
  const int *stride, const double *phi0, const double *theta,
  const double *weight, psht_geom_info **geom_info);

/*! Deallocates the geometry information in \a info. */
void psht_destroy_geom_info (psht_geom_info *info);

/* \} */

/*! \defgroup sjoblistgroup Functions for dealing with single precision job lists
\note All pointers to maps or a_lm that are passed to the job-adding functions
must not be de-allocated until after the last call of execute_jobs() for
a particular job list! This is because PSHT does not copy the input data,
but only stores the pointers to the supplied maps and a_lm.
 */
/*! \{ */

/*! Creates a new joblist object. */
void pshts_make_joblist (pshts_joblist **joblist);
/*! Removes all jobs in \a joblist. */
void pshts_clear_joblist (pshts_joblist *joblist);
/*! Deallocates the given joblist object. */
void pshts_destroy_joblist (pshts_joblist *joblist);

/*! Adds a new scalar alm2map job to \a joblist, which reads data from \a alm
    and writes data to \a map. If \a add_output is 0, \a map will be
    overwritten, else the result will be added to \a map. */
void pshts_add_job_alm2map (pshts_joblist *joblist, const pshts_cmplx *alm,
  float *map, int add_output);
/*! Adds a new scalar map2alm job to \a joblist, which reads data from \a map
    and writes data to \a alm. If \a add_output is 0, \a alm will be
    overwritten, else the result will be added to \a alm. */
void pshts_add_job_map2alm (pshts_joblist *joblist, const float *map,
  pshts_cmplx *alm, int add_output);
/*! Adds a new polarised alm2map job to \a joblist, which reads data from
    \a almT, \a almG and \a almC and writes data to \a mapT, \a mapQ and
    \a mapU. If \a add_output is 0, the output maps will be
    overwritten, else the result will be added to the output maps. */
void pshts_add_job_alm2map_pol (pshts_joblist *joblist,
  const pshts_cmplx *almT, const pshts_cmplx *almG, const pshts_cmplx *almC,
  float *mapT, float *mapQ, float *mapU, int add_output);
/*! Adds a new polarised map2alm job to \a joblist, which reads data from
    \a mapT, \a mapQ and \a mapU and writes data to \a almT, \a almG and
    \a almC. If \a add_output is 0, the output a_lm will be
    overwritten, else the result will be added to the output a_lm. */
void pshts_add_job_map2alm_pol (pshts_joblist *joblist,
  const float *mapT, const float *mapQ, const float *mapU,
  pshts_cmplx *almT, pshts_cmplx *almG, pshts_cmplx *almC, int add_output);
/*! Adds a new spin alm2map job to \a joblist, which reads data from
    \a alm1 and \a alm2 and writes data to \a map1 and map2.
    \a spin must be 1, 2, or 3. If \a add_output is 0,
    the output maps will be overwritten, else the result will be
    added to the output maps. */
void pshts_add_job_alm2map_spin (pshts_joblist *joblist,
  const pshts_cmplx *alm1, const pshts_cmplx *alm2, float *map1, float *map2,
  int spin, int add_output);
/*! Adds a new spin map2alm job to \a joblist, which reads data from
    \a map1 and \a map2 and writes data to \a alm1 and \a alm2.
    \a spin must be 1, 2, or 3. If \a add_output is 0,
    the output a_lm will be overwritten, else the result will be added
    to the output a_lm. */
void pshts_add_job_map2alm_spin (pshts_joblist *joblist, const float *map1,
  const float *map2, pshts_cmplx *alm1, pshts_cmplx *alm2, int spin,
  int add_output);
/*! Adds a new job to \a joblist, which reads data from
    \a alm and writes maps of the first derivatives to \a mapdtheta and
    \a mapdphi, respectively. If \a add_output is 0,
    the output maps will be overwritten, else the result will be added
    to the output maps. */
void pshts_add_job_alm2map_deriv1 (pshts_joblist *joblist,
  const pshts_cmplx *alm, float *mapdtheta, float *mapdphi, int add_output);

/*! Executes the jobs in \a joblist, using \a geom_info as map geometry
    and \a alm_info as structure of the a_lm coefficients.
    \note The map geometry and the a_lm structure have to be supplied to this
    function only, since this information is not needed by PSHT anywhere else.
    However, it is the user's responsibility to ensure that the input arrays
    (specified by calls to the job-adding functions) are consistent with the
    specified geometry and a_lm structure, and that the output arrays are
    large enough to hold the produced results.
 */
void pshts_execute_jobs (pshts_joblist *joblist,
  const psht_geom_info *geom_info, const psht_alm_info *alm_info);

/* \} */

/*! \defgroup djoblistgroup Functions for dealing with double precision job lists
\note All pointers to maps or a_lm that are passed to the job-adding functions
must not be de-allocated until after the last call of execute_jobs() for
a particular job list! This is because PSHT does not copy the input data,
but only stores the pointers to the supplied maps and a_lm.
*/
/*! \{ */

/*! Creates a new joblist object. */
void pshtd_make_joblist (pshtd_joblist **joblist);
/*! Removes all jobs in \a joblist. */
void pshtd_clear_joblist (pshtd_joblist *joblist);
/*! Deallocates the given joblist object. */
void pshtd_destroy_joblist (pshtd_joblist *joblist);

/*! Adds a new scalar alm2map job to \a joblist, which reads data from \a alm
    and writes data to \a map. If \a add_output is 0, \a map will be
    overwritten, else the result will be added to \a map. */
void pshtd_add_job_alm2map (pshtd_joblist *joblist, const pshtd_cmplx *alm,
  double *map, int add_output);
/*! Adds a new scalar map2alm job to \a joblist, which reads data from \a map
    and writes data to \a alm. If \a add_output is 0, \a alm will be
    overwritten, else the result will be added to \a alm. */
void pshtd_add_job_map2alm (pshtd_joblist *joblist, const double *map,
  pshtd_cmplx *alm, int add_output);
/*! Adds a new polarised alm2map job to \a joblist, which reads data from
    \a almT, \a almG and \a almC and writes data to \a mapT, \a mapQ and
    \a mapU. If \a add_output is 0, the output maps will be
    overwritten, else the result will be added to the output maps. */
void pshtd_add_job_alm2map_pol (pshtd_joblist *joblist,
  const pshtd_cmplx *almT, const pshtd_cmplx *almG, const pshtd_cmplx *almC,
  double *mapT, double *mapQ, double *mapU, int add_output);
/*! Adds a new polarised map2alm job to \a joblist, which reads data from
    \a mapT, \a mapQ and \a mapU and writes data to \a almT, \a almG and
    \a almC. If \a add_output is 0, the output a_lm will be
    overwritten, else the result will be added to the output a_lm. */
void pshtd_add_job_map2alm_pol (pshtd_joblist *joblist,
  const double *mapT, const double *mapQ, const double *mapU,
  pshtd_cmplx *almT, pshtd_cmplx *almG, pshtd_cmplx *almC, int add_output);
/*! Adds a new spin alm2map job to \a joblist, which reads data from
    \a alm1 and \a alm2 and writes data to \a map1 and map2.
    \a spin must be 1, 2, or 3. If \a add_output is 0,
    the output maps will be overwritten, else the result will be
    added to the output maps. */
void pshtd_add_job_alm2map_spin (pshtd_joblist *joblist,
  const pshtd_cmplx *alm1, const pshtd_cmplx *alm2, double *map1, double *map2,
  int spin, int add_output);
/*! Adds a new spin map2alm job to \a joblist, which reads data from
    \a map1 and \a map2 and writes data to \a alm1 and \a alm2.
    \a spin must be 1, 2, or 3. If \a add_output is 0,
    the output a_lm will be overwritten, else the result will be added
    to the output a_lm. */
void pshtd_add_job_map2alm_spin (pshtd_joblist *joblist, const double *map1,
  const double *map2, pshtd_cmplx *alm1, pshtd_cmplx *alm2, int spin,
  int add_output);
/*! Adds a new job to \a joblist, which reads data from
    \a alm and writes maps of the first derivatives to \a mapdtheta and
    \a mapdphi, respectively. If \a add_output is 0,
    the output maps will be overwritten, else the result will be added
    to the output maps. */
void pshtd_add_job_alm2map_deriv1 (pshtd_joblist *joblist,
  const pshtd_cmplx *alm, double *mapdtheta, double *mapdphi, int add_output);

/*! Executes the jobs in \a joblist, using \a geom_info as map geometry
    and \a alm_info as structure of the a_lm coefficients.
    \note The map geometry and the a_lm structure have to be supplied to this
    function only, since this information is not needed by PSHT anywhere else.
    However, it is the user's responsibility to ensure that the input arrays
    (specified by calls to the job-adding functions) are consistent with the
    specified geometry and a_lm structure, and that the output arrays are
    large enough to hold the produced results.
 */
void pshtd_execute_jobs (pshtd_joblist *joblist,
  const psht_geom_info *geom_info, const psht_alm_info *alm_info);

/* \} */

#ifdef __cplusplus
}
#endif

#endif

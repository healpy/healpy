# Wrapper around the geometry methods in Healpix_base class

import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
cimport cython

cdef extern from "healpix_base.h":
    ctypedef long int64

cdef extern from "stddef.h":
    ctypedef long ptrdiff_t

ctypedef size_t tsize
ctypedef ptrdiff_t tdiff

cdef extern from "rangeset.h":
    cdef cppclass rangeset[T]:
        tsize size()
        int64 nval()
        T ivbegin(tdiff i)
        T ivend(tdiff i)

cdef extern from "arr.h":
    cdef cppclass arr[T]:
        pass
    cdef cppclass fix_arr[T, int]:
        pass

cdef extern from "vec3.h":
    cdef cppclass vec3:
        double x, y, z
        vec3()
        vec3(double xc, double yc, double zc)

cdef extern from "pointing.h":
    cdef cppclass pointing:
        pointing()
        pointing(vec3 inp)
        double theta
        double phi

cdef extern from "healpix_base.h":
    ctypedef int val_4 "4"
    ctypedef int val_8 "8"
    cdef enum Healpix_Ordering_Scheme:
        RING, NEST
    cdef cppclass nside_dummy:
        pass
    cdef nside_dummy SET_NSIDE

    cdef cppclass T_Healpix_Base[I]:
       int nside2order(I nside)
       I npix2nside(I npix)
       T_Healpix_Base()
       T_Healpix_Base(int order, Healpix_Ordering_Scheme scheme)
       T_Healpix_Base(I nside, Healpix_Ordering_Scheme scheme,
                     nside_dummy)
       double ring2z(I ring)
       I pix2ring(I pix)
       I xyf2pix(int ix, int iy, int face_num)
       void pix2xyf(I pix, int &ix, int &iy, int &face_num)
       I nest2ring (I pix)
       I ring2nest (I pix)
       I nest2peano (I pix)
       I peano2nest (I pix)
       I zphi2pix (double z, double phi)
       I ang2pix (pointing &ang)
       I vec2pix (vec3 &vec)
       void pix2zphi (I pix, double &z, double &phi)
       pointing pix2ang (I pix)
       vec3 pix2vec (I pix)
       void get_ring_info (I ring, I &startpix, I &ringpix,
                           double &costheta, double &sintheta, bool &shifted)
       void get_ring_info2 (I ring, I &startpix, I &ringpix,
                            double &theta, bool &shifted)

       void get_ring_info_small (I ring, I &startpix, I &ringpix,
                                 bool &shifted)
       void neighbors (I pix, fix_arr[I,val_8] &result)
       void get_interpol (pointing &ptg, fix_arr[I,val_4] &pix,
                          fix_arr[double,val_4] &wgt)

       int Order()
       I Nside()
       I Npix()
       Healpix_Ordering_Scheme Scheme()
       bool conformable (T_Healpix_Base &other)
       void swap (T_Healpix_Base &other)
       double max_pixrad()
       double max_pixrad(I ring)
       arr[int] swap_cycles()


def ringinfo(nside, np.ndarray[int, ndim=1] ring not None):
    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef Healpix_Ordering_Scheme scheme = NEST
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    num = ring.shape[0]
    cdef np.ndarray[int64, ndim=1] startpix = np.empty(num, dtype=np.int64)
    cdef np.ndarray[int64, ndim=1] ringpix = np.empty(num, dtype=np.int64)
    cdef np.ndarray[double, ndim=1] costheta = np.empty(num, dtype=np.float)
    cdef np.ndarray[double, ndim=1] sintheta = np.empty(num, dtype=np.float)
    cdef np.ndarray[bool, ndim=1] shifted = np.empty(num, dtype=np.bool)
    for i in range(num):
        hb.get_ring_info(ring[i], startpix[i], ringpix[i], costheta[i], sintheta[i], shifted[i])
    return startpix, ringpix, costheta, sintheta, shifted

def pix2ring(nside, np.ndarray[int64, ndim=1] pix not None, nest=False):
    if not isnsideok(nside):
        raise ValueError('Wrong nside value, must be a power of 2')
    cdef Healpix_Ordering_Scheme scheme
    if nest:
        scheme = NEST
    else:
        scheme = RING
    cdef T_Healpix_Base[int64] hb = T_Healpix_Base[int64](nside, scheme, SET_NSIDE)
    num = pix.shape[0]
    cdef np.ndarray[int64, ndim=1] ring = np.empty(num, dtype=np.int64)
    for i in range(num):
        ring[i] = hb.pix2ring(pix[i])
    return ring


cdef bool isnsideok(int nside):
    if nside < 0 or nside != 2**int(round(np.log2(nside))):
        return False
    else:
        return True



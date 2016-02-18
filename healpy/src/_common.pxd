import numpy as np
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector

ctypedef unsigned size_t
ctypedef size_t tsize

cdef extern from "stddef.h":
    ctypedef long ptrdiff_t
ctypedef ptrdiff_t tdiff

cdef extern from "datatypes.h":
    ctypedef int int64

cdef extern from "xcomplex.h":
    cdef cppclass xcomplex[T]:
       T re, im

cdef extern from "arr.h":    
    cdef cppclass arr[T]:
        arr()
        arr(T *ptr, tsize sz)
        void allocAndFill (tsize sz, T &inival)
        tsize size()
        T &operator[] (tsize n)
    cdef cppclass fix_arr[T, int]:
        pass

cdef extern from "rangeset.h":
    cdef cppclass rangeset[T]:
        tsize size()
        int64 nval()
        T ivbegin(tdiff i)
        T ivend(tdiff i)

cdef extern from "vec3.h":
    cdef cppclass vec3:
        double x, y, z
        vec3()
        vec3(double xc, double yc, double zc)
 
cdef extern from "pointing.h":
    cdef cppclass pointing:
        pointing()
        pointing(vec3 inp)
        void normalize()
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
       T_Healpix_Base(I nside, Healpix_Ordering_Scheme scheme, nside_dummy)
       void Set(int order, Healpix_Ordering_Scheme scheme)
       void SetNside(I nside, Healpix_Ordering_Scheme scheme)
       double ring2z(I ring)
       I pix2ring(I pix)
       I xyf2pix(int ix, int iy, int face_num)
       void pix2xyf(I pix, int &ix, int &iy, int &face_num)
       I nest2ring(I pix)
       I ring2nest(I pix)
       I nest2peano(I pix)
       I peano2nest(I pix)
       I zphi2pix(double z, double phi)
       I ang2pix(pointing &ang)
       I vec2pix(vec3 &vec)
       void pix2zphi(I pix, double &z, double &phi)
       pointing pix2ang(I pix)
       vec3 pix2vec(I pix)
       void get_ring_info(I ring, I &startpix, I &ringpix, double &costheta,
                          double &sintheta, bool &shifted)
       void get_ring_info2(I ring, I &startpix, I &ringpix, double &theta,
                           bool &shifted)
       void get_ring_info_small(I ring, I &startpix, I &ringpix, bool &shifted)
       void neighbors(I pix, fix_arr[I,val_8] &result)
       void get_interpol(pointing &ptg, fix_arr[I,val_4] &pix,
                         fix_arr[double,val_4] &wgt)
       int Order()
       I Nside()
       I Npix()
       Healpix_Ordering_Scheme Scheme()
       bool conformable(T_Healpix_Base &other)
       void swap(T_Healpix_Base &other)
       double max_pixrad()
       double max_pixrad(I ring)
       void boundaries(I pix, size_t step, vector[vec3] &out)
       arr[int] swap_cycles()
       void query_disc(pointing ptg, double radius, rangeset[I]& pixset) 
       void query_disc_inclusive(pointing ptg, double radius,
                                 rangeset[I]& pixset, int fact)
       void query_polygon(vector[pointing] vert, rangeset[I]& pixset) except +
       void query_polygon_inclusive(vector[pointing] vert,
                                    rangeset[I]& pixset, int fact) except +
       void query_strip(double theta1, double theta2, bool inclusive,
                        rangeset[I]& pixset)

cdef extern from "healpix_map.h":
    cdef cppclass Healpix_Map[T]:
        Healpix_Map()
        void Set(arr[T] &data, Healpix_Ordering_Scheme scheme)
        T average()
        void Add(T x)

cdef extern from "alm.h":
    cdef cppclass Alm[T]:
        Alm()
        Alm(int lmax_, int mmax_)
        void Set (int lmax_, int mmax_)
        void Set (arr[T] &data, int lmax_, int mmax_)
        tsize Num_Alms (int l, int m)

cdef inline Healpix_Map[double]* ndarray2map(np.ndarray[np.float64_t, ndim=1, mode='c'] array, Healpix_Ordering_Scheme scheme) except *:
    """ View a contiguous ndarray as a Healpix Map. """
    # To ensure that the output map is a view of the input array, the latter
    # is forced to be contiguous, of correct type and dimensions (otherwise, an
    # exception is raised).
    cdef arr[double] *a = new arr[double](&array[0], array.size)
    cdef Healpix_Map[double] *map = new Healpix_Map[double]()
    map.Set(a[0], scheme)
    del a # a does not own its buffer, so it won't be deallocated
    return map

cdef inline Alm[xcomplex[double]]* ndarray2alm(np.ndarray[np.complex128_t, ndim=1, mode='c'] array, int lmax, int mmax) except *:
    """ View a contiguous ndarray as an Alm. """
    cdef arr[xcomplex[double]] *a = new arr[xcomplex[double]](<xcomplex[double]*>&array[0], array.size)
    cdef Alm[xcomplex[double]] *alm = new Alm[xcomplex[double]]()
    alm.Set(a[0], lmax, mmax)
    del a
    return alm

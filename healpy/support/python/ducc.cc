#include "ducc0/infra/string_utils.cc"
#include "ducc0/infra/threading.cc"
#include "ducc0/math/pointing.cc"
#include "ducc0/math/geom_utils.cc"
#include "ducc0/math/space_filling.cc"
#include "ducc0/sht/sht.cc"
#include "ducc0/healpix/healpix_tables.cc"
#include "ducc0/healpix/healpix_base.cc"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <string>

#include "ducc0/healpix/healpix_base.h"
#include "ducc0/math/constants.h"
#include "ducc0/infra/string_utils.h"
#include "ducc0/math/geom_utils.h"
#include "ducc0/bindings/pybind_utils.h"

namespace ducc0 {

using namespace std;

namespace py = pybind11;

using shape_t = fmav_info::shape_t;

template<typename To1, typename Ti1, typename Func>
  py::object doStuff_i1o1(const py::array_t<Ti1> &in1_, Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto shp = in1.shape();
  auto out1_ = make_Pyarr<To1>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  in1.prepend_dim();
  out1.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<To1,1> iout1(out1);
  while (!iin1.done())
    {
    func(iin1, iout1);
    iin1.inc(); iout1.inc();
    }
  return (out1_.ndim()==0) ? py::cast(out1(0)) : py::object(out1_);
  }

template<typename To1, typename Ti1, typename Ti2, typename Func>
  py::object doStuff_i2o1(const py::array_t<Ti1> &in1_, const py::array_t<Ti2> &in2_, Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto in2 = to_fmav<Ti2>(in2_);
  auto shp = in2.bcast_shape(in1.shape());
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  auto out1_ = make_Pyarr<To1>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  out1.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<Ti2,1> iin2(in2);
  MavIter<To1,1> iout1(out1);
  while (!iin1.done())
    {
    func(iin1, iin2, iout1);
    iin1.inc(); iin2.inc(); iout1.inc();
    }
  return (out1_.ndim()==0) ? py::cast(out1(0)) : py::object(out1_);
  }

template<typename To1, typename To2, typename Ti1, typename Ti2, typename Func>
  py::tuple doStuff_i2o2(const py::array_t<Ti1> &in1_, const py::array_t<Ti2> &in2_, Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto in2 = to_fmav<Ti2>(in2_);
  auto shp = in2.bcast_shape(in1.shape());
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  auto out1_ = make_Pyarr<To1>(shp);
  auto out2_ = make_Pyarr<To2>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  auto out2 = to_fmav<To2>(out2_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  out1.prepend_dim();
  out2.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<Ti2,1> iin2(in2);
  MavIter<To1,1> iout1(out1);
  MavIter<To2,1> iout2(out2);
  while (!iin1.done())
    {
    func(iin1, iin2, iout1, iout2);
    iin1.inc(); iin2.inc(); iout1.inc(); iout2.inc();
    }
  if (out1_.ndim()==0)
    return py::make_tuple(py::cast(out1(0)), py::cast(out2(0)));
  else
    return py::make_tuple(out1_, out2_);
  }

template<typename To1, typename To2, typename To3, typename Ti1, typename Ti2, typename Func>
  py::tuple doStuff_i2o3(const py::array_t<Ti1> &in1_, const py::array_t<Ti2> &in2_, Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto in2 = to_fmav<Ti2>(in2_);
  auto shp = in2.bcast_shape(in1.shape());
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  auto out1_ = make_Pyarr<To1>(shp);
  auto out2_ = make_Pyarr<To2>(shp);
  auto out3_ = make_Pyarr<To3>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  auto out2 = to_fmav<To2>(out2_,true);
  auto out3 = to_fmav<To3>(out3_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  out1.prepend_dim();
  out2.prepend_dim();
  out3.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<Ti2,1> iin2(in2);
  MavIter<To1,1> iout1(out1);
  MavIter<To2,1> iout2(out2);
  MavIter<To3,1> iout3(out3);
  while (!iin1.done())
    {
    func(iin1, iin2, iout1, iout2, iout3);
    iin1.inc(); iin2.inc(); iout1.inc(); iout2.inc(); iout3.inc();
    }
  if (out1_.ndim()==0)
    return py::make_tuple(py::cast(out1(0)), py::cast(out2(0)), py::cast(out3(0)));
  else
    return py::make_tuple(out1_, out2_, out3_);
  }

template<typename To1, typename Ti1, typename Ti2, typename Ti3,
  typename Func>
  py::object doStuff_i3o1(const py::array_t<Ti1> &in1_, const py::array_t<Ti2> &in2_, const py::array_t<Ti3> &in3_, Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto in2 = to_fmav<Ti2>(in2_);
  auto in3 = to_fmav<Ti3>(in3_);
  auto shp = in3.bcast_shape(in2.bcast_shape(in1.shape()));
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  in3.bcast_to_shape(shp);
  auto out1_ = make_Pyarr<To1>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  in3.prepend_dim();
  out1.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<Ti2,1> iin2(in2);
  MavIter<Ti3,1> iin3(in3);
  MavIter<To1,1> iout1(out1);
  while (!iin1.done())
    {
    func(iin1, iin2, iin3, iout1);
    iin1.inc(); iin2.inc(); iin3.inc(); iout1.inc();
    }
  return (out1_.ndim()==0) ? py::cast(out1(0)) : py::object(out1_);
  }

template<typename To1, typename Ti1, typename Ti2, typename Ti3, typename Ti4,
  typename Func>
  py::object doStuff_i4o1(const py::array_t<Ti1> &in1_, const py::array_t<Ti2> &in2_, const py::array_t<Ti3> &in3_, const py::array_t<Ti4> &in4_,
  Func func)
  {
  auto in1 = to_fmav<Ti1>(in1_);
  auto in2 = to_fmav<Ti2>(in2_);
  auto in3 = to_fmav<Ti3>(in3_);
  auto in4 = to_fmav<Ti4>(in4_);
  auto shp = in4.bcast_shape(in3.bcast_shape(in2.bcast_shape(in1.shape())));
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  in3.bcast_to_shape(shp);
  in4.bcast_to_shape(shp);
  auto out1_ = make_Pyarr<To1>(shp);
  auto out1 = to_fmav<To1>(out1_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  in3.prepend_dim();
  in4.prepend_dim();
  out1.prepend_dim();
  MavIter<Ti1,1> iin1(in1);
  MavIter<Ti2,1> iin2(in2);
  MavIter<Ti3,1> iin3(in3);
  MavIter<Ti4,1> iin4(in4);
  MavIter<To1,1> iout1(out1);
  while (!iin1.done())
    {
    func(iin1, iin2, iin3, iin4, iout1);
    iin1.inc(); iin2.inc(); iin3.inc(); iin4.inc(); iout1.inc();
    }
  return (out1_.ndim()==0) ? py::cast(out1(0)) : py::object(out1_);
  }

template<typename Func>
  py::tuple doStuff_interpol(const py::array_t<int64_t> &in1_, const py::array_t<double> &in2_, const py::array_t<double> &in3_, Func func)
  {
  auto in1 = to_fmav<int64_t>(in1_);
  auto in2 = to_fmav<double>(in2_);
  auto in3 = to_fmav<double>(in3_);
  auto shp = in3.bcast_shape(in2.bcast_shape(in1.shape()));
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  in3.bcast_to_shape(shp);
  auto shp2 = shp;
  shp2.push_back(4);
  auto out1_ = make_Pyarr<int64_t>(shp2);
  auto out1 = to_fmav<int64_t>(out1_,true);
  auto out2_ = make_Pyarr<double>(shp2);
  auto out2 = to_fmav<double>(out2_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  in3.prepend_dim();
  out1.prepend_dim();
  out2.prepend_dim();
  MavIter<int64_t,1> iin1(in1);
  MavIter<double,1> iin2(in2);
  MavIter<double,1> iin3(in3);
  MavIter<int64_t,2> iout1(out1);
  MavIter<double,2> iout2(out2);
  while (!iin1.done())
    {
    func(iin1, iin2, iin3, iout1, iout2);
    iin1.inc(); iin2.inc(); iin3.inc(); iout1.inc(); iout2.inc();
    }
  return py::make_tuple(out1_, out2_);
  }

template<typename Func>
  py::array doStuff_neighbors(const py::array_t<int64_t> &in1_, const py::array_t<int64_t> &in2_, Func func)
  {
  auto in1 = to_fmav<int64_t>(in1_);
  auto in2 = to_fmav<int64_t>(in2_);
  auto shp = in2.bcast_shape(in1.shape());
  in1.bcast_to_shape(shp);
  in2.bcast_to_shape(shp);
  auto shp2 = shp;
  shp2.push_back(8);
  auto out1_ = make_Pyarr<int64_t>(shp2);
  auto out1 = to_fmav<int64_t>(out1_,true);
  in1.prepend_dim();
  in2.prepend_dim();
  out1.prepend_dim();
  MavIter<int64_t,1> iin1(in1);
  MavIter<int64_t,1> iin2(in2);
  MavIter<int64_t,2> iout1(out1);
  while (!iin1.done())
    {
    func(iin1, iin2, iout1);
    iin1.inc(); iin2.inc(); iout1.inc();
    }
  return out1_;
  }

py::object ang2pix(const py::array_t<int64_t> &nside, const py::array_t<double> &theta, const py::array_t<double> &phi, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i3o1<int64_t>(nside, theta, phi, [&]
    (const auto &nside_, const auto &theta_, const auto &phi_, auto &pix)
    {
    for (size_t i=0; i<pix.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      pix.v(i) = base.ang2pix(pointing(theta_(i),phi_(i)));
      }
    });
  }

py::object pix2ang(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &ipix, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i2o2<double, double>(nside, ipix, [&]
    (const auto &nside_, const auto &ipix_, auto &theta, auto &phi)
    {
    for (size_t i=0; i<theta.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      auto tmp = base.pix2ang(ipix_(i));
      theta.v(i) = tmp.theta;
      phi.v(i) = tmp.phi;
      }
    });
  }

py::object xyf2pix(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &x,
  const py::array_t<int64_t> &y, const py::array_t<int64_t> &f, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i4o1<int64_t>(nside, x, y, f, [&]
    (const auto &nside_, const auto &x_, const auto &y_, const auto &f_, auto &ipix)
    {
    for (size_t i=0; i<ipix.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      ipix.v(i) = base.xyf2pix(x_(i), y_(i), f_(i));
      }
    });
  }

py::tuple pix2xyf(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &ipix, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i2o3<int64_t, int64_t, int64_t>(nside, ipix, [&]
    (const auto &nside_, const auto &ipix_, auto &x, auto &y, auto &f)
    {
    for (size_t i=0; i<x.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      int tx, ty, tf;
      base.pix2xyf(ipix_(i), tx, ty, tf);
      x.v(i) = tx;
      y.v(i) = ty;
      f.v(i) = tf;
      }
    });
  }

py::object vec2pix(const py::array_t<int64_t> &nside, const py::array_t<double> &x,
  const py::array_t<double> &y, const py::array_t<double> &z, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i4o1<int64_t>(nside, x, y, z, [&]
    (const auto &nside_, const auto &x_, const auto &y_, const auto &z_, auto &ipix)
    {
    for (size_t i=0; i<ipix.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      ipix.v(i) = base.vec2pix(vec3(x_(i), y_(i), z_(i)));
      }
    });
  }

py::tuple pix2vec(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &ipix, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_i2o3<double, double, double>(nside, ipix, [&]
    (const auto &nside_, const auto &ipix_, auto &x, auto &y, auto &z)
    {
    for (size_t i=0; i<x.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      auto tmp = base.pix2vec(ipix_(i));
      x.v(i) = tmp.x;
      y.v(i) = tmp.y;
      z.v(i) = tmp.z;
      }
    });
  }
py::object nest2ring(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &inest)
  {
  Healpix_Base2 base(1, NEST);
  return doStuff_i2o1<int64_t>(nside, inest, [&]
    (const auto &nside_, const auto &inest_, auto &iring)
    {
    for (size_t i=0; i<iring.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      iring.v(i) = base.nest2ring(inest_(i));
      }
    });
  }

py::object ring2nest(const py::array_t<int64_t> &nside, const py::array_t<int64_t> &iring)
  {
  Healpix_Base2 base(1, NEST);
  return doStuff_i2o1<int64_t>(nside, iring, [&]
    (const auto &nside_, const auto &iring_, auto &inest)
    {
    for (size_t i=0; i<inest.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      inest.v(i) = base.ring2nest(iring_(i));
      }
    });
  }

py::object max_pixrad(const py::array_t<int64_t> &nside)
  {
  Healpix_Base2 base(1, NEST);
  return doStuff_i1o1<double>(nside, [&]
    (const auto &nside_, auto &maxrad)
    {
    for (size_t i=0; i<maxrad.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      maxrad.v(i) = base.max_pixrad();
      }
    });
  }

py::tuple get_interpol(const py::array_t<int64_t> &nside,
  const py::array_t<double> &theta, const py::array_t<double> &phi, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_interpol(nside, theta, phi, [&]
    (const auto &nside_, auto &theta_, auto &phi_, auto &pix, auto &wgt)
    {
    array<int64_t,4> lpix;
    array<double,4> lwgt;
    for (size_t i=0; i<pix.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      pointing tmp(theta_(i), phi_(i));
      tmp.normalize();
      base.get_interpol(tmp, lpix, lwgt);
      pix.v(i,0) = lpix[0];
      pix.v(i,1) = lpix[1];
      pix.v(i,2) = lpix[2];
      pix.v(i,3) = lpix[3];
      wgt.v(i,0) = lwgt[0];
      wgt.v(i,1) = lwgt[1];
      wgt.v(i,2) = lwgt[2];
      wgt.v(i,3) = lwgt[3];
      }
    });
  }

py::array get_neighbors(const py::array_t<int64_t> &nside,
  const py::array_t<int64_t> &ipix, bool nest)
  {
  Healpix_Base2 base(1, nest ? NEST : RING);
  return doStuff_neighbors(nside, ipix, [&]
    (const auto &nside_, auto &ipix_, auto &neigh)
    {
    array<int64_t,8> lneigh;
    for (size_t i=0; i<neigh.shape(0); ++i)
      {
      if (base.Nside()!=nside_(i))
        base.SetNside(nside_(i), base.Scheme());
      base.neighbors(ipix_(i), lneigh);
      for (size_t j=0; j<8; ++j)
        neigh.v(i,j) = lneigh[j];
      }
    });
  }

}

using namespace ducc0;

PYBIND11_MODULE(_support, m)
  {
  using namespace pybind11::literals;

  m.attr("UNSEEN") = -1.6375e30;
  m.def("ang2pix", &ang2pix, "nside"_a, "theta"_a, "phi"_a, "nest"_a);
  m.def("pix2ang", &pix2ang, "nside"_a, "ipix"_a, "nest"_a);
  m.def("vec2pix", &vec2pix, "nside"_a, "x"_a, "y"_a, "z"_a, "nest"_a);
  m.def("pix2vec", &pix2vec, "nside"_a, "ipix"_a, "nest"_a);
  m.def("xyf2pix", &xyf2pix, "nside"_a, "x"_a, "y"_a, "f"_a, "nest"_a);
  m.def("pix2xyf", &pix2xyf, "nside"_a, "ipix"_a, "nest"_a);
  m.def("nest2ring", &nest2ring, "nside"_a, "ipix"_a);
  m.def("ring2nest", &ring2nest, "nside"_a, "ipix"_a);
  m.def("max_pixrad", &max_pixrad, "nside"_a);
  m.def("get_interpol", &get_interpol, "nside"_a, "theta"_a, "phi"_a, "nest"_a);
  m.def("get_neighbors", &get_neighbors, "nside"_a, "ipix"_a, "nest"_a);
  }

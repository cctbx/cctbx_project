/*
 * pixel_to_millimeter.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/pixel_to_millimeter.h>

namespace dxtbx { namespace model { namespace boost_python {

  using namespace boost::python;

  struct PxMmStrategyPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const PxMmStrategy& obj) {
      return boost::python::make_tuple();
    }
  };

  struct ParallaxCorrectedPxMmStrategyPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const ParallaxCorrectedPxMmStrategy& obj) {
      return boost::python::make_tuple(obj.mu(), obj.t0());
    }
  };

  struct OffsetParallaxCorrectedPxMmStrategyPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const OffsetParallaxCorrectedPxMmStrategy& obj) {
      return boost::python::make_tuple(obj.mu(), obj.t0(), obj.dx(), obj.dy());
    }
  };

  void export_pixel_to_millimeter()
  {
    class_<PxMmStrategy, boost::noncopyable>("PxMmStrategy", no_init)
      .def("to_millimeter", &PxMmStrategy::to_millimeter, (
        arg("panel"), arg("xy")))
      .def("to_pixel", &PxMmStrategy::to_pixel, (
        arg("panel"), arg("xy")))
      .def("__str__", &PxMmStrategy::strategy_name);

    class_<SimplePxMmStrategy, bases<PxMmStrategy> >("SimplePxMmStrategy")
      .def_pickle(PxMmStrategyPickleSuite())
      .def("__str__", &PxMmStrategy::strategy_name);

    class_<ParallaxCorrectedPxMmStrategy, bases<SimplePxMmStrategy> >(
      "ParallaxCorrectedPxMmStrategy", no_init)
      .def(init<double, double>((arg("mu"), arg("t0"))))
      .def("mu",&ParallaxCorrectedPxMmStrategy::mu)
      .def("t0",&ParallaxCorrectedPxMmStrategy::t0)
      .def_pickle(ParallaxCorrectedPxMmStrategyPickleSuite())
      .def("__str__", &PxMmStrategy::strategy_name);

    class_<OffsetParallaxCorrectedPxMmStrategy, bases<ParallaxCorrectedPxMmStrategy> >(
      "OffsetParallaxCorrectedPxMmStrategy", no_init)
      .def(init<
          double,
          double,
          scitbx::af::const_ref<double, scitbx::af::c_grid<2> >,
          scitbx::af::const_ref<double, scitbx::af::c_grid<2> > >(
            (arg("mu"),
             arg("t0"),
             arg("dx"),
             arg("dy"))))
      .def("dx", &OffsetParallaxCorrectedPxMmStrategy::dx)
      .def("dy", &OffsetParallaxCorrectedPxMmStrategy::dy)
      .def_pickle(OffsetParallaxCorrectedPxMmStrategyPickleSuite())
      .def("__str__", &PxMmStrategy::strategy_name);

    register_ptr_to_python<shared_ptr<PxMmStrategy> >();
    register_ptr_to_python<shared_ptr<SimplePxMmStrategy> >();
    register_ptr_to_python<shared_ptr<ParallaxCorrectedPxMmStrategy> >();
    register_ptr_to_python<shared_ptr<OffsetParallaxCorrectedPxMmStrategy> >();
  }

}}} // namespace = dxtbx::model::boost_python

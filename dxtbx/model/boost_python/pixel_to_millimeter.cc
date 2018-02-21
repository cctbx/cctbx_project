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

  static
  OffsetPxMmStrategy* OffsetPxMmStrategy_init(
      scitbx::af::versa< double, scitbx::af::flex_grid<> > dx,
      scitbx::af::versa< double, scitbx::af::flex_grid<> > dy) {

    DXTBX_ASSERT(dx.accessor().all().size() == 2);
    DXTBX_ASSERT(dy.accessor().all().size() == 2);
    DXTBX_ASSERT(dx.accessor().all().all_eq(dy.accessor().all()));

    std::size_t ysize = dx.accessor().all()[0];
    std::size_t xsize = dx.accessor().all()[1];

    scitbx::af::c_grid<2> grid(ysize, xsize);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dx2(dx.handle(), grid);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dy2(dy.handle(), grid);

    return  new OffsetPxMmStrategy(dx2, dy2);
  }

  static
  OffsetParallaxCorrectedPxMmStrategy* OffsetParallaxCorrectedPxMmStrategy_init(
      double mu,
      double t0,
      scitbx::af::versa< double, scitbx::af::flex_grid<> > dx,
      scitbx::af::versa< double, scitbx::af::flex_grid<> > dy) {

    DXTBX_ASSERT(dx.accessor().all().size() == 2);
    DXTBX_ASSERT(dy.accessor().all().size() == 2);
    DXTBX_ASSERT(dx.accessor().all().all_eq(dy.accessor().all()));

    std::size_t ysize = dx.accessor().all()[0];
    std::size_t xsize = dx.accessor().all()[1];

    scitbx::af::c_grid<2> grid(ysize, xsize);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dx2(dx.handle(), grid);
    scitbx::af::versa< double, scitbx::af::c_grid<2> > dy2(dy.handle(), grid);

    return  new OffsetParallaxCorrectedPxMmStrategy(mu, t0, dx2, dy2);
  }


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

  struct OffsetPxMmStrategyPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const OffsetPxMmStrategy& obj) {
      return boost::python::make_tuple(obj.dx(), obj.dy());
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
      .def("name", &PxMmStrategy::name)
      .def("__str__", &PxMmStrategy::strategy_name)
      ;

    class_<SimplePxMmStrategy, bases<PxMmStrategy> >("SimplePxMmStrategy")
      .def_pickle(PxMmStrategyPickleSuite())
      ;

    class_<ParallaxCorrectedPxMmStrategy, bases<SimplePxMmStrategy> >(
      "ParallaxCorrectedPxMmStrategy", no_init)
      .def(init<double, double>((arg("mu"), arg("t0"))))
      .def("mu",&ParallaxCorrectedPxMmStrategy::mu)
      .def("t0",&ParallaxCorrectedPxMmStrategy::t0)
      .def_pickle(ParallaxCorrectedPxMmStrategyPickleSuite())
      ;

    class_<OffsetPxMmStrategy, bases<PxMmStrategy> >(
      "OffsetPxMmStrategy", no_init)
      .def("__init__",
          make_constructor(
          &OffsetPxMmStrategy_init,
          default_call_policies(), (
            (arg("dx"),
             arg("dy")))))
      .def("dx", &OffsetPxMmStrategy::dx)
      .def("dy", &OffsetPxMmStrategy::dy)
      .def_pickle(OffsetPxMmStrategyPickleSuite())
      ;

    class_<OffsetParallaxCorrectedPxMmStrategy, bases<ParallaxCorrectedPxMmStrategy> >(
      "OffsetParallaxCorrectedPxMmStrategy", no_init)
      .def("__init__",
          make_constructor(
          &OffsetParallaxCorrectedPxMmStrategy_init,
          default_call_policies(), (
            (arg("mu"),
             arg("t0"),
             arg("dx"),
             arg("dy")))))
      .def("dx", &OffsetParallaxCorrectedPxMmStrategy::dx)
      .def("dy", &OffsetParallaxCorrectedPxMmStrategy::dy)
      .def_pickle(OffsetParallaxCorrectedPxMmStrategyPickleSuite())
      ;

    register_ptr_to_python<shared_ptr<PxMmStrategy> >();
    register_ptr_to_python<shared_ptr<SimplePxMmStrategy> >();
    register_ptr_to_python<shared_ptr<ParallaxCorrectedPxMmStrategy> >();
    register_ptr_to_python<shared_ptr<OffsetPxMmStrategy> >();
    register_ptr_to_python<shared_ptr<OffsetParallaxCorrectedPxMmStrategy> >();
  }

}}} // namespace = dxtbx::model::boost_python

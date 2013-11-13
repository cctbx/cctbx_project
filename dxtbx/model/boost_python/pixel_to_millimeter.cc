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
      return boost::python::make_tuple(obj.attenuation_length());
    }  
  };

  void export_pixel_to_millimeter()
  {
    class_<PxMmStrategy, boost::noncopyable>("PxMmStrategy", no_init)
      .def("to_millimeter", &PxMmStrategy::to_millimeter, (
        arg("panel"), arg("xy")))
      .def("to_pixel", &PxMmStrategy::to_pixel, (
        arg("panel"), arg("xy")));
          
    class_<SimplePxMmStrategy, bases<PxMmStrategy> >("SimplePxMmStrategy")
      .def_pickle(PxMmStrategyPickleSuite());

    class_<ParallaxCorrectedPxMmStrategy, bases<SimplePxMmStrategy> >(
      "ParallaxCorrectedPxMmStrategy", no_init)
      .def(init<double>((arg("la"))))
      .def("attenuation_length", 
        &ParallaxCorrectedPxMmStrategy::attenuation_length)
      .def_pickle(ParallaxCorrectedPxMmStrategyPickleSuite());

    register_ptr_to_python<shared_ptr<PxMmStrategy> >();
    register_ptr_to_python<shared_ptr<SimplePxMmStrategy> >();
    register_ptr_to_python<shared_ptr<ParallaxCorrectedPxMmStrategy> >();
  }

}}} // namespace = dxtbx::model::boost_python

#include <smtbx/refinement/weighting_schemes.h>

#include <boost/python/class.hpp>


namespace smtbx { namespace refinement { namespace boost_python {

  template<template<typename FloatType> class WeightingScheme>
  struct weighting_scheme_class
    : boost::python::class_< WeightingScheme<double> >
  {
    typedef WeightingScheme<double> wt;
    typedef boost::python::class_<wt> base_t;

    weighting_scheme_class(char const *name)
      : base_t(name, boost::python::no_init)
    {
      using namespace boost::python;
      def("__call__", &wt::operator(),
          (arg("fo_sq"), arg("sigma"), arg("fc_sq")));
    }
  };

  struct mainstream_shelx_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      weighting_scheme_class<
        mainstream_shelx_weighting
      >("mainstream_shelx_weighting")
        .def(init<optional<double, double> >((arg("a"), arg("b"))))
        ;
    }
  };

  struct unit_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      weighting_scheme_class<unit_weighting>("unit_weighting")
        .def(init<>())
        ;
    }
  };


  void wrap_weighting_schemes() {
    mainstream_shelx_weighting_wrapper::wrap();
    unit_weighting_wrapper::wrap();
  }


}}}

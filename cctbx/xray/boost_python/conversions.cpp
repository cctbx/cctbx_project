#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/conversions.h>
#include <boost/python/class.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103000
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#endif

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  template<template<class> class FsqAsF>
  struct array_f_sq_as_f_wrappers
  {
    typedef array_f_sq_as_f<FsqAsF> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
#if BOOST_VERSION >= 103000
      typedef return_value_policy<return_by_value> rbv;
#endif
      class_<w_t>(python_name, no_init)
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  optional<double const&> >())
        .def(init<af::const_ref<double> const&>())
#if BOOST_VERSION >= 103000
        .add_property("f", make_getter(&w_t::f, rbv()))
        .add_property("sigma_f", make_getter(&w_t::sigma_f, rbv()))
#else
        .def_readonly("f", &w_t::f)
        .def_readonly("sigma_f", &w_t::sigma_f)
#endif
      ;
    }
  };

  struct array_f_as_f_sq_wrappers
  {
    typedef array_f_as_f_sq<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
#if BOOST_VERSION >= 103000
      typedef return_value_policy<return_by_value> rbv;
#endif
      class_<w_t>("array_f_as_f_sq", no_init)
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .def(init<af::const_ref<double> const&>())
#if BOOST_VERSION >= 103000
        .add_property("f_sq", make_getter(&w_t::f_sq, rbv()))
        .add_property("sigma_f_sq", make_getter(&w_t::sigma_f_sq, rbv()))
#else
        .def_readonly("f_sq", &w_t::f_sq)
        .def_readonly("sigma_f_sq", &w_t::sigma_f_sq)
#endif
      ;
    }
  };

} // namespace <anoymous>

  void wrap_conversions()
  {
    array_f_sq_as_f_wrappers<cctbx::xray::f_sq_as_f_xtal_3_7>::wrap(
      "array_f_sq_as_f_xtal_3_7");
    array_f_sq_as_f_wrappers<cctbx::xray::f_sq_as_f_crystals>::wrap(
      "array_f_sq_as_f_crystals");
    array_f_as_f_sq_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python

/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

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

  struct array_f_sq_as_f_wrappers
  {
    typedef array_f_sq_as_f<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
#if BOOST_VERSION >= 103000
      typedef return_value_policy<return_by_value> rbv;
#endif
      class_<w_t>("array_f_sq_as_f", no_init)
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
    array_f_sq_as_f_wrappers::wrap();
    array_f_as_f_sq_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python

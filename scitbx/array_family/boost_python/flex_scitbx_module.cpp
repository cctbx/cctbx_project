/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/arraytbx/flexmodule.cpp (rwgk)
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#include <scitbx/math/linear_regression.h>
#include <scitbx/array_family/boost_python/ref_c_grid_flex_conversions.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/array_family/boost_python/small_conversions.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_grid();
  void wrap_flex_bool();
  void wrap_flex_size_t();
  void wrap_flex_int();
  void wrap_flex_long();
  void wrap_flex_float();
  void wrap_flex_double();
  void wrap_flex_complex_double();
  void wrap_flex_std_string();

namespace {

  struct linear_regression_core_wrappers
  {
    typedef math::linear_regression_core<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("linear_regression_core", no_init)
        .def("is_well_defined", &w_t::is_well_defined)
        .def("y_intercept", &w_t::y_intercept)
        .def("slope", &w_t::slope)
        .def("cc", &w_t::cc)
      ;
    }
  };

  struct linear_regression_wrappers
  {
    typedef math::linear_regression<> w_t;
    typedef w_t::float_type float_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<math::linear_regression_core<> > >(
            "linear_regression", no_init)
        .def(init<af::const_ref<float_t> const&,
                  af::const_ref<float_t> const&,
                  optional<float_t const&> >());
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    wrap_flex_grid();
    register_flex_grid_default_index_type_conversions();

    linear_regression_core_wrappers::wrap();
    linear_regression_wrappers::wrap();

    wrap_flex_bool();
    wrap_flex_size_t();
    wrap_flex_int();
    wrap_flex_long();
    wrap_flex_float();
    wrap_flex_double();
    wrap_flex_complex_double();
    wrap_flex_std_string();

    default_ref_c_grid_flex_conversions<int>();
    default_ref_c_grid_flex_conversions<long>();
    default_ref_c_grid_flex_conversions<float>();
    default_ref_c_grid_flex_conversions<double>();
    default_ref_c_grid_flex_conversions<std::complex<double> >();
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(flex_scitbx_ext)
{
  scitbx::af::boost_python::init_module();
}

/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/arraytbx/flexmodule.cpp (rwgk)
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/linear_regression.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
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

  void register_scitbx_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_fixed_size<int3>();
    tuple_mapping_fixed_size<int9>();
    tuple_mapping_fixed_size<long3>();
    tuple_mapping_fixed_size<double2>();
    tuple_mapping_fixed_size<double3>();
    tuple_mapping_fixed_size<double6>();
    tuple_mapping_fixed_size<double9>();

    tuple_mapping_fixed_capacity<flex_grid_default_index_type>();
  }

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
        .def("correlation", &w_t::correlation)
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

    register_scitbx_tuple_mappings();

    wrap_flex_grid();

    wrap_flex_bool();
    wrap_flex_size_t();
    wrap_flex_int();
    wrap_flex_long();
    wrap_flex_float();
    wrap_flex_double();
    wrap_flex_complex_double();
    wrap_flex_std_string();

    default_c_grid_flex_conversions<int>();
    default_c_grid_flex_conversions<long>();
    default_c_grid_flex_conversions<float>();
    default_c_grid_flex_conversions<double>();
    default_c_grid_flex_conversions<std::complex<double> >();

    linear_regression_core_wrappers::wrap();
    linear_regression_wrappers::wrap();
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(flex_scitbx_ext)
{
  scitbx::af::boost_python::init_module();
}

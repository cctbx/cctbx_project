#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/linear_regression.h>
#include <scitbx/math/linear_correlation.h>
#include <scitbx/math/gaussian/sum.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/boost_python/slice.h>
#include <boost/python/module.hpp>
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
  void wrap_flex_vec3_double();

  void wrap_flex_sort();
  void wrap_flex_histogram();
  void wrap_flex_mean_and_variance();
  void wrap_flex_linear_interpolation();

  void wrap_loops();

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

    tuple_mapping_fixed_size<tiny<std::size_t,3> >();

    tuple_mapping_fixed_size<tiny<int,24> >(); // scitbx/math/golay.h

    tuple_mapping_fixed_size<vec3<int> >();
    tuple_mapping_fixed_size<mat3<int> >();
    tuple_mapping_fixed_size<vec3<double> >();
    tuple_mapping_fixed_size<mat3<double> >();
    tuple_mapping_fixed_size<sym_mat3<double> >();

    tuple_mapping_fixed_capacity<flex_grid_default_index_type>();
    tuple_mapping_fixed_capacity<
      small<double, math::gaussian::sum<double>::max_n_terms> >();
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

  struct linear_correlation_wrappers
  {
    typedef math::linear_correlation<> w_t;
    typedef w_t::float_type float_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("linear_correlation", no_init)
        .def(init<af::const_ref<float_t> const&,
                  af::const_ref<float_t> const&,
                  optional<float_t const&> >())
        .def("is_well_defined", &w_t::is_well_defined)
        .def("n", &w_t::n)
        .def("mean_x", &w_t::mean_x)
        .def("mean_y", &w_t::mean_y)
        .def("numerator", &w_t::numerator)
        .def("sum_denominator_x", &w_t::sum_denominator_x)
        .def("sum_denominator_y", &w_t::sum_denominator_y)
        .def("denominator", &w_t::denominator)
        .def("coefficient", &w_t::coefficient)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    register_scitbx_tuple_mappings();

    scitbx::boost_python::slice_from_python();
    wrap_flex_grid();

    wrap_flex_bool();
    wrap_flex_size_t();
    wrap_flex_int();
    wrap_flex_long();
    wrap_flex_float();
    wrap_flex_double();
    wrap_flex_complex_double();
    wrap_flex_std_string();
    wrap_flex_vec3_double();

    default_c_grid_flex_conversions<int>();
    default_c_grid_flex_conversions<long>();
    default_c_grid_flex_conversions<float>();
    default_c_grid_flex_conversions<double>();
    default_c_grid_flex_conversions<std::complex<double> >();

    wrap_flex_sort();
    wrap_flex_histogram();
    wrap_flex_mean_and_variance();
    wrap_flex_linear_interpolation();

    wrap_loops();

    linear_regression_core_wrappers::wrap();
    linear_regression_wrappers::wrap();
    linear_correlation_wrappers::wrap();
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_array_family_flex_ext)
{
  scitbx::af::boost_python::init_module();
}

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/linear_regression.h>
#include <scitbx/math/linear_correlation.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/sym_mat2.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <scitbx/array_family/boost_python/owning_ref_conversions.h>
#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/boost_python/slice.h>
#include <boost_adaptbx/optional_conversions.h>
#include <boost_adaptbx/optional_copy.h>
#include <boost_adaptbx/type_id_eq.h>
#include <boost/optional.hpp>
#include <boost/rational.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <vector>

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
  void wrap_flex_vec2_double();
  void wrap_flex_sym_mat3_double();
  void wrap_flex_tiny_size_t_2();

  void wrap_flex_random();
  void wrap_flex_sort();
  void wrap_flex_histogram();
  void wrap_flex_mean_and_variance();
  void wrap_flex_median_statistics();
  void wrap_flex_linear_interpolation();

  void wrap_loops();
  void wrap_empty_container_sizes();

namespace {

  void register_scitbx_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_fixed_size<tiny<bool, 3> >();
    tuple_mapping_fixed_size<tiny<int, 2> >();
    tuple_mapping_fixed_size<int3>();
    tuple_mapping_fixed_size<int9>();
    tuple_mapping_fixed_size<long3>();
    tuple_mapping_fixed_size<double2>();
    tuple_mapping_fixed_size<double3>();
    tuple_mapping_fixed_size<double4>();
    tuple_mapping_fixed_size<double6>();
    tuple_mapping_fixed_size<double9>();
    tuple_mapping_fixed_size<tiny<std::string, 2> >();
    tuple_mapping_fixed_size<tiny<std::string, 3> >();
    tuple_mapping_fixed_size<tiny<std::string, 4> >();

    tuple_mapping_fixed_capacity<small<int, 3> >();
    tuple_mapping_fixed_capacity<small<int, 10> >();
    tuple_mapping_fixed_capacity<small<unsigned, 2> >();
    tuple_mapping_fixed_capacity<small<unsigned, 3> >();
    tuple_mapping_fixed_capacity<small<unsigned, 6> >();
    tuple_mapping_fixed_capacity<small<std::size_t, 5> >();
    tuple_mapping_fixed_capacity<small<double, 3> >();
    tuple_mapping_fixed_capacity<small<double, 6> >();
    // scitbx/math/gaussian/sum.h SCITBX_MATH_GAUSSIAN_SUM_MAX_N_TERMS
    tuple_mapping_fixed_capacity<small<double, 10> >();

    tuple_mapping_fixed_size<tiny<int, 12> >();
    tuple_mapping_fixed_size<tiny<int, 24> >(); // scitbx/math/golay.h

    tuple_mapping_fixed_size<tiny<double, 12> >();

    tuple_mapping_fixed_size<tiny<unsigned, 2> >();
    tuple_mapping_fixed_size<tiny<unsigned, 3> >();
    tuple_mapping_fixed_size<tiny<unsigned, 4> >();
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
    tuple_mapping_fixed_size<tiny<std::size_t, 2> >();
    tuple_mapping_fixed_size<tiny<std::size_t, 3> >();
    tuple_mapping_fixed_size<tiny<std::size_t, 4> >();
#endif

    tuple_mapping_fixed_size<vec3<unsigned int> >();
    tuple_mapping_fixed_size<vec3<unsigned long> >();
    tuple_mapping_fixed_size<vec3<int> >();
    tuple_mapping_fixed_size<mat3<int> >();
    tuple_mapping_fixed_size<vec3<double> >();
    tuple_mapping_fixed_size<mat3<double> >();
    tuple_mapping_fixed_size<sym_mat3<double> >();
    tuple_mapping_fixed_size<vec2<unsigned int> >();
    tuple_mapping_fixed_size<vec2<unsigned long> >();
    tuple_mapping_fixed_size<vec2<int> >();
    tuple_mapping_fixed_size<mat2<int> >();
    tuple_mapping_fixed_size<vec2<double> >();
    tuple_mapping_fixed_size<mat2<double> >();
    tuple_mapping_fixed_size<sym_mat2<double> >();

    tuple_mapping_fixed_size<mat3<boost::rational<int> > >();

    tuple_mapping_fixed_size<tiny<vec3<double>, 2> >();
    tuple_mapping_fixed_size<tiny<vec3<double>, 3> >();
    tuple_mapping_fixed_size<tiny<vec3<double>, 4> >();
    tuple_mapping_fixed_size<tiny<vec2<double>, 2> >();

    tuple_mapping_fixed_capacity<flex_grid_default_index_type>();
    tuple_mapping_fixed_capacity<small<vec3<int>, 3> >();
    tuple_mapping_fixed_capacity<small<vec2<int>, 2> >();
  }

  af::shared<std::size_t>
  slice_indices(
    std::size_t array_size,
    scitbx::boost_python::slice const& python_slice)
  {
    scitbx::boost_python::adapted_slice a_sl(python_slice, array_size);
    af::shared<std::size_t> result(af::reserve(a_sl.size));
    for(long i=a_sl.start;i!=a_sl.stop;i+=a_sl.step) {
      result.push_back(static_cast<std::size_t>(i));
    }
    return result;
  }

  struct linear_regression_core_wrappers
  {
    typedef math::linear_regression_core<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      using boost::python::arg;
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
      using boost::python::arg;
      class_<w_t, bases<math::linear_regression_core<> > >(
            "linear_regression", no_init)
        .def(init<
          af::const_ref<float_t> const&,
          af::const_ref<float_t> const&,
          float_t const&>((
            arg("x"),
            arg("y"),
            arg("epsilon")=1e-15)))
        .def(init<
          af::const_ref<float_t> const&,
          af::const_ref<float_t> const&,
          af::const_ref<float_t> const&,
          float_t const&>((
            arg("x"),
            arg("y"),
            arg("weights"),
            arg("epsilon")=1e-15)))
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
      using boost::python::arg;
      class_<w_t>("linear_correlation", no_init)
        .def(init<
          af::const_ref<float_t> const&,
          af::const_ref<float_t> const&,
          float_t const&>(
            (arg("x"),
             arg("y"),
             arg("epsilon")=1.e-15)))
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

  boost::python::tuple
  integer_offsets_vs_pointers(
    af::ref<double> const& data,
    af::ref<std::size_t> const& permutation,
    unsigned n_repeats,
    bool use_pointers,
    int use_iterators)
  {
    SCITBX_ASSERT(permutation.size() == data.size());
    std::vector<double*> data_pointers;
    data_pointers.reserve(data.size());
    for(unsigned i=0;i<data.size();i++) {
      data_pointers.push_back(&data[permutation[i]]);
    }
    const char* result_type = 0;
    double result_value = 0;
    if (!use_pointers) {
      if (use_iterators == 0) {
        result_type = "d[p[i]]";
        for(unsigned i_repeat=0;i_repeat<n_repeats;i_repeat++) {
          for(unsigned i=0;i<data.size();i++) {
            result_value += data[permutation[i]];
          }
          for(unsigned i=0;i<data.size();i++) {
            result_value -= data[permutation[i]];
          }
        }
      }
      else if (use_iterators == 1) {
        result_type = "d[*p++]";
        for(unsigned i_repeat=0;i_repeat<n_repeats;i_repeat++) {
          std::size_t* p_end = permutation.end();
          std::size_t* p = permutation.begin();
          while (p != p_end) {
            result_value += data[*p++];
          }
          p = permutation.begin();
          while (p != p_end) {
            result_value -= data[*p++];
          }
        }
      }
      else {
        throw std::runtime_error("use_iterators: value error");
      }
    }
    else {
      if (use_iterators == 0) {
        result_type = "*dp[i]";
        for(unsigned i_repeat=0;i_repeat<n_repeats;i_repeat++) {
          for(unsigned i=0;i<data.size();i++) {
            result_value += *data_pointers[i];
          }
          for(unsigned i=0;i<data.size();i++) {
            result_value -= *data_pointers[i];
          }
        }
      }
      else if (use_iterators == 1) {
        result_type = "**dpi++";
        for(unsigned i_repeat=0;i_repeat<n_repeats;i_repeat++) {
          std::vector<double*>::const_iterator dp_end = data_pointers.end();
          std::vector<double*>::const_iterator dp = data_pointers.begin();
          while (dp != dp_end) {
            result_value += **dp++;
          }
          dp = data_pointers.begin();
          while (dp != dp_end) {
            result_value -= **dp++;
          }
        }
      }
      else if (use_iterators == 2) {
        result_type = "**dpp++";
        for(unsigned i_repeat=0;i_repeat<n_repeats;i_repeat++) {
          double** dp = &*data_pointers.begin();
          double** dp_end = dp + data_pointers.size();
          while (dp != dp_end) {
            result_value += **dp++;
          }
          dp = &*data_pointers.begin();
          while (dp != dp_end) {
            result_value -= **dp++;
          }
        }
      }
      else {
        throw std::runtime_error("use_iterators: value error");
      }
    }
    return boost::python::make_tuple(result_type, result_value);
  }

} // namespace <anonymous>

  struct cost_of_m_handle_in_af_shared
  {
    af::shared<double> input, result;

    cost_of_m_handle_in_af_shared(af::shared<double> const &data)
      : input(data),
        result(data.size(), af::init_functor_null<double>())
    {}

    const char*
    operator()(
      unsigned n_repeats,
      unsigned test_id)
    {
      if (test_id == 0) {
        for(int n=0; n < n_repeats; ++n) {
          for(std::size_t i=1; i < input.size(); ++i) {
            result[i] = input[i] - input[i-1];
          }
        }
        return "size+begin inside  loop";
      }
      if (test_id == 1) {
        for(int n=0; n < n_repeats; ++n) {
          double *r = result.begin();
          for(std::size_t i=1; i < input.size(); ++i) {
            r[i] = input[i] - input[i-1];
          }
        }
        return "     begin outside loop";
      }
      for(int n=0; n < n_repeats; ++n) {
        double *r = result.begin();
        std::size_t s = input.size();
        for(std::size_t i=1; i < s; ++i) {
          r[i] = input[i] - input[i-1];
        }
      }
      return "size+begin outside loop";
    }
  };

  // Testing argument passing from Python to C++
  // from flex.double to various scitbx::af array types
  struct flex_argument_passing
  {
    double x[3];

    flex_argument_passing() { x[0] = 1.5; x[1] = 2.5; x[2] = 3.5;}

    // This template pattern ensures that this function is as easy to use
    // from C++ (instantiating it with af::shared<double> for example)
    // as it is easy to wrap in a Boost.Python binding
    // which will transparently pass a flex.double through
    // (instantiating it with af::flex_1d<double>, c.f. below)
    template<template<class> class SharedArray1D>
    void easy_versa_flex_grid_as_reference(SharedArray1D<double> a) {
      a.extend(x, x+3);
      check(a);
      a.push_back(4.5);
      a.insert(&a[1], 0.5);
      SCITBX_ASSERT(a.begin() == &a[0]);
      SCITBX_ASSERT(a.end() == &a[5]);
      SCITBX_ASSERT(a.ref().size() == 5);
      SCITBX_ASSERT(a.ref()[2] == 2.5);
    }

    template<class A>
    void check(A &a) {
      SCITBX_ASSERT(a.size() == 3);
      SCITBX_ASSERT(a[0] == x[0]);
      SCITBX_ASSERT(a[1] == x[1]);
      SCITBX_ASSERT(a[2] == x[2]);
    }

    void shared_as_reference_fails(shared<double> &a) {
      a.extend(x, x+3);
      check(a);
    }

    void shared_as_value_fails(shared<double> a) {
      a.extend(x, x+3);
      check(a);
    }

    void versa_flex_grid_as_reference_succeeds(versa<double, flex_grid<> > &a){
      versa<double, flex_grid<> >::base_array_type b = a.as_base_array();
      b.extend(x, x+3);
      a.resize(flex_grid<>(b.size()));
      check(a);
    }

    void versa_flex_grid_as_value_fails(versa<double, flex_grid<> > a){
      versa<double, flex_grid<> >::base_array_type b = a.as_base_array();
      b.extend(x, x+3);
      a.resize(flex_grid<>(b.size()));
      check(a);
    }
  };

  int tst_c_grid_flex_conversion(
    const_ref<int, c_grid_periodic<3> > const &a,
    int i, int j, int k)
  {
    return a(i,j,k);
  }

namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;

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
    wrap_flex_vec3_double();
    wrap_flex_vec2_double();
    wrap_flex_sym_mat3_double();
    wrap_flex_tiny_size_t_2();

    default_c_grid_flex_conversions<int>();
    default_c_grid_flex_conversions<long>();
    default_c_grid_flex_conversions<float>();
    default_c_grid_flex_conversions<double>();
    default_c_grid_flex_conversions<std::complex<double> >();
    def("tst_c_grid_flex_conversion", tst_c_grid_flex_conversion);

    ref_owning_shared_conversions<double>();
    ref_owning_shared_conversions< std::complex<double> >();

    {
      namespace oc = boost_adaptbx::optional_conversions;
      oc::to_and_from_python<boost::optional<vec3<double> > >();
      oc::to_and_from_python<boost::optional<af::shared<double> > >();
      oc::to_and_from_python<
       boost_adaptbx::optional_container<af::small<int, 10> > >();
      oc::to_and_from_python<
       boost_adaptbx::optional_container<af::small<double, 6> > >();
    }

    wrap_flex_random();
    wrap_flex_sort();
    wrap_flex_histogram();
    wrap_flex_mean_and_variance();
    wrap_flex_median_statistics();
    wrap_flex_linear_interpolation();

    wrap_loops();
    wrap_empty_container_sizes();

    def("slice_indices", slice_indices, (
      arg("array_size"), arg("python_slice")));

    linear_regression_core_wrappers::wrap();
    linear_regression_wrappers::wrap();
    linear_correlation_wrappers::wrap();

    def("integer_offsets_vs_pointers", integer_offsets_vs_pointers);

    {
      typedef cost_of_m_handle_in_af_shared wt;
      class_<wt>("cost_of_m_handle_in_af_shared", no_init)
        .def(init<af::shared<double> const &>((arg("data"))))
        .add_property("result", &wt::result)
        .def("__call__", &wt::operator(),
             (arg("n_repeats"), arg("test_id")))
        ;
    }
    {
      typedef flex_argument_passing wt;
      typedef void (wt::*easy_t)(flex_1d<double>);
      class_<wt>("flex_argument_passing")
        .def("easy_versa_flex_grid_as_reference",
             (easy_t)&wt::easy_versa_flex_grid_as_reference)
        .def("shared_as_reference_fails", &wt::shared_as_reference_fails)
        .def("shared_as_value_fails", &wt::shared_as_value_fails)
        .def("versa_flex_grid_as_reference_succeeds",
             &wt::versa_flex_grid_as_reference_succeeds)
        .def("versa_flex_grid_as_value_fails",
             &wt::versa_flex_grid_as_value_fails)
        ;
    }
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_array_family_flex_ext)
{
  scitbx::af::boost_python::init_module();
}
